# R/publication_extensions.R
# Publication-oriented extensions for LogicComm.

#' Extract Spatial Coordinates
#'
#' Extracts a two-dimensional coordinate matrix from a matrix/data frame or,
#' when available, from a Seurat object. The returned object has rows as cells or
#' spots and columns named x and y.
#'
#' @param object Coordinate matrix/data frame or Seurat object.
#' @param coord_cols Optional coordinate columns to use.
#' @param reduction Optional Seurat reduction name used as a coordinate source.
#' @param image Optional Seurat image name passed to Seurat::GetTissueCoordinates().
#' @return Numeric matrix with columns x and y.
#' @export
extract_spatial_coordinates <- function(object,
                                        coord_cols = NULL,
                                        reduction = NULL,
                                        image = NULL) {
  if (inherits(object, "Seurat")) {
    if (!is.null(reduction)) {
      emb <- tryCatch(object@reductions[[reduction]]@cell.embeddings,
                      error = function(e) NULL)
      if (is.null(emb)) stop("Could not extract Seurat reduction: ", reduction)
      return(.prepare_xy_coordinates(emb, coord_cols = coord_cols))
    }

    coords <- NULL
    if (requireNamespace("Seurat", quietly = TRUE)) {
      coords <- tryCatch(
        Seurat::GetTissueCoordinates(object, image = image),
        error = function(e) NULL
      )
    }
    if (!is.null(coords) && nrow(coords) > 0) {
      return(.prepare_xy_coordinates(coords, coord_cols = coord_cols))
    }

    for (red in c("spatial", "umap", "tsne", "pca")) {
      emb <- tryCatch(object@reductions[[red]]@cell.embeddings,
                      error = function(e) NULL)
      if (!is.null(emb) && ncol(emb) >= 2) {
        return(.prepare_xy_coordinates(emb, coord_cols = coord_cols))
      }
    }
    stop("No spatial coordinates or 2D reduction were found in the Seurat object.")
  }

  .prepare_xy_coordinates(object, coord_cols = coord_cols)
}

#' Build a Spatial Neighborhood Graph
#'
#' Builds a sparse cell/spot adjacency graph from spatial coordinates. The graph
#' can be passed directly to IdentifyLogicConsensus() or
#' summarize_celltype_communication().
#'
#' @param coords Coordinate matrix/data frame or Seurat object.
#' @param coord_cols Optional coordinate columns.
#' @param mode "knn" for k-nearest neighbors or "radius" for a fixed-radius graph.
#' @param k Number of neighbors in kNN mode.
#' @param radius Radius in coordinate units for radius mode.
#' @param directed Whether to keep only directed i -> neighbor edges. If FALSE,
#'   edges are symmetrized by maximum weight.
#' @param distance_weight "binary", "inverse", or "gaussian".
#' @param sigma Gaussian kernel bandwidth. If NULL, the median selected edge
#'   distance is used.
#' @param verbose Print progress messages.
#' @return Sparse adjacency matrix with cells/spots as rows and columns.
#' @export
build_spatial_graph <- function(coords,
                                coord_cols = NULL,
                                mode = c("knn", "radius"),
                                k = 6,
                                radius = NULL,
                                directed = FALSE,
                                distance_weight = c("binary", "inverse", "gaussian"),
                                sigma = NULL,
                                verbose = TRUE) {
  mode <- match.arg(mode)
  distance_weight <- match.arg(distance_weight)
  xy <- extract_spatial_coordinates(coords, coord_cols = coord_cols)
  n <- nrow(xy)
  if (n < 2) stop("At least two coordinates are required to build a spatial graph.")
  if (mode == "knn") {
    if (!is.numeric(k) || length(k) != 1 || k < 1) stop("k must be a positive integer.")
    k <- as.integer(min(k, n - 1L))
  } else {
    if (is.null(radius) || !is.numeric(radius) || length(radius) != 1 || radius <= 0) {
      stop("radius must be a positive scalar in radius mode.")
    }
  }

  if (isTRUE(verbose)) {
    msg <- if (mode == "knn") sprintf("[Spatial] Building %d-NN graph for %d coordinates.", k, n) else
      sprintf("[Spatial] Building radius graph (r = %.3g) for %d coordinates.", radius, n)
    message(msg)
  }

  from <- integer(0)
  to <- integer(0)
  dist_vec <- numeric(0)
  for (i in seq_len(n)) {
    d <- sqrt(rowSums(sweep(xy, 2, xy[i, ], "-")^2))
    d[i] <- Inf
    if (mode == "knn") {
      ord <- order(d, na.last = NA)
      nn <- ord[seq_len(min(k, length(ord)))]
      nn <- nn[is.finite(d[nn])]
    } else {
      nn <- which(is.finite(d) & d <= radius)
    }
    if (length(nn) > 0) {
      from <- c(from, rep.int(i, length(nn)))
      to <- c(to, nn)
      dist_vec <- c(dist_vec, d[nn])
    }
  }
  if (length(from) == 0) stop("No spatial edges were created. Increase k or radius.")
  if (!isTRUE(directed)) {
    from0 <- from
    to0 <- to
    dist0 <- dist_vec
    from <- c(from0, to0)
    to <- c(to0, from0)
    dist_vec <- c(dist0, dist0)
  }
  w <- .distance_to_edge_weight(dist_vec, distance_weight = distance_weight, sigma = sigma)
  df <- data.frame(i = from, j = to, x = w, stringsAsFactors = FALSE)
  df <- df[df$i != df$j & is.finite(df$x) & df$x > 0, , drop = FALSE]
  if (nrow(df) == 0) stop("No positive non-self spatial edges remain.")
  agg <- stats::aggregate(x ~ i + j, df, max)
  out <- Matrix::sparseMatrix(
    i = agg$i, j = agg$j, x = agg$x,
    dims = c(n, n), dimnames = list(rownames(xy), rownames(xy))
  )
  attr(out, "graph_type") <- "spatial"
  attr(out, "spatial_mode") <- mode
  attr(out, "distance_weight") <- distance_weight
  attr(out, "directed") <- isTRUE(directed)
  out
}

#' Summarize Spatial Logic Communication
#'
#' Convenience wrapper that builds a spatial graph and then calls
#' summarize_celltype_communication().
#'
#' @param reo_mat Binary REO matrix or LogicCommREOResult.
#' @param coords Spatial coordinates or Seurat object.
#' @param cell_labels Named cell-type labels.
#' @param lr_db LR database.
#' @param coord_cols Optional coordinate columns.
#' @param spatial_mode \code{"knn"} or \code{"radius"} spatial graph construction.
#' @param k Number of neighbors in kNN mode.
#' @param radius Radius threshold in radius mode.
#' @param directed Whether to keep directed spatial edges.
#' @param distance_weight \code{"binary"}, \code{"inverse"}, or \code{"gaussian"}.
#' @param sigma Gaussian distance kernel bandwidth.
#' @param graph_symmetrize How to symmetrize the spatial graph before scoring.
#' @param edge_weight_mode Optional edge weighting mode passed to communication scoring.
#' @param verbose Print progress messages.
#' @param ... Additional arguments passed to summarize_celltype_communication().
#' @return LogicCommCellTypeComm object with spatial_graph and spatial_coords.
#' @export
summarize_spatial_communication <- function(reo_mat,
                                            coords,
                                            cell_labels,
                                            lr_db = lr_pairs_human,
                                            coord_cols = NULL,
                                            spatial_mode = c("knn", "radius"),
                                            k = 6,
                                            radius = NULL,
                                            directed = FALSE,
                                            distance_weight = c("binary", "inverse", "gaussian"),
                                            sigma = NULL,
                                            graph_symmetrize = c("none", "or", "max"),
                                            edge_weight_mode = NULL,
                                            verbose = TRUE,
                                            ...) {
  if (is.list(reo_mat) && !is.null(reo_mat$logic)) reo_logic <- reo_mat$logic else reo_logic <- reo_mat
  xy <- extract_spatial_coordinates(coords, coord_cols = coord_cols)
  if (is.null(rownames(xy))) stop("Spatial coordinates must have cell/spot rownames.")
  if (is.null(colnames(reo_logic))) stop("reo_mat must have cell colnames.")
  if (identical(rownames(xy), as.character(seq_len(nrow(xy)))) && nrow(xy) == ncol(reo_logic)) {
    rownames(xy) <- colnames(reo_logic)
  }
  missing <- setdiff(colnames(reo_logic), rownames(xy))
  if (length(missing) > 0) stop("Spatial coordinates are missing ", length(missing), " REO cells/spots. Example: ", missing[1])
  xy <- xy[colnames(reo_logic), , drop = FALSE]
  distance_weight <- match.arg(distance_weight)
  spatial_mode <- match.arg(spatial_mode)
  graph_symmetrize <- match.arg(graph_symmetrize)
  if (is.null(edge_weight_mode)) edge_weight_mode <- if (distance_weight == "binary") "binary" else "weighted"
  graph <- build_spatial_graph(
    xy, mode = spatial_mode, k = k, radius = radius, directed = directed,
    distance_weight = distance_weight, sigma = sigma, verbose = verbose
  )
  ct <- summarize_celltype_communication(
    reo_mat = reo_logic, cell_labels = cell_labels, knn_mat = graph, lr_db = lr_db,
    mode = "neighborhood", graph_symmetrize = graph_symmetrize,
    edge_weight_mode = edge_weight_mode, verbose = verbose, ...
  )
  ct$spatial_graph <- graph
  ct$spatial_coords <- xy
  ct$params$spatial_mode <- spatial_mode
  ct$params$spatial_distance_weight <- distance_weight
  class(ct) <- "LogicCommCellTypeComm"
  ct
}

#' Plot Spatial Logic States for One L-R Pair
#'
#' @param reo_mat Binary REO matrix or LogicCommREOResult.
#' @param coords Spatial coordinates or Seurat object.
#' @param lr_pair L-R pair to display.
#' @param lr_db LR database.
#' @param coord_cols Optional coordinate columns.
#' @param cell_labels Optional named cell-type labels.
#' @param sender_type Optional sender cell type to highlight.
#' @param receiver_type Optional receiver cell type to highlight.
#' @param pt_size Point size.
#' @param title Optional title.
#' @return ggplot2 object.
#' @export
plot_spatial_logic <- function(reo_mat,
                               coords,
                               lr_pair,
                               lr_db = lr_pairs_human,
                               coord_cols = NULL,
                               cell_labels = NULL,
                               sender_type = NULL,
                               receiver_type = NULL,
                               pt_size = 1.2,
                               title = NULL) {
  if (is.list(reo_mat) && !is.null(reo_mat$logic)) reo_mat <- reo_mat$logic
  xy <- extract_spatial_coordinates(coords, coord_cols = coord_cols)
  if (identical(rownames(xy), as.character(seq_len(nrow(xy)))) && nrow(xy) == ncol(reo_mat)) {
    rownames(xy) <- colnames(reo_mat)
  }
  common <- intersect(colnames(reo_mat), rownames(xy))
  if (length(common) == 0) stop("No overlapping cells between reo_mat and coords.")
  xy <- xy[common, , drop = FALSE]
  reo_sub <- reo_mat[, common, drop = FALSE]
  lr_row <- lr_db[lr_db$lr_pair == lr_pair, , drop = FALSE]
  if (nrow(lr_row) == 0) stop("lr_pair not found in lr_db: ", lr_pair)
  lig <- .resolve_complex_logic(.lr_list_entry(lr_row, "ligand_genes", 1), reo_sub)
  rec <- .resolve_complex_logic(.lr_list_entry(lr_row, "receptor_genes", 1), reo_sub)
  status <- rep("Inactive", length(common))
  status[lig] <- "Ligand active"
  status[rec] <- "Receptor active"
  status[lig & rec] <- "Both active"
  if (!is.null(cell_labels)) {
    if (is.null(names(cell_labels))) {
      if (length(cell_labels) != length(colnames(reo_mat))) stop("Unnamed cell_labels must match ncol(reo_mat).")
      names(cell_labels) <- colnames(reo_mat)
    }
    labs <- as.character(cell_labels[common])
    if (!is.null(sender_type)) status[labs %in% sender_type & lig] <- "Selected sender ligand active"
    if (!is.null(receiver_type)) status[labs %in% receiver_type & rec] <- "Selected receiver receptor active"
    if (!is.null(sender_type) && !is.null(receiver_type)) {
      status[labs %in% sender_type & labs %in% receiver_type & lig & rec] <- "Selected sender/receiver both active"
    }
  }
  df <- data.frame(x = xy[, 1], y = xy[, 2], status = status, stringsAsFactors = FALSE)
  if (is.null(title)) title <- paste("Spatial L-R logic:", lr_pair)
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = status)) +
    ggplot2::geom_point(size = pt_size, alpha = 0.85) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, x = "spatial x", y = "spatial y", color = "Logic state") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

#' Summarize Communication Dynamics Along Pseudotime or Time Bins
#'
#' @param reo_mat Binary REO matrix or LogicCommREOResult.
#' @param pseudotime Named numeric/factor vector.
#' @param cell_labels Named cell-type labels.
#' @param lr_db LR database.
#' @param knn_mat Optional graph. If supplied, each bin is scored on the induced
#'   subgraph.
#' @param n_bins Number of bins for numeric pseudotime.
#' @param bin_method "quantile", "equal_width", or "categorical".
#' @param min_cells_per_bin Minimum cells required in a bin.
#' @param verbose Print progress messages.
#' @param ... Passed to summarize_celltype_communication().
#' @return LogicCommDynamics object.
#' @export
summarize_communication_dynamics <- function(reo_mat,
                                             pseudotime,
                                             cell_labels,
                                             lr_db = lr_pairs_human,
                                             knn_mat = NULL,
                                             n_bins = 5,
                                             bin_method = c("quantile", "equal_width", "categorical"),
                                             min_cells_per_bin = 20,
                                             verbose = TRUE,
                                             ...) {
  if (is.list(reo_mat) && !is.null(reo_mat$logic)) reo_mat <- reo_mat$logic
  stopifnot(inherits(reo_mat, c("dgCMatrix", "sparseMatrix", "matrix")))
  if (is.null(colnames(reo_mat))) stop("reo_mat must have cell colnames.")
  cells <- colnames(reo_mat)
  pseudotime <- .align_named_vector(pseudotime, cells, "pseudotime")
  cell_labels <- .align_named_vector(cell_labels, cells, "cell_labels")
  bin_method <- match.arg(bin_method)
  bins <- .make_time_bins(pseudotime, n_bins = n_bins, method = bin_method)
  valid_bins <- names(which(table(bins) >= min_cells_per_bin))
  if (length(valid_bins) == 0) stop("No pseudotime/time bin has at least min_cells_per_bin cells.")
  if (!is.null(knn_mat)) knn_mat <- .validate_and_align_knn(knn_mat, cells)

  bin_results <- vector("list", length(valid_bins))
  names(bin_results) <- valid_bins
  pair_list <- role_list <- lr_list <- vector("list", length(valid_bins))
  names(pair_list) <- names(role_list) <- names(lr_list) <- valid_bins

  for (b in valid_bins) {
    keep <- names(bins)[bins == b]
    if (isTRUE(verbose)) message(sprintf("[Dynamics] Scoring bin %s (%d cells).", b, length(keep)))
    reo_b <- reo_mat[, keep, drop = FALSE]
    labels_b <- cell_labels[keep]
    knn_b <- if (!is.null(knn_mat)) knn_mat[keep, keep, drop = FALSE] else NULL
    ct_b <- summarize_celltype_communication(
      reo_b, cell_labels = labels_b, knn_mat = knn_b, lr_db = lr_db,
      verbose = FALSE, ...
    )
    bin_results[[b]] <- ct_b
    pair_list[[b]] <- cbind(bin = b, bin_order = match(b, valid_bins), n_bin_cells = length(keep), ct_b$pair_summary)
    role_list[[b]] <- cbind(bin = b, bin_order = match(b, valid_bins), n_bin_cells = length(keep), ct_b$role_summary)
    lr_list[[b]] <- cbind(bin = b, bin_order = match(b, valid_bins), n_bin_cells = length(keep), ct_b$lr_table)
  }

  out <- list(
    bin_results = bin_results,
    pair_trajectory = do.call(rbind, pair_list),
    role_trajectory = do.call(rbind, role_list),
    lr_trajectory = do.call(rbind, lr_list),
    bins = bins,
    params = list(n_bins = n_bins, bin_method = bin_method,
                  min_cells_per_bin = min_cells_per_bin)
  )
  class(out) <- "LogicCommDynamics"
  out
}

#' Plot Communication Dynamics
#'
#' @param dynamics Output from summarize_communication_dynamics().
#' @param level "pair", "role", or "lr".
#' @param metric Metric to plot.
#' @param sender Optional sender filter.
#' @param receiver Optional receiver filter.
#' @param lr_pair Optional LR-pair filter for level = "lr".
#' @param cell_type Optional cell-type filter for level = "role".
#' @param top_n If no specific feature is selected, show top_n features by total metric.
#' @param title Optional plot title.
#' @return ggplot2 object.
#' @export
plot_communication_dynamics <- function(dynamics,
                                        level = c("pair", "role", "lr"),
                                        metric = NULL,
                                        sender = NULL,
                                        receiver = NULL,
                                        lr_pair = NULL,
                                        cell_type = NULL,
                                        top_n = 10,
                                        title = NULL) {
  stopifnot(inherits(dynamics, "LogicCommDynamics"))
  level <- match.arg(level)
  if (level == "pair") {
    df <- dynamics$pair_trajectory
    metric <- metric %||% "sum_lcs"
    if (!is.null(sender)) df <- df[df$sender_type %in% sender, , drop = FALSE]
    if (!is.null(receiver)) df <- df[df$receiver_type %in% receiver, , drop = FALSE]
    if (!metric %in% names(df)) stop("metric not found: ", metric)
    df$feature <- paste(df$sender_type, df$receiver_type, sep = " -> ")
  } else if (level == "role") {
    df <- dynamics$role_trajectory
    metric <- metric %||% "hub_score"
    if (!is.null(cell_type)) df <- df[df$cell_type %in% cell_type, , drop = FALSE]
    if (!metric %in% names(df)) stop("metric not found: ", metric)
    df$feature <- df$cell_type
  } else {
    df <- dynamics$lr_trajectory
    metric <- metric %||% "lcs"
    if (!is.null(sender)) df <- df[df$sender_type %in% sender, , drop = FALSE]
    if (!is.null(receiver)) df <- df[df$receiver_type %in% receiver, , drop = FALSE]
    if (!is.null(lr_pair)) df <- df[df$lr_pair %in% lr_pair, , drop = FALSE]
    if (!metric %in% names(df)) stop("metric not found: ", metric)
    df$feature <- paste(df$sender_type, df$receiver_type, df$lr_pair, sep = " -> ")
  }
  df <- df[is.finite(df[[metric]]), , drop = FALSE]
  if (nrow(df) == 0) stop("No rows remain for the selected dynamic plot.")
  if (is.null(sender) && is.null(receiver) && is.null(lr_pair) && is.null(cell_type)) {
    totals <- tapply(df[[metric]], df$feature, sum, na.rm = TRUE)
    keep <- names(sort(totals, decreasing = TRUE))[seq_len(min(top_n, length(totals)))]
    df <- df[df$feature %in% keep, , drop = FALSE]
  }
  df$value <- df[[metric]]
  bin_map <- unique(df[order(df$bin_order), c("bin_order", "bin")])
  if (is.null(title)) title <- paste("LogicComm dynamics:", level, metric)
  ggplot2::ggplot(df, ggplot2::aes(x = bin_order, y = value, color = feature, group = feature)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(breaks = bin_map$bin_order, labels = as.character(bin_map$bin)) +
    ggplot2::labs(title = title, x = "Pseudotime/time bin", y = metric, color = "Feature") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}

#' Fit Sample-Level Binomial/Quasibinomial Models for Cell-Type Communication
#'
#' Fits a per-feature generalized linear model using active edge counts as
#' successes and available edge opportunities as trials. The statistical unit is
#' the biological sample.
#'
#' @param sample_ct_list Named list of LogicCommCellTypeComm objects.
#' @param sample_metadata Data frame with one row per sample.
#' @param design Model formula, e.g. ~ group + batch.
#' @param coef Coefficient name or index to report. Defaults to the last
#'   non-intercept coefficient.
#' @param level "celltype_lr" or "celltype_pair".
#' @param min_samples Minimum samples with usable counts for a feature.
#' @param min_total_trials Minimum total trials across usable samples.
#' @param fdr_method Multiple-testing method.
#' @param verbose Print progress messages.
#' @return Data frame of model statistics.
#' @export
fit_celltype_comm_glm <- function(sample_ct_list,
                                  sample_metadata,
                                  design = ~ group,
                                  coef = NULL,
                                  level = c("celltype_lr", "celltype_pair"),
                                  min_samples = 3,
                                  min_total_trials = 20,
                                  fdr_method = "BH",
                                  verbose = TRUE) {
  level <- match.arg(level)
  if (is.null(names(sample_ct_list)) || any(!nzchar(names(sample_ct_list)))) {
    stop("sample_ct_list must be a named list.")
  }
  if (is.null(rownames(sample_metadata))) {
    if (nrow(sample_metadata) != length(sample_ct_list)) stop("sample_metadata needs rownames or one row per sample.")
    rownames(sample_metadata) <- names(sample_ct_list)
  }
  sample_names <- intersect(names(sample_ct_list), rownames(sample_metadata))
  if (length(sample_names) < min_samples) stop("Fewer than min_samples have both communication objects and metadata.")
  sample_ct_list <- sample_ct_list[sample_names]
  sample_metadata <- sample_metadata[sample_names, , drop = FALSE]

  sample_tables <- lapply(sample_ct_list, .extract_comm_count_table, level = level)
  all_features <- Reduce(union, lapply(sample_tables, function(x) x$feature))
  X <- stats::model.matrix(design, data = sample_metadata)
  coef_names <- colnames(X)
  if (is.null(coef)) {
    candidate <- setdiff(coef_names, "(Intercept)")
    if (length(candidate) == 0) stop("design must contain at least one non-intercept coefficient.")
    coef <- candidate[length(candidate)]
  } else if (is.numeric(coef)) {
    coef <- coef_names[coef]
  }
  if (!coef %in% coef_names) stop("coef not found in model matrix: ", coef)

  if (isTRUE(verbose)) {
    message(sprintf("[GLM] Fitting %d %s features across %d samples; coefficient = %s.",
                    length(all_features), level, length(sample_names), coef))
  }

  out_list <- vector("list", length(all_features))
  for (f in seq_along(all_features)) {
    feature <- all_features[f]
    success <- total <- lcs <- rep(NA_real_, length(sample_names))
    names(success) <- names(total) <- names(lcs) <- sample_names
    for (sname in sample_names) {
      tab <- sample_tables[[sname]]
      ii <- match(feature, tab$feature)
      if (!is.na(ii)) {
        success[sname] <- tab$success[ii]
        total[sname] <- tab$total[ii]
        lcs[sname] <- tab$lcs[ii]
      }
    }
    keep <- is.finite(success) & is.finite(total) & total > 0 & success >= 0 & success <= total
    if (sum(keep) < min_samples || sum(total[keep], na.rm = TRUE) < min_total_trials) next
    dat <- sample_metadata[keep, , drop = FALSE]
    dat$success <- success[keep]
    dat$failure <- pmax(total[keep] - success[keep], 0)
    rhs <- paste(deparse(design), collapse = "")
    form <- stats::as.formula(paste("cbind(success, failure)", rhs))
    fit <- tryCatch(suppressWarnings(stats::glm(form, data = dat, family = stats::quasibinomial())),
                    error = function(e) NULL)
    if (is.null(fit)) next
    cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (is.null(cf) || !coef %in% rownames(cf)) next
    row <- cf[coef, ]
    parsed <- .parse_celltype_feature_names(feature)
    out_list[[f]] <- data.frame(
      feature = feature,
      sender_type = parsed$sender_type,
      receiver_type = parsed$receiver_type,
      lr_pair = parsed$feature,
      estimate = as.numeric(row[["Estimate"]]),
      std_error = as.numeric(row[["Std. Error"]]),
      statistic = as.numeric(row[[3]]),
      p_value = as.numeric(row[[4]]),
      odds_ratio = exp(as.numeric(row[["Estimate"]])),
      n_samples = sum(keep),
      total_success = sum(success[keep], na.rm = TRUE),
      total_trials = sum(total[keep], na.rm = TRUE),
      mean_lcs = mean(lcs[keep], na.rm = TRUE),
      coef = coef,
      level = level,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, out_list)
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame(feature = character(0), estimate = numeric(0), p_value = numeric(0), fdr = numeric(0)))
  }
  out$fdr <- stats::p.adjust(out$p_value, method = fdr_method)
  out <- out[order(out$fdr, out$p_value), , drop = FALSE]
  class(out) <- c("LogicCommGLMResult", "data.frame")
  attr(out, "design") <- design
  attr(out, "coef") <- coef
  out
}

#' Plot a GLM Differential Communication Volcano
#'
#' @param glm_result Output from fit_celltype_comm_glm().
#' @param top_n_label Number of top features to label.
#' @param title Optional plot title.
#' @return ggplot2 object.
#' @export
plot_celltype_glm_volcano <- function(glm_result,
                                      top_n_label = 15,
                                      title = NULL) {
  stopifnot(is.data.frame(glm_result))
  df <- as.data.frame(glm_result)
  if (!all(c("estimate", "fdr") %in% names(df))) stop("glm_result must contain estimate and fdr columns.")
  df$neglog10 <- -log10(pmax(df$fdr, .Machine$double.xmin))
  if (!"sender_type" %in% names(df)) df$sender_type <- NA_character_
  df$label <- .compact_feature_label(df$feature)
  ord <- order(df$fdr, -abs(df$estimate), na.last = NA)
  df$to_label <- FALSE
  if (length(ord) > 0 && top_n_label > 0) df$to_label[ord[seq_len(min(top_n_label, length(ord)))]] <- TRUE
  if (is.null(title)) title <- "Sample-level GLM differential communication"
  ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = neglog10)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey65") +
    ggplot2::geom_point(ggplot2::aes(color = sender_type), alpha = 0.8) +
    ggrepel::geom_text_repel(data = df[df$to_label, , drop = FALSE],
                             ggplot2::aes(label = label), size = 3, max.overlaps = 12,
                             min.segment.length = 0) +
    ggplot2::labs(title = title, x = "GLM coefficient (log odds)", y = "-log10(FDR)", color = "Sender") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

#' Summarize Publication-Level Communication Findings
#'
#' @param ct_comm Output from summarize_celltype_communication().
#' @param top_n Number of rows per summary table.
#' @return LogicCommFindings list with top events, roles, pathways, and QC.
#' @export
summarize_communication_findings <- function(ct_comm, top_n = 10) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  pair <- ct_comm$pair_summary
  if (!"communication_support_label" %in% names(pair)) {
    pair$communication_support_label <- ifelse(pair$n_active_lr > 0 & pair$n_local_active == 0 & pair$n_distal_candidate > 0,
                                                "global_only_candidate",
                                                ifelse(pair$n_distal_candidate > 0, "mixed_local_global", "local_graph_supported"))
  }
  if (!"local_support_fraction_active" %in% names(pair)) {
    pair$local_support_fraction_active <- ifelse(pair$n_active_lr > 0, pair$n_local_active / pair$n_active_lr, NA_real_)
  }
  pair <- pair[is.finite(pair$sum_lcs) & pair$n_active_lr > 0, , drop = FALSE]
  support_rank <- match(pair$communication_support_label,
                        c("local_graph_supported", "mixed_local_global", "global_only_candidate", "inactive"))
  pair <- pair[order(support_rank, -pair$sum_lcs, -pair$n_active_lr, na.last = TRUE), , drop = FALSE]
  global_only_pair <- pair[pair$communication_support_label == "global_only_candidate", , drop = FALSE]
  global_only_pair <- global_only_pair[order(-global_only_pair$sum_lcs, -global_only_pair$n_active_lr), , drop = FALSE]
  primary_pair <- pair[pair$communication_support_label != "global_only_candidate", , drop = FALSE]
  if (nrow(primary_pair) == 0) primary_pair <- pair
  lr <- ct_comm$lr_table[ct_comm$lr_table$active %in% TRUE & !is.na(ct_comm$lr_table$lcs), , drop = FALSE]
  lr <- lr[order(lr$lcs, decreasing = TRUE), , drop = FALSE]
  path <- ct_comm$pathway_summary[order(ct_comm$pathway_summary$sum_lcs, decreasing = TRUE), , drop = FALSE]
  roles <- ct_comm$role_summary[order(ct_comm$role_summary$hub_score, decreasing = TRUE), , drop = FALSE]
  qc <- data.frame(
    n_cells = ct_comm$params$n_cells %||% length(ct_comm$cell_labels),
    n_cell_types = ct_comm$params$n_cell_types %||% nrow(ct_comm$role_summary),
    n_lr_pairs = ct_comm$params$n_lr_pairs %||% length(unique(ct_comm$lr_table$lr_pair)),
    n_active_celltype_lr_events = sum(ct_comm$lr_table$active, na.rm = TRUE),
    n_global_only_celltype_pairs = if ("communication_support_label" %in% names(ct_comm$pair_summary)) sum(ct_comm$pair_summary$communication_support_label == "global_only_candidate", na.rm = TRUE) else NA_integer_,
    median_edges_per_celltype_pair = stats::median(ct_comm$pair_summary$n_edges, na.rm = TRUE),
    n_low_communication_celltypes = sum(ct_comm$role_summary$low_communication, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  guide <- data.frame(
    biological_question = c(
      "Which locally supported or mixed cell-type pairs communicate most strongly?",
      "Which distal/global candidates lack local active edge support?",
      "Which L-R axes drive a selected sender-to-receiver interaction?",
      "Which cell types behave as senders, receivers, mediators, or influencers?",
      "Are conclusions supported by enough edge opportunities?"
    ),
    recommended_table = c("top_celltype_pairs", "distal_candidate_pairs", "top_lr_events", "role_summary", "qc + pair_summary"),
    recommended_metric = c("sum_lcs, n_active_lr, local_support_fraction_active, and communication_support_label", "sum_lcs, n_edges, local_support_fraction_active, and permutation diagnostics", "lcs, n_active_edges, and communication_range", "role scores, role_separation_label, communication_evidence_label, and role_reliability_label", "n_edges, communication_support_label, and bootstrap/permutation diagnostics"),
    stringsAsFactors = FALSE
  )
  out <- list(
    top_celltype_pairs = utils::head(primary_pair, top_n),
    distal_candidate_pairs = utils::head(global_only_pair, top_n),
    top_lr_events = utils::head(lr, top_n),
    top_pathways = utils::head(path, top_n),
    role_summary = utils::head(roles, top_n),
    qc = qc,
    interpretation_guide = guide
  )
  class(out) <- "LogicCommFindings"
  out
}

#' Write a Markdown Communication Report
#'
#' @param ct_comm Output from summarize_celltype_communication().
#' @param file Output markdown file.
#' @param title Report title.
#' @param top_n Number of rows per table.
#' @return Invisibly returns file path.
#' @export
write_communication_report <- function(ct_comm,
                                       file = "LogicComm_report.md",
                                       title = "LogicComm communication report",
                                       top_n = 10) {
  findings <- summarize_communication_findings(ct_comm, top_n = top_n)
  lines <- c(
    paste0("# ", title),
    "",
    "## Quality-control summary",
    .df_to_md(findings$qc),
    "",
    "## Top sender-to-receiver cell-type pairs with local or mixed support",
    .df_to_md(findings$top_celltype_pairs[, intersect(c("sender_type", "receiver_type", "communication_support_label", "local_support_fraction_active", "n_edges", "n_active_lr", "n_local_active", "n_distal_candidate", "sum_lcs", "mean_lcs_active", "top_pathway"), names(findings$top_celltype_pairs)), drop = FALSE]),
    "",
    "## Global-only distal candidate pairs",
    "These pairs have active global/distal evidence but no local active edge support in the supplied graph; treat them as candidates, not direct neighborhood interactions.",
    "",
    .df_to_md(findings$distal_candidate_pairs[, intersect(c("sender_type", "receiver_type", "communication_support_label", "local_support_fraction_active", "n_edges", "n_active_lr", "sum_lcs", "mean_lcs_active", "top_pathway"), names(findings$distal_candidate_pairs)), drop = FALSE]),
    "",
    "## Top cell-type-resolved L-R events",
    .df_to_md(findings$top_lr_events[, intersect(c("sender_type", "receiver_type", "lr_pair", "pathway", "n_edges", "n_active_edges", "lcs", "ligand_active_frac_sender", "receptor_active_frac_receiver"), names(findings$top_lr_events)), drop = FALSE]),
    "",
    "## Top pathways",
    .df_to_md(findings$top_pathways[, intersect(c("sender_type", "receiver_type", "pathway", "n_active_lr", "sum_lcs", "mean_lcs_active"), names(findings$top_pathways)), drop = FALSE]),
    "",
    "## Cell-type roles",
    .df_to_md(findings$role_summary[, intersect(c("cell_type", "dominant_role", "secondary_role", "hub_score", "sender_role_score", "receiver_role_score", "mediator_role_score", "influencer_role_score", "role_separation_label", "communication_evidence_label", "role_reliability_label", "dominant_role_strict"), names(findings$role_summary)), drop = FALSE]),
    "",
    "## Interpretation guide",
    .df_to_md(findings$interpretation_guide),
    "",
    "## Reporting cautions",
    "- LogicComm LCS is a rank-logic consensus score, not raw expression magnitude.",
    "- A KNN/SNN edge is a transcriptomic neighborhood opportunity, not necessarily physical contact unless a spatial graph is used.",
    "- Treat biological samples, not individual cells, as independent units for group comparison.",
    "- Inspect n_edges, active edge support, bootstrap intervals, and permutation nulls before emphasizing a rare cell-type pair.",
    "- Report global-only distal candidates separately from local or mixed graph-supported communication.",
    "- Separate broad or identity-associated axes from pair-specific mechanistic candidates using score_communication_specificity().",
    "- Interpret degenerate positive permutation nulls as structural-null flags, not ordinary Gaussian z-scores."
  )
  writeLines(lines, con = file)
  invisible(normalizePath(file, mustWork = FALSE))
}

#' Plot Communication Quality-Control Landscape
#'
#' @param ct_comm Output from summarize_celltype_communication().
#' @param title Optional title.
#' @return ggplot2 object.
#' @export
plot_communication_qc <- function(ct_comm, title = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  df <- ct_comm$pair_summary
  if (is.null(title)) title <- "Communication evidence and support QC"
  ggplot2::ggplot(df, ggplot2::aes(x = n_edges, y = sum_lcs)) +
    ggplot2::geom_point(ggplot2::aes(size = n_active_lr, color = mean_edge_support_fraction_active), alpha = 0.8) +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::scale_size_continuous(name = "Active LR events") +
    ggplot2::scale_color_gradient(low = "grey80", high = "firebrick", na.value = "grey80", name = "Mean support") +
    ggplot2::labs(title = title, x = "Cell-type-pair edge opportunities (log10)", y = "Summed active LCS") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

#' Plot Cell-Type Role Confidence
#'
#' @param ct_comm Output from summarize_celltype_communication().
#' @param title Optional title.
#' @return ggplot2 object.
#' @export
plot_role_confidence <- function(ct_comm, title = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  df <- ct_comm$role_summary
  df <- df[order(df$hub_score, decreasing = TRUE), , drop = FALSE]
  df$cell_type <- factor(df$cell_type, levels = rev(df$cell_type))
  if (is.null(title)) title <- "Cell-type communication role reliability"
  ggplot2::ggplot(df, ggplot2::aes(x = role_confidence, y = cell_type, fill = dominant_role)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::labs(title = title, x = "Role reliability score", y = "Cell type", fill = "Dominant role") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

# Print helpers ---------------------------------------------------------------

#' Print method for LogicCommDynamics
#' @param x LogicCommDynamics object.
#' @param ... Ignored.
#' @export
print.LogicCommDynamics <- function(x, ...) {
  cat(sprintf("LogicCommDynamics | %d bins | %d pair rows | %d role rows\n",
              length(x$bin_results), nrow(x$pair_trajectory), nrow(x$role_trajectory)))
  cat("Bins:\n")
  print(table(x$bins))
  invisible(x)
}

#' Print method for LogicCommFindings
#' @param x LogicCommFindings object.
#' @param ... Ignored.
#' @export
print.LogicCommFindings <- function(x, ...) {
  cat("LogicCommFindings\n")
  cat("QC:\n")
  print(x$qc)
  cat("\nTop cell-type pairs:\n")
  print(utils::head(x$top_celltype_pairs[, intersect(c("sender_type", "receiver_type", "sum_lcs", "n_active_lr", "top_pathway"), names(x$top_celltype_pairs)), drop = FALSE], 5))
  invisible(x)
}

# Internal helpers ------------------------------------------------------------

.prepare_xy_coordinates <- function(coords, coord_cols = NULL) {
  if (is.null(coords)) stop("coords must not be NULL.")
  coords <- as.data.frame(coords)
  if (nrow(coords) == 0) stop("coords has no rows.")
  if (is.null(coord_cols)) {
    candidates <- list(
      c("x", "y"), c("X", "Y"), c("spatial_1", "spatial_2"),
      c("imagecol", "imagerow"), c("pxl_col_in_fullres", "pxl_row_in_fullres"),
      c("col", "row")
    )
    for (cc in candidates) {
      if (all(cc %in% names(coords))) {
        coord_cols <- cc
        break
      }
    }
    if (is.null(coord_cols)) {
      numeric_cols <- names(coords)[vapply(coords, is.numeric, logical(1))]
      if (length(numeric_cols) < 2) stop("Could not identify two numeric coordinate columns.")
      coord_cols <- numeric_cols[1:2]
    }
  }
  if (length(coord_cols) != 2 || !all(coord_cols %in% names(coords))) {
    stop("coord_cols must name two columns in coords.")
  }
  xy <- as.matrix(coords[, coord_cols, drop = FALSE])
  storage.mode(xy) <- "double"
  if (any(!is.finite(xy))) stop("Coordinates contain non-finite values.")
  colnames(xy) <- c("x", "y")
  rownames(xy) <- rownames(coords)
  if (is.null(rownames(xy))) rownames(xy) <- paste0("Cell", seq_len(nrow(xy)))
  xy
}

.distance_to_edge_weight <- function(d, distance_weight = c("binary", "inverse", "gaussian"), sigma = NULL) {
  distance_weight <- match.arg(distance_weight)
  d <- as.numeric(d)
  if (distance_weight == "binary") return(rep(1, length(d)))
  if (distance_weight == "inverse") return(1 / pmax(d, .Machine$double.eps))
  if (is.null(sigma)) sigma <- stats::median(d[is.finite(d) & d > 0], na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 0) sigma <- 1
  exp(-(d^2) / (2 * sigma^2))
}

.align_named_vector <- function(x, cells, label) {
  if (is.null(names(x))) {
    if (length(x) != length(cells)) stop(label, " must be named or have length ncol(reo_mat).")
    names(x) <- cells
  }
  missing <- setdiff(cells, names(x))
  if (length(missing) > 0) stop(label, " is missing cells. Example: ", missing[1])
  out <- x[cells]
  names(out) <- cells
  out
}

.make_time_bins <- function(pseudotime, n_bins = 5, method = c("quantile", "equal_width", "categorical")) {
  method <- match.arg(method)
  if (method == "categorical" || !is.numeric(pseudotime)) {
    bins <- factor(as.character(pseudotime), levels = unique(as.character(pseudotime)))
    names(bins) <- names(pseudotime)
    return(bins)
  }
  pt <- as.numeric(pseudotime)
  if (method == "quantile") {
    probs <- seq(0, 1, length.out = n_bins + 1L)
    br <- unique(stats::quantile(pt, probs = probs, na.rm = TRUE, names = FALSE))
  } else {
    br <- unique(seq(min(pt, na.rm = TRUE), max(pt, na.rm = TRUE), length.out = n_bins + 1L))
  }
  if (length(br) < 2) stop("Could not create at least two time-bin breaks.")
  bins <- cut(pt, breaks = br, include.lowest = TRUE, labels = paste0("Bin", seq_len(length(br) - 1L)))
  names(bins) <- names(pseudotime)
  bins
}

.extract_comm_count_table <- function(ct, level = c("celltype_lr", "celltype_pair")) {
  stopifnot(inherits(ct, "LogicCommCellTypeComm"))
  level <- match.arg(level)
  df <- ct$lr_table
  # Build active counts in the SAME universe as the n_edges denominator.
  # lcs_unweighted = active / n_edges in both neighborhood mode (active
  # neighborhood edges / neighborhood edges) and global mode (active global
  # pairs / possible global pairs), so lcs_unweighted * n_edges recovers the
  # consistent active count. The previous code used n_active_edges, which for
  # distal candidates holds a global cell-count product against a neighborhood
  # n_edges denominator and was then silently clamped to 100% by pmin().
  n_edges <- as.numeric(df$n_edges)
  active_consistent <- as.numeric(df$lcs_unweighted) * n_edges
  scored <- is.finite(active_consistent) & is.finite(n_edges) & n_edges > 0
  if (level == "celltype_lr") {
    out <- data.frame(
      feature = paste(df$sender_type, df$receiver_type, df$lr_pair, sep = "|"),
      success = round(active_consistent),
      total = n_edges,
      lcs = as.numeric(df$lcs),
      stringsAsFactors = FALSE
    )
  } else {
    # n_edges is constant across LR pairs within a sender-receiver group, so the
    # pair-level binomial sums consistent per-LR (active, opportunity) counts.
    key <- paste(df$sender_type, df$receiver_type, sep = "|")[scored]
    agg_success <- tapply(round(active_consistent[scored]), key, sum)
    agg_total <- tapply(n_edges[scored], key, sum)
    agg_lcs <- tapply(as.numeric(df$lcs_unweighted[scored]), key, mean)
    feats <- names(agg_total)
    out <- data.frame(
      feature = feats,
      success = as.numeric(agg_success[feats]),
      total = as.numeric(agg_total[feats]),
      lcs = as.numeric(agg_lcs[feats]),
      stringsAsFactors = FALSE
    )
  }
  out <- out[is.finite(out$total) & out$total > 0 & is.finite(out$success) & out$success >= 0, , drop = FALSE]
  out$success <- pmin(out$success, out$total)
  out
}

.df_to_md <- function(df, digits = 4) {
  if (is.null(df) || nrow(df) == 0) return("No rows.")
  df <- as.data.frame(df)
  for (nm in names(df)) {
    if (is.numeric(df[[nm]])) df[[nm]] <- signif(df[[nm]], digits)
    df[[nm]] <- as.character(df[[nm]])
    df[[nm]][is.na(df[[nm]])] <- ""
  }
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

#' Print method for LogicCommGLMResult
#' @param x LogicCommGLMResult object.
#' @param ... Ignored.
#' @export
print.LogicCommGLMResult <- function(x, ...) {
  cat(sprintf("LogicCommGLMResult | %d features | coefficient: %s\n",
              nrow(x), attr(x, "coef") %||% "NA"))
  if (nrow(x) > 0) {
    print(utils::head(as.data.frame(x)[, intersect(c("feature", "estimate", "odds_ratio", "p_value", "fdr", "n_samples"), names(x)), drop = FALSE], 10))
  }
  invisible(x)
}
