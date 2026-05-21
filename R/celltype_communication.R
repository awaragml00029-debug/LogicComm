# R/celltype_communication.R

#' Summarize Cell-Type-Resolved Logic Communication
#'
#' Aggregates ligand-receptor Logic Consensus Scores (LCS) by directed
#' sender-cell-type to receiver-cell-type relationships. When a graph is supplied,
#' the result reports both local graph-supported evidence and global cell-type
#' communication potential so distal candidates are not discarded solely because
#' they are absent from a Seurat NN/SNN graph.
#'
#' @details
#' Self-loops and same-cell-type communication are intentionally separated.
#' With the default \\code{remove_self_edges = TRUE}, diagonal graph entries
#' such as cell A -> cell A are removed before neighborhood scoring, which avoids
#' treating within-cell ligand/receptor co-expression as direct intercellular
#' communication evidence. With the default \\code{include_self = TRUE}, sender
#' and receiver cell types are still allowed to be identical, so rows such as
#' Tumor -> Tumor or Macrophage -> Macrophage are retained and should be
#' interpreted as autocrine-like, homotypic, or within-cell-type signaling
#' potential. Global candidate scores provide complementary cell-type-level
#' potential, including signals that may be distal or not represented in the
#' local KNN/SNN graph.
#'
#' @param reo_mat Binary REO matrix (genes x cells).
#' @param cell_labels Named character/factor vector mapping cells to cell types.
#' @param seurat_obj Optional Seurat object.
#' @param label_col Optional Seurat metadata column for labels.
#' @param knn_mat Optional KNN adjacency matrix.
#' @param lr_db LogicComm LR database.
#' @param graph_name Optional Seurat graph name.
#' @param lcs_threshold Minimum cell-type LCS.
#' @param min_edges Minimum local/global opportunities required for active calls.
#' @param min_active_edges Minimum active local/global support required.
#' @param mode Scoring mode.
#' @param remove_self_edges Remove diagonal cell-level self-loops before local
#'   graph scoring. This does not remove same-cell-type communication.
#' @param graph_symmetrize Graph symmetrization.
#' @param edge_weight_mode Edge weighting.
#' @param min_role_hub_quantile Hub score quantile.
#' @param min_role_event_count Minimum event count.
#' @param include_self Include same-cell-type sender/receiver pairs, interpreted
#'   as autocrine-like or homotypic signaling potential.
#' @param verbose Print progress.
#'
#' @return A list of class \code{LogicCommCellTypeComm}.
#' @export
summarize_celltype_communication <- function(reo_mat,
                                             cell_labels = NULL,
                                             seurat_obj = NULL,
                                             label_col = NULL,
                                             knn_mat = NULL,
                                             lr_db = lr_pairs_human,
                                             graph_name = NULL,
                                             lcs_threshold = 0.01,
                                             min_edges = 20,
                                             min_active_edges = 1,
                                             mode = c("auto", "neighborhood", "global"),
                                             remove_self_edges = TRUE,
                                             graph_symmetrize = c("none", "or", "max"),
                                             edge_weight_mode = c("binary", "weighted"),
                                             min_role_hub_quantile = 0.2,
                                             min_role_event_count = 5,
                                             include_self = TRUE,
                                             verbose = TRUE) {
  if (is.list(reo_mat) && !is.null(reo_mat$logic)) reo_mat <- reo_mat$logic
  mode <- match.arg(mode)
  graph_symmetrize <- match.arg(graph_symmetrize)
  edge_weight_mode <- match.arg(edge_weight_mode)
  .validate_lr_db_for_celltype(lr_db)

  labels <- .resolve_celltype_labels(cell_labels, seurat_obj, label_col, colnames(reo_mat))
  cell_types <- sort(unique(labels))
  n_by_type <- table(factor(labels, levels = cell_types))
  celltype_sizes <- data.frame(cell_type = names(n_by_type), n_cells = as.integer(n_by_type), stringsAsFactors = FALSE)

  if (identical(mode, "auto")) {
    mode <- if (!is.null(knn_mat) || !is.null(seurat_obj)) "neighborhood" else "global"
  }

  group_defs <- expand.grid(sender_type = cell_types, receiver_type = cell_types, stringsAsFactors = FALSE)
  if (!include_self) group_defs <- group_defs[group_defs$sender_type != group_defs$receiver_type, , drop = FALSE]
  group_defs$key <- paste(group_defs$sender_type, group_defs$receiver_type, sep = "|||")
  n_groups <- nrow(group_defs)

  sender_n <- as.numeric(n_by_type[group_defs$sender_type])
  receiver_n <- as.numeric(n_by_type[group_defs$receiver_type])
  n_global_possible_by_group <- sender_n * receiver_n

  edge_i <- edge_j <- integer(0)
  edge_w <- numeric(0)
  edge_group_id <- integer(0)
  n_edges_by_group <- rep(0L, n_groups)
  edge_weight_sum_by_group <- rep(0, n_groups)

  if (mode == "neighborhood") {
    if (is.null(knn_mat)) knn_mat <- .extract_knn(seurat_obj, graph_name)
    knn_mat <- .celltype_validate_and_align_knn(knn_mat, colnames(reo_mat))
    knn_mat <- .celltype_symmetrize_sparse_graph(knn_mat, graph_symmetrize)
    knn_edges <- .celltype_sparse_to_edges(knn_mat, remove_self_edges, edge_weight_mode)

    edge_i <- knn_edges$i
    edge_j <- knn_edges$j
    edge_w <- knn_edges$w

    if (!include_self && length(edge_i) > 0) {
      keep <- labels[edge_i] != labels[edge_j]
      edge_i <- edge_i[keep]
      edge_j <- edge_j[keep]
      edge_w <- edge_w[keep]
    }

    if (length(edge_i) > 0) {
      edge_keys <- paste(labels[edge_i], labels[edge_j], sep = "|||")
      edge_group_id <- match(edge_keys, group_defs$key)
      valid_edge <- !is.na(edge_group_id)
      edge_group_id <- edge_group_id[valid_edge]
      edge_i <- edge_i[valid_edge]
      edge_j <- edge_j[valid_edge]
      edge_w <- edge_w[valid_edge]
      n_edges_by_group <- tabulate(edge_group_id, nbins = n_groups)
      edge_weight_sum_by_group <- .sum_by_group(edge_w, edge_group_id, n_groups)
    }
  } else {
    n_edges_by_group <- as.integer(n_global_possible_by_group)
    edge_weight_sum_by_group <- as.numeric(n_global_possible_by_group)
  }

  if (verbose) message("[CellTypeComm] Resolving unique complex logic...")
  unique_complexes <- unique(c(lr_db$ligand_genes, lr_db$receptor_genes))
  complex_keys <- vapply(unique_complexes, .complex_key, character(1))
  complex_logic_map <- lapply(unique_complexes, function(gs) .resolve_complex_logic(gs, reo_mat))
  names(complex_logic_map) <- complex_keys
  complex_type_sums <- lapply(complex_logic_map, .sum_logic_by_type, labels = labels, cell_types = cell_types)

  if (verbose) message(sprintf("[CellTypeComm] Scoring %d L-R pairs for %d CT pairs...", nrow(lr_db), n_groups))

  lr_table_list <- lapply(seq_len(nrow(lr_db)), function(idx) {
    lig_genes <- lr_db$ligand_genes[[idx]]
    rec_genes <- lr_db$receptor_genes[[idx]]

    if (!any(lig_genes %in% rownames(reo_mat)) || !any(rec_genes %in% rownames(reo_mat))) return(NULL)

    lig_key <- .complex_key(lig_genes)
    rec_key <- .complex_key(rec_genes)
    lig_logic <- complex_logic_map[[lig_key]]
    rec_logic <- complex_logic_map[[rec_key]]

    lcs_neighborhood <- rep(NA_real_, n_groups)
    lcs_unweighted <- rep(NA_real_, n_groups)
    n_active_neighborhood <- rep(0L, n_groups)
    active_edge_weight_sum <- rep(0, n_groups)

    if (mode == "neighborhood" && length(edge_i) > 0) {
      active_events <- lig_logic[edge_i] & rec_logic[edge_j]
      if (any(active_events)) {
        n_active_neighborhood <- tabulate(edge_group_id[active_events], nbins = n_groups)
        active_edge_weight_sum <- .sum_by_group(edge_w[active_events], edge_group_id[active_events], n_groups)
      }
      lcs_unweighted <- ifelse(n_edges_by_group > 0, n_active_neighborhood / n_edges_by_group, NA_real_)
      denom_nb <- if (edge_weight_mode == "weighted") edge_weight_sum_by_group else n_edges_by_group
      num_nb <- if (edge_weight_mode == "weighted") active_edge_weight_sum else n_active_neighborhood
      lcs_neighborhood <- ifelse(denom_nb > 0, num_nb / denom_nb, NA_real_)
    }

    lig_sum <- complex_type_sums[[lig_key]]
    rec_sum <- complex_type_sums[[rec_key]]
    ligand_active_n_sender <- as.numeric(lig_sum[group_defs$sender_type])
    receptor_active_n_receiver <- as.numeric(rec_sum[group_defs$receiver_type])
    ligand_active_frac_sender <- ligand_active_n_sender / pmax(sender_n, 1)
    receptor_active_frac_receiver <- receptor_active_n_receiver / pmax(receiver_n, 1)
    n_active_global <- ligand_active_n_sender * receptor_active_n_receiver
    lcs_global <- n_active_global / pmax(n_global_possible_by_group, 1)
    if (mode == "global") lcs_unweighted <- lcs_global

    local_denominator_ok <- mode == "neighborhood" & !is.na(lcs_neighborhood) & n_edges_by_group >= min_edges
    local_active <- local_denominator_ok & lcs_neighborhood >= lcs_threshold & n_active_neighborhood >= min_active_edges

    global_denominator_ok <- n_global_possible_by_group >= min_edges
    global_candidate_active <- global_denominator_ok & !is.na(lcs_global) & lcs_global >= lcs_threshold & n_active_global >= min_active_edges

    distal_candidate <- mode == "neighborhood" & global_candidate_active & !local_active
    candidate_active <- local_active | global_candidate_active

    lcs_final <- ifelse(local_active, lcs_neighborhood,
                        ifelse(global_candidate_active, lcs_global,
                               if (mode == "neighborhood") lcs_neighborhood else lcs_global))
    lcs_primary_mode <- ifelse(local_active, "local",
                               ifelse(distal_candidate, "distal_global",
                                      ifelse(global_candidate_active, "global", "unscored")))
    comm_range <- ifelse(local_active & global_candidate_active, "paracrine",
                         ifelse(local_active & !global_candidate_active, "juxtacrine",
                                ifelse(global_candidate_active & !local_active, "distal/endocrine", "inactive")))

    n_active_edges <- ifelse(local_active, n_active_neighborhood,
                             ifelse(global_candidate_active, n_active_global, n_active_neighborhood))

    data.frame(
      sender_type = group_defs$sender_type,
      receiver_type = group_defs$receiver_type,
      lr_pair = .lr_scalar(lr_db, "lr_pair", idx, paste(lig_genes, rec_genes, sep = "_")),
      ligand = .lr_scalar(lr_db, "ligand", idx, paste(lig_genes, collapse = "+")),
      receptor = .lr_scalar(lr_db, "receptor", idx, paste(rec_genes, collapse = "+")),
      pathway = .lr_scalar(lr_db, "pathway", idx, "Unknown"),
      lcs = as.numeric(lcs_final),
      lcs_neighborhood = as.numeric(lcs_neighborhood),
      lcs_global = as.numeric(lcs_global),
      lcs_unweighted = as.numeric(lcs_unweighted),
      lcs_primary_mode = as.character(lcs_primary_mode),
      communication_range = as.character(comm_range),
      n_edges = as.integer(n_edges_by_group),
      edge_weight_sum = as.numeric(edge_weight_sum_by_group),
      n_active_edges = as.numeric(n_active_edges),
      n_active_neighborhood = as.integer(n_active_neighborhood),
      active_edge_weight_sum = as.numeric(active_edge_weight_sum),
      n_global_possible = as.numeric(n_global_possible_by_group),
      n_active_global = as.numeric(n_active_global),
      sender_cell_count = as.integer(sender_n),
      receiver_cell_count = as.integer(receiver_n),
      ligand_active_n_sender = as.numeric(ligand_active_n_sender),
      receptor_active_n_receiver = as.numeric(receptor_active_n_receiver),
      ligand_active_frac_sender = as.numeric(ligand_active_frac_sender),
      receptor_active_frac_receiver = as.numeric(receptor_active_frac_receiver),
      local_active = as.logical(local_active),
      global_candidate_active = as.logical(global_candidate_active),
      distal_candidate = as.logical(distal_candidate),
      candidate_active = as.logical(candidate_active),
      active = as.logical(candidate_active),
      stringsAsFactors = FALSE
    )
  })

  lr_table <- do.call(rbind, lr_table_list[!vapply(lr_table_list, is.null, logical(1))])
  if (is.null(lr_table)) lr_table <- .empty_lr_table()

  pair_summary <- .summarize_celltype_pairs(lr_table, group_defs)
  pathway_summary <- .summarize_celltype_pathways(lr_table)

  adjacency_strength <- .pair_summary_to_matrix(pair_summary, cell_types, "sum_lcs")
  adjacency_count <- .pair_summary_to_matrix(pair_summary, cell_types, "n_active_lr")
  adjacency_active_edge_support <- .pair_summary_to_matrix(pair_summary, cell_types, "sum_active_edges")

  role_summary <- .communication_role_summary(
    adjacency_strength,
    adjacency_count,
    active_edge_mat = adjacency_active_edge_support,
    lr_table = lr_table,
    min_role_hub_quantile = min_role_hub_quantile,
    min_role_event_count = min_role_event_count
  )

  res <- list(
    lr_table = lr_table,
    pair_summary = pair_summary,
    pathway_summary = pathway_summary,
    role_summary = role_summary,
    celltype_sizes = celltype_sizes,
    adjacency_strength = adjacency_strength,
    adjacency_count = adjacency_count,
    adjacency_active_edge_support = adjacency_active_edge_support,
    cell_labels = labels,
    lr_db = lr_db,
    params = list(
      mode = mode,
      lcs_threshold = lcs_threshold,
      min_edges = min_edges,
      min_active_edges = min_active_edges,
      n_cells = length(labels),
      n_cell_types = length(cell_types),
      n_lr_pairs = nrow(lr_db)
    )
  )
  class(res) <- "LogicCommCellTypeComm"

  if (verbose) message(sprintf("[CellTypeComm] Done. %d active L-R events.", sum(lr_table$active, na.rm = TRUE)))
  res
}

#' Convert cell-type communication summaries to named LCS vectors
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param level Feature level.
#' @param metric Numeric metric to extract.
#' @param active_only Whether to keep active rows for LR-level features.
#' @return Named numeric vector.
#' @export
celltype_comm_to_lcs <- function(ct_comm,
                                 level = c("celltype_lr", "celltype_pair", "pathway_pair"),
                                 metric = NULL,
                                 active_only = FALSE) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  level <- match.arg(level)
  if (level == "celltype_lr") {
    df <- ct_comm$lr_table
    if (isTRUE(active_only) && "active" %in% names(df)) df <- df[df$active %in% TRUE, , drop = FALSE]
    metric <- if (is.null(metric)) "lcs" else metric
    keys <- paste(df$sender_type, df$receiver_type, df$lr_pair, sep = "|")
  } else if (level == "celltype_pair") {
    df <- ct_comm$pair_summary
    metric <- if (is.null(metric)) "sum_lcs" else metric
    keys <- paste(df$sender_type, df$receiver_type, sep = "|")
  } else {
    df <- ct_comm$pathway_summary
    metric <- if (is.null(metric)) "sum_lcs" else metric
    keys <- paste(df$sender_type, df$receiver_type, df$pathway, sep = "|")
  }
  if (is.null(df) || nrow(df) == 0) return(stats::setNames(numeric(0), character(0)))
  if (!metric %in% names(df)) stop("metric not found: ", metric)
  if (!is.numeric(df[[metric]])) stop("metric must be numeric: ", metric)
  stats::setNames(as.numeric(df[[metric]]), keys)
}

.resolve_celltype_labels <- function(cell_labels, seurat_obj, label_col, cell_names) {
  if (inherits(cell_labels, "Seurat")) {
    seurat_obj <- cell_labels
    cell_labels <- NULL
  }
  if (is.null(cell_labels)) {
    if (is.null(seurat_obj)) stop("cell_labels or seurat_obj must be provided for cell-type-resolved communication.")
    if (!is.null(label_col)) {
      if (!label_col %in% colnames(seurat_obj@meta.data)) stop("label_col '", label_col, "' not found in Seurat metadata.")
      cell_labels <- seurat_obj@meta.data[[label_col]]
      names(cell_labels) <- rownames(seurat_obj@meta.data)
    } else {
      cell_labels <- SeuratObject::Idents(seurat_obj)
      names(cell_labels) <- colnames(seurat_obj)
    }
  }
  label_names <- names(cell_labels)
  cell_labels <- as.character(cell_labels)
  if (!is.null(label_names)) names(cell_labels) <- label_names
  if (is.null(names(cell_labels))) {
    if (length(cell_labels) != length(cell_names)) stop("Unnamed cell_labels must have length ncol(reo_mat).")
    names(cell_labels) <- cell_names
  }

  missing_cells <- setdiff(cell_names, names(cell_labels))
  if (length(missing_cells) > 0) {
    stop(
      "cell_labels is missing labels for ", length(missing_cells),
      " cells in reo_mat; first missing cells: ",
      paste(utils::head(missing_cells, 5), collapse = ", "),
      "."
    )
  }

  labels <- cell_labels[cell_names]
  invalid <- is.na(labels) | !nzchar(labels)
  if (any(invalid)) {
    invalid_cells <- cell_names[invalid]
    stop(
      "cell_labels contains ", sum(invalid),
      " missing or empty label values; first affected cells: ",
      paste(utils::head(invalid_cells, 5), collapse = ", "),
      ". Remove these cells or assign valid cell type labels before summarizing communication."
    )
  }
  labels
}

.sum_by_group <- function(values, group, n_groups) {
  out <- numeric(n_groups)
  if (length(values) == 0 || length(group) == 0) return(out)
  ss <- tapply(values, group, sum, na.rm = TRUE)
  out[as.integer(names(ss))] <- as.numeric(ss)
  out
}

.sum_logic_by_type <- function(log_vec, labels, cell_types) {
  vals <- tapply(as.numeric(log_vec), factor(labels, levels = cell_types), sum, na.rm = TRUE)
  out <- as.numeric(vals)
  names(out) <- cell_types
  out[is.na(out)] <- 0
  out
}

.summarize_celltype_pairs <- function(lr_table, group_defs) {
  if (nrow(lr_table) == 0) return(.empty_pair_summary(group_defs))
  keys <- paste(lr_table$sender_type, lr_table$receiver_type, sep = "|||")
  sp <- split(seq_len(nrow(lr_table)), keys)

  out <- lapply(sp, function(ii) {
    sub <- lr_table[ii, , drop = FALSE]
    active_sub <- sub[sub$active %in% TRUE, , drop = FALSE]
    finite_sub <- sub[is.finite(sub$lcs), , drop = FALSE]
    support_frac <- if (nrow(active_sub) > 0) active_sub$n_active_edges / pmax(active_sub$n_edges, 1) else numeric(0)
    top_lr <- if (nrow(active_sub) > 0) active_sub$lr_pair[which.max(active_sub$lcs)] else NA_character_
    top_pathway <- if (nrow(active_sub) > 0) {
      pw_sum <- tapply(active_sub$lcs, active_sub$pathway, sum, na.rm = TRUE)
      names(pw_sum)[which.max(pw_sum)]
    } else NA_character_
    dominant_range <- if (nrow(active_sub) > 0) {
      ranges <- table(active_sub$communication_range)
      names(ranges)[which.max(ranges)]
    } else NA_character_

    data.frame(
      sender_type = sub$sender_type[1],
      receiver_type = sub$receiver_type[1],
      n_edges = as.integer(sub$n_edges[1]),
      edge_weight_sum = as.numeric(sub$edge_weight_sum[1]),
      n_scored_lr = nrow(finite_sub),
      n_active_lr = nrow(active_sub),
      active_lr_event_count = nrow(active_sub),
      n_local_active = sum(sub$local_active %in% TRUE),
      n_global_candidate_active = sum(sub$global_candidate_active %in% TRUE),
      n_distal_candidate = sum(sub$distal_candidate %in% TRUE),
      n_juxtacrine = sum(active_sub$communication_range == "juxtacrine"),
      n_paracrine = sum(active_sub$communication_range == "paracrine"),
      n_distal = sum(active_sub$communication_range == "distal/endocrine"),
      sum_lcs = sum(active_sub$lcs, na.rm = TRUE),
      sum_lcs_all = sum(finite_sub$lcs, na.rm = TRUE),
      mean_lcs_active = if (nrow(active_sub) > 0) mean(active_sub$lcs, na.rm = TRUE) else NA_real_,
      mean_lcs_all = if (nrow(finite_sub) > 0) mean(finite_sub$lcs, na.rm = TRUE) else NA_real_,
      sum_active_edges = sum(active_sub$n_active_edges, na.rm = TRUE),
      sum_active_edge_weight = sum(active_sub$active_edge_weight_sum, na.rm = TRUE),
      active_edge_weight_sum = sum(active_sub$active_edge_weight_sum, na.rm = TRUE),
      mean_edge_support_fraction_active = if (length(support_frac) > 0) mean(support_frac, na.rm = TRUE) else NA_real_,
      dominant_communication_range = dominant_range,
      top_lr_pair = top_lr,
      top_pathway = top_pathway,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  missing_keys <- setdiff(group_defs$key, paste(out$sender_type, out$receiver_type, sep = "|||"))
  if (length(missing_keys) > 0) {
    out <- rbind(out, .empty_pair_summary(group_defs[group_defs$key %in% missing_keys, , drop = FALSE]))
  }
  out[order(out$sender_type, out$receiver_type), , drop = FALSE]
}

.summarize_celltype_pathways <- function(lr_table) {
  df <- lr_table[lr_table$active %in% TRUE, , drop = FALSE]
  if (nrow(df) == 0) return(.empty_pathway_summary())

  keys <- paste(df$sender_type, df$receiver_type, df$pathway, sep = "|||")
  sp <- split(seq_len(nrow(df)), keys)

  out <- lapply(sp, function(ii) {
    sub <- df[ii, , drop = FALSE]
    support_frac <- sub$n_active_edges / pmax(sub$n_edges, 1)
    ranges <- table(sub$communication_range)
    top_lr <- sub$lr_pair[which.max(sub$lcs)]
    data.frame(
      sender_type = sub$sender_type[1],
      receiver_type = sub$receiver_type[1],
      pathway = sub$pathway[1],
      n_active_lr = nrow(sub),
      n_local_active = sum(sub$local_active %in% TRUE),
      n_global_candidate_active = sum(sub$global_candidate_active %in% TRUE),
      n_distal_candidate = sum(sub$distal_candidate %in% TRUE),
      sum_lcs = sum(sub$lcs, na.rm = TRUE),
      sum_lcs_all = sum(sub$lcs, na.rm = TRUE),
      mean_lcs_active = mean(sub$lcs, na.rm = TRUE),
      sum_active_edges = sum(sub$n_active_edges, na.rm = TRUE),
      sum_active_edge_weight = sum(sub$active_edge_weight_sum, na.rm = TRUE),
      mean_edge_support_fraction_active = mean(support_frac, na.rm = TRUE),
      dominant_communication_range = names(ranges)[which.max(ranges)],
      top_lr_pair = top_lr,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  out[order(out$sender_type, out$receiver_type, out$pathway), , drop = FALSE]
}

.pair_summary_to_matrix <- function(pair_summary, cell_types, metric) {
  mat <- matrix(0, nrow = length(cell_types), ncol = length(cell_types), dimnames = list(cell_types, cell_types))
  if (nrow(pair_summary) == 0 || !metric %in% names(pair_summary)) return(mat)
  for (i in seq_len(nrow(pair_summary))) {
    s <- pair_summary$sender_type[i]
    r <- pair_summary$receiver_type[i]
    val <- pair_summary[[metric]][i]
    if (!is.na(val) && s %in% rownames(mat) && r %in% colnames(mat)) mat[s, r] <- val
  }
  mat
}

.extract_knn <- function(seurat_obj, graph_name) {
  if (is.null(seurat_obj)) return(NULL)
  if (is.null(graph_name)) {
    graph_name <- grep("_nn$", names(seurat_obj@graphs), value = TRUE)[1]
    if (is.na(graph_name)) graph_name <- grep("snn$", names(seurat_obj@graphs), value = TRUE)[1]
    if (is.na(graph_name)) graph_name <- names(seurat_obj@graphs)[1]
  }
  if (is.null(graph_name) || !graph_name %in% names(seurat_obj@graphs)) stop("Could not find graph in Seurat object. Please specify graph_name.")
  seurat_obj@graphs[[graph_name]]
}

.celltype_validate_and_align_knn <- function(knn_mat, cell_names) {
  if (is.null(knn_mat)) return(NULL)
  if (is.null(dim(knn_mat)) || nrow(knn_mat) != ncol(knn_mat)) stop("knn_mat must be a square cells x cells adjacency matrix.")
  if (!is.null(rownames(knn_mat)) && !is.null(colnames(knn_mat))) {
    if (!all(cell_names %in% rownames(knn_mat)) || !all(cell_names %in% colnames(knn_mat))) stop("KNN matrix is missing cells from REO matrix.")
    knn_mat <- knn_mat[cell_names, cell_names, drop = FALSE]
  } else if (nrow(knn_mat) != length(cell_names)) {
    stop("knn_mat dimensions must match ncol(reo_mat) when names are absent.")
  } else {
    rownames(knn_mat) <- colnames(knn_mat) <- cell_names
  }
  if (!inherits(knn_mat, "sparseMatrix")) knn_mat <- Matrix::Matrix(knn_mat, sparse = TRUE)
  knn_mat
}

.celltype_symmetrize_sparse_graph <- function(m, mode) {
  if (is.null(m) || mode == "none") return(m)
  if (!inherits(m, "sparseMatrix")) m <- Matrix::Matrix(m, sparse = TRUE)
  m_t <- Matrix::t(m)
  if (mode == "or") return((m + m_t) > 0)
  if (mode == "max") {
    s1 <- Matrix::summary(m)
    s2 <- Matrix::summary(m_t)
    df <- rbind(
      data.frame(i = s1$i, j = s1$j, x = s1$x),
      data.frame(i = s2$i, j = s2$j, x = s2$x)
    )
    if (!nrow(df)) return(m)
    agg <- stats::aggregate(x ~ i + j, df, max)
    return(Matrix::sparseMatrix(i = agg$i, j = agg$j, x = agg$x,
                                dims = dim(m), dimnames = dimnames(m)))
  }
  m
}

.celltype_sparse_to_edges <- function(m, remove_self, weight_mode) {
  if (is.null(m)) return(list(i = integer(0), j = integer(0), w = numeric(0)))
  if (!inherits(m, "sparseMatrix")) m <- Matrix::Matrix(m, sparse = TRUE)
  if (isTRUE(remove_self) && nrow(m) > 0) {
    diag(m) <- 0
    m <- Matrix::drop0(m)
  }
  summary_m <- Matrix::summary(m)
  if (nrow(summary_m) == 0) return(list(i = integer(0), j = integer(0), w = numeric(0)))
  w <- if (weight_mode == "binary") rep(1, nrow(summary_m)) else as.numeric(summary_m$x)
  w[!is.finite(w)] <- 0
  keep <- w > 0
  list(i = as.integer(summary_m$i[keep]), j = as.integer(summary_m$j[keep]), w = as.numeric(w[keep]))
}

.validate_lr_db_for_celltype <- function(lr_db) {
  required <- c("lr_pair", "ligand_genes", "receptor_genes")
  missing <- setdiff(required, names(lr_db))
  if (length(missing) > 0) stop("lr_db is missing required columns: ", paste(missing, collapse = ", "))
  if (!is.list(lr_db$ligand_genes) || !is.list(lr_db$receptor_genes)) stop("lr_db ligand_genes and receptor_genes must be list columns.")
}

.complex_key <- function(gs) paste(sort(unique(as.character(gs))), collapse = "|")

.lr_scalar <- function(lr_db, col, idx, default = NA_character_) {
  if (!col %in% names(lr_db)) return(default)
  val <- lr_db[[col]][idx]
  if (length(val) == 0 || is.na(val) || !nzchar(as.character(val))) default else as.character(val)
}

.empty_lr_table <- function() {
  data.frame(
    sender_type = character(), receiver_type = character(), lr_pair = character(), ligand = character(), receptor = character(), pathway = character(),
    lcs = numeric(), lcs_neighborhood = numeric(), lcs_global = numeric(), lcs_unweighted = numeric(), lcs_primary_mode = character(), communication_range = character(),
    n_edges = integer(), edge_weight_sum = numeric(), n_active_edges = numeric(), n_active_neighborhood = integer(), active_edge_weight_sum = numeric(),
    n_global_possible = numeric(), n_active_global = numeric(), sender_cell_count = integer(), receiver_cell_count = integer(), ligand_active_n_sender = numeric(),
    receptor_active_n_receiver = numeric(), ligand_active_frac_sender = numeric(), receptor_active_frac_receiver = numeric(), local_active = logical(),
    global_candidate_active = logical(), distal_candidate = logical(), candidate_active = logical(), active = logical(), stringsAsFactors = FALSE
  )
}

.empty_pair_summary <- function(group_defs) {
  data.frame(
    sender_type = group_defs$sender_type,
    receiver_type = group_defs$receiver_type,
    n_edges = integer(nrow(group_defs)),
    edge_weight_sum = numeric(nrow(group_defs)),
    n_scored_lr = integer(nrow(group_defs)),
    n_active_lr = integer(nrow(group_defs)),
    active_lr_event_count = integer(nrow(group_defs)),
    n_local_active = integer(nrow(group_defs)),
    n_global_candidate_active = integer(nrow(group_defs)),
    n_distal_candidate = integer(nrow(group_defs)),
    n_juxtacrine = integer(nrow(group_defs)),
    n_paracrine = integer(nrow(group_defs)),
    n_distal = integer(nrow(group_defs)),
    sum_lcs = numeric(nrow(group_defs)),
    sum_lcs_all = numeric(nrow(group_defs)),
    mean_lcs_active = rep(NA_real_, nrow(group_defs)),
    mean_lcs_all = rep(NA_real_, nrow(group_defs)),
    sum_active_edges = numeric(nrow(group_defs)),
    sum_active_edge_weight = numeric(nrow(group_defs)),
    active_edge_weight_sum = numeric(nrow(group_defs)),
    mean_edge_support_fraction_active = rep(NA_real_, nrow(group_defs)),
    dominant_communication_range = rep(NA_character_, nrow(group_defs)),
    top_lr_pair = rep(NA_character_, nrow(group_defs)),
    top_pathway = rep(NA_character_, nrow(group_defs)),
    stringsAsFactors = FALSE
  )
}

.empty_pathway_summary <- function() {
  data.frame(
    sender_type = character(), receiver_type = character(), pathway = character(), n_active_lr = integer(), n_local_active = integer(),
    n_global_candidate_active = integer(), n_distal_candidate = integer(), sum_lcs = numeric(), sum_lcs_all = numeric(), mean_lcs_active = numeric(),
    sum_active_edges = numeric(), sum_active_edge_weight = numeric(), mean_edge_support_fraction_active = numeric(), dominant_communication_range = character(),
    top_lr_pair = character(), stringsAsFactors = FALSE
  )
}

#' @export
print.LogicCommCellTypeComm <- function(x, ...) {
  cat(sprintf("LogicCommCellTypeComm | %d cells | %d cell types | %d active L-R events\n",
              length(x$cell_labels),
              nrow(x$role_summary),
              sum(x$lr_table$active, na.rm = TRUE)))
  invisible(x)
}
