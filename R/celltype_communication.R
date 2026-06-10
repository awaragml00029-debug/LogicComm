# R/celltype_communication.R

#' Summarize Cell-Type-Resolved Logic Communication
#'
#' Aggregates ligand-receptor Logic Consensus Scores (LCS) by directed
#' sender-cell-type to receiver-cell-type relationships. Communication is scored
#' at the cell-type level from REO co-expression: for each sender -> receiver
#' cell-type pair, the LCS is the fraction of the pair's opportunity universe
#' (the sender x receiver cell-count product) in which the ligand is active in
#' the sender cell and the receptor is active in the receiver cell.
#'
#' @details
#' LogicComm scores communication at the cell-type level and does not aggregate
#' over a per-cell KNN/SNN neighborhood graph. For dissociated scRNA-seq that
#' graph lives in expression space rather than physical space, so it cannot
#' license spatial juxtacrine/paracrine distance claims; cell-type co-expression
#' is the honest, comparable unit. With the default \code{include_self = TRUE},
#' sender and receiver cell types are allowed to be identical, so rows such as
#' Tumor -> Tumor or Macrophage -> Macrophage are retained and interpreted as
#' autocrine-like or homotypic signaling potential.
#'
#' Legacy neighborhood arguments (\code{knn_mat}, \code{graph_name},
#' \code{mode}, \code{remove_self_edges}, \code{graph_symmetrize},
#' \code{edge_weight_mode}) are accepted via \code{...} for backward
#' compatibility but are deprecated and ignored with a warning.
#'
#' @param reo_mat Binary REO matrix (genes x cells).
#' @param cell_labels Named character/factor vector mapping cells to cell types.
#' @param seurat_obj Optional Seurat object (used only to resolve cell labels).
#' @param label_col Optional Seurat metadata column for labels.
#' @param lr_db LogicComm LR database.
#' @param lcs_threshold Minimum cell-type LCS for an axis to be called active.
#' @param min_edges Minimum opportunity universe (sender x receiver cell-count
#'   product) required for an active call.
#' @param min_active_edges Minimum co-expressing support required for an active call.
#' @param min_expr_frac Minimum fraction of sender cells that must express the
#'   ligand and of receiver cells that must express the receptor for a
#'   sender-receiver L-R axis to be called active. Prevents a few high-expressing
#'   cells from calling an axis active when the gene is detected in a negligible
#'   fraction of the cell type. Default: \code{0.1}; set \code{0} to disable.
#' @param min_role_hub_quantile Hub score quantile.
#' @param min_role_event_count Minimum event count.
#' @param include_self Include same-cell-type sender/receiver pairs, interpreted
#'   as autocrine-like or homotypic signaling potential.
#' @param lcs_weighting Communication score. \code{"binary"} (default) is the
#'   co-expression fraction product (binary REO). \code{"rank"} weights by REO
#'   intensity -- the prevalence-weighted within-cell rank of the ligand in the
#'   sender type (fraction expressing times mean rank among expressers) times the
#'   same for the receptor -- which keeps the prevalence signal that separates
#'   specific from ubiquitous axes while adding dynamic range. \code{"rank"}
#'   requires the input to carry the rank matrix, i.e. build it with
#'   \code{calc_REO_matrix(..., return_rank = TRUE)}. The \code{active} call is
#'   unchanged (binary), so only the reported \code{lcs} differs.
#' @param verbose Print progress.
#' @param ... Deprecated neighborhood arguments, accepted for backward
#'   compatibility and ignored with a warning.
#'
#' @return A list of class \code{LogicCommCellTypeComm}.
#' @export
summarize_celltype_communication <- function(reo_mat,
                                             cell_labels = NULL,
                                             seurat_obj = NULL,
                                             label_col = NULL,
                                             lr_db = lr_pairs_human,
                                             lcs_threshold = 0.01,
                                             min_edges = 20,
                                             min_active_edges = 1,
                                             min_expr_frac = 0.1,
                                             min_role_hub_quantile = 0.2,
                                             min_role_event_count = 5,
                                             include_self = TRUE,
                                             lcs_weighting = c("binary", "rank"),
                                             verbose = TRUE,
                                             ...) {
  lcs_weighting <- match.arg(lcs_weighting)
  rank_mat <- NULL
  if (is.list(reo_mat) && !is.null(reo_mat$logic)) {
    rank_mat <- reo_mat$rank
    reo_mat <- reo_mat$logic
  }
  .deprecate_neighborhood_args(list(...))
  if (lcs_weighting == "rank" && is.null(rank_mat)) {
    stop("lcs_weighting = 'rank' needs the within-cell REO rank matrix; build the ",
         "input with calc_REO_matrix(..., return_rank = TRUE).", call. = FALSE)
  }
  stopifnot(is.numeric(min_expr_frac), length(min_expr_frac) == 1,
            min_expr_frac >= 0, min_expr_frac <= 1)
  .validate_lr_db_for_celltype(lr_db)

  labels <- .resolve_celltype_labels(cell_labels, seurat_obj, label_col, colnames(reo_mat))
  invalid_labels <- is.na(labels) | !nzchar(labels)
  if (any(invalid_labels)) {
    n_filtered <- sum(invalid_labels)
    filtered_cells <- names(labels)[invalid_labels]
    if (n_filtered == length(labels)) {
      stop(
        "All cell_labels are missing or empty; assign valid cell type labels before summarizing communication."
      )
    }
    if (isTRUE(verbose)) {
      message(
        "Filtered ", n_filtered,
        " cells with missing or empty cell type labels before communication summary; first affected cells: ",
        paste(utils::head(filtered_cells, 5), collapse = ", "),
        "."
      )
    }
    labels <- labels[!invalid_labels]
    reo_mat <- reo_mat[, names(labels), drop = FALSE]
  }
  if (!is.null(rank_mat)) rank_mat <- rank_mat[, colnames(reo_mat), drop = FALSE]
  cell_types <- sort(unique(labels))
  n_by_type <- table(factor(labels, levels = cell_types))
  celltype_sizes <- data.frame(cell_type = names(n_by_type), n_cells = as.integer(n_by_type), stringsAsFactors = FALSE)

  group_defs <- expand.grid(sender_type = cell_types, receiver_type = cell_types, stringsAsFactors = FALSE)
  if (!include_self) group_defs <- group_defs[group_defs$sender_type != group_defs$receiver_type, , drop = FALSE]
  group_defs$key <- paste(group_defs$sender_type, group_defs$receiver_type, sep = "|||")
  n_groups <- nrow(group_defs)

  sender_n <- as.numeric(n_by_type[group_defs$sender_type])
  receiver_n <- as.numeric(n_by_type[group_defs$receiver_type])
  n_global_possible_by_group <- sender_n * receiver_n

  # Cell-type co-expression scoring. The LCS denominator is the cell-type-pair
  # opportunity (sender x receiver cell-count product). LogicComm scores
  # communication at the cell-type level and no longer aggregates over a
  # per-cell KNN/SNN neighborhood graph: for dissociated scRNA-seq that graph
  # lives in expression space, not physical space, so it cannot license spatial
  # juxtacrine/paracrine distance claims.
  n_edges_by_group <- as.integer(n_global_possible_by_group)
  edge_weight_sum_by_group <- as.numeric(n_global_possible_by_group)

  if (verbose) message("[CellTypeComm] Resolving unique complex logic...")
  unique_complexes <- unique(c(lr_db$ligand_genes, lr_db$receptor_genes))
  complex_keys <- vapply(unique_complexes, .complex_key, character(1))
  complex_logic_map <- lapply(unique_complexes, function(gs) .resolve_complex_logic(gs, reo_mat))
  names(complex_logic_map) <- complex_keys
  complex_type_sums <- lapply(complex_logic_map, .sum_logic_by_type, labels = labels, cell_types = cell_types)
  complex_type_meanrank <- NULL
  if (lcs_weighting == "rank") {
    complex_rank_map <- lapply(unique_complexes, function(gs) .resolve_complex_rank(gs, rank_mat))
    names(complex_rank_map) <- complex_keys
    complex_type_meanrank <- stats::setNames(
      Map(function(lg, rk) .mean_rank_by_type(rk, lg, labels, cell_types),
          complex_logic_map, complex_rank_map),
      complex_keys)
  }

  if (verbose) message(sprintf("[CellTypeComm] Scoring %d L-R pairs for %d CT pairs...", nrow(lr_db), n_groups))

  lr_table_list <- lapply(seq_len(nrow(lr_db)), function(idx) {
    lig_genes <- lr_db$ligand_genes[[idx]]
    rec_genes <- lr_db$receptor_genes[[idx]]

    if (!any(lig_genes %in% rownames(reo_mat)) || !any(rec_genes %in% rownames(reo_mat))) return(NULL)

    lig_key <- .complex_key(lig_genes)
    rec_key <- .complex_key(rec_genes)
    lig_sum <- complex_type_sums[[lig_key]]
    rec_sum <- complex_type_sums[[rec_key]]
    ligand_active_n_sender <- as.numeric(lig_sum[group_defs$sender_type])
    receptor_active_n_receiver <- as.numeric(rec_sum[group_defs$receiver_type])
    ligand_active_frac_sender <- ligand_active_n_sender / pmax(sender_n, 1)
    receptor_active_frac_receiver <- receptor_active_n_receiver / pmax(receiver_n, 1)
    n_active_global <- ligand_active_n_sender * receptor_active_n_receiver
    lcs_global <- n_active_global / pmax(n_global_possible_by_group, 1)

    # Communication score. "binary" (default) is the co-expression fraction
    # product. "rank" weights by REO intensity: the prevalence-weighted within-cell
    # rank of the ligand in the sender type (fraction expressing x mean rank among
    # expressers) times the same for the receptor -- this keeps the prevalence
    # signal that separates specific from ubiquitous axes while adding dynamic
    # range, so a ligand at the 99th within-cell percentile separates from the 51st.
    if (lcs_weighting == "rank") {
      lig_str <- as.numeric(complex_type_meanrank[[lig_key]][group_defs$sender_type])
      rec_str <- as.numeric(complex_type_meanrank[[rec_key]][group_defs$receiver_type])
      lcs_value <- lig_str * rec_str
    } else {
      lcs_value <- lcs_global
    }

    # Cell-type expressing-fraction gate: the ligand must be active in at least
    # min_expr_frac of sender cells and the receptor in at least min_expr_frac of
    # receiver cells. This prevents a few high-expressing cells from calling an
    # axis "active" when the gene is detected in a negligible fraction of the
    # cell type (e.g. ambient CD8A in a few Treg cells). The active call uses the
    # binary co-expression fraction, so the active set is identical under either
    # weighting; only the reported lcs score differs.
    expr_ok <- ligand_active_frac_sender >= min_expr_frac &
               receptor_active_frac_receiver >= min_expr_frac

    active <- n_global_possible_by_group >= min_edges & !is.na(lcs_global) &
              lcs_global >= lcs_threshold & n_active_global >= min_active_edges & expr_ok

    # Cell-type co-expression scoring. `lcs` is the reported score (binary
    # co-expression fraction by default, or the rank-weighted score when
    # lcs_weighting = "rank"); `lcs_global` is always the binary co-expression
    # fraction (ligand-active fraction of sender x receptor-active fraction of
    # receiver) and `lcs_unweighted` is kept as its alias. `n_edges` /
    # `edge_weight_sum` are the cell-type-pair opportunity universe
    # (sender x receiver cell-count product); `n_active_edges` is the
    # co-expressing support count. There is no per-cell neighborhood graph, so
    # no local/distal/range fields are produced.
    lcs_unweighted <- lcs_global
    n_active_edges <- n_active_global

    data.frame(
      sender_type = group_defs$sender_type,
      receiver_type = group_defs$receiver_type,
      lr_pair = .lr_scalar(lr_db, "lr_pair", idx, paste(lig_genes, rec_genes, sep = "_")),
      ligand = .lr_scalar(lr_db, "ligand", idx, paste(lig_genes, collapse = "+")),
      receptor = .lr_scalar(lr_db, "receptor", idx, paste(rec_genes, collapse = "+")),
      pathway = .lr_scalar(lr_db, "pathway", idx, "Unknown"),
      lcs = as.numeric(lcs_value),
      lcs_global = as.numeric(lcs_global),
      lcs_unweighted = as.numeric(lcs_unweighted),
      n_edges = as.integer(n_edges_by_group),
      edge_weight_sum = as.numeric(edge_weight_sum_by_group),
      n_active_edges = as.numeric(n_active_edges),
      n_global_possible = as.numeric(n_global_possible_by_group),
      n_active_global = as.numeric(n_active_global),
      sender_cell_count = as.integer(sender_n),
      receiver_cell_count = as.integer(receiver_n),
      ligand_active_n_sender = as.numeric(ligand_active_n_sender),
      receptor_active_n_receiver = as.numeric(receptor_active_n_receiver),
      ligand_active_frac_sender = as.numeric(ligand_active_frac_sender),
      receptor_active_frac_receiver = as.numeric(receptor_active_frac_receiver),
      active = as.logical(active),
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
      mode = "celltype",
      lcs_weighting = lcs_weighting,
      include_self = include_self,
      lcs_threshold = lcs_threshold,
      min_edges = min_edges,
      min_active_edges = min_active_edges,
      min_expr_frac = min_expr_frac,
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

  cell_labels[cell_names]
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

# Per-cell rank of a (possibly multi-subunit) complex is resolved by the single
# canonical .resolve_complex_rank() defined in R/modulators.R (conservative "min"
# aggregation across subunits by default, matching the AND logic of
# .resolve_complex_logic). It returns NA for cells where a required subunit is
# absent; .mean_rank_by_type() below only reads ranks at cells where the complex
# logic is active (all subunits present, hence finite), so the result is
# unaffected.

# Prevalence-weighted complex rank per cell type: the summed complex rank over
# the cells where the complex is active, divided by ALL cells of the type. This
# is (fraction expressing) x (mean within-cell rank among expressers) -- it keeps
# the prevalence signal that separates specific from ubiquitous axes (a gene
# expressed in one type scores ~0 elsewhere) while adding REO intensity. Averaging
# rank over expressers ONLY would discard prevalence and inflate ubiquitous axes.
.mean_rank_by_type <- function(rank_vec, logic_vec, labels, cell_types) {
  active <- as.logical(logic_vec)
  fl <- factor(labels, levels = cell_types)
  num <- tapply(ifelse(active, as.numeric(rank_vec), 0), fl, sum, na.rm = TRUE)
  den <- as.numeric(table(fl))   # all cells of the type, not just expressers
  out <- as.numeric(num) / pmax(den, 1)
  out[!is.finite(out)] <- 0
  names(out) <- cell_types
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

    data.frame(
      sender_type = sub$sender_type[1],
      receiver_type = sub$receiver_type[1],
      n_edges = as.integer(sub$n_edges[1]),
      edge_weight_sum = as.numeric(sub$edge_weight_sum[1]),
      n_scored_lr = nrow(finite_sub),
      n_active_lr = nrow(active_sub),
      active_lr_event_count = nrow(active_sub),
      sum_lcs = sum(active_sub$lcs, na.rm = TRUE),
      sum_lcs_all = sum(finite_sub$lcs, na.rm = TRUE),
      mean_lcs_active = if (nrow(active_sub) > 0) mean(active_sub$lcs, na.rm = TRUE) else NA_real_,
      mean_lcs_all = if (nrow(finite_sub) > 0) mean(finite_sub$lcs, na.rm = TRUE) else NA_real_,
      sum_active_edges = sum(active_sub$n_active_edges, na.rm = TRUE),
      mean_edge_support_fraction_active = if (length(support_frac) > 0) mean(support_frac, na.rm = TRUE) else NA_real_,
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
    top_lr <- sub$lr_pair[which.max(sub$lcs)]
    data.frame(
      sender_type = sub$sender_type[1],
      receiver_type = sub$receiver_type[1],
      pathway = sub$pathway[1],
      n_active_lr = nrow(sub),
      sum_lcs = sum(sub$lcs, na.rm = TRUE),
      sum_lcs_all = sum(sub$lcs, na.rm = TRUE),
      mean_lcs_active = mean(sub$lcs, na.rm = TRUE),
      sum_active_edges = sum(sub$n_active_edges, na.rm = TRUE),
      mean_edge_support_fraction_active = mean(support_frac, na.rm = TRUE),
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

# Transitional deprecation shim. LogicComm scores communication at the cell-type
# level; the former per-cell KNN/SNN neighborhood graph (which lives in
# expression space, not physical space) has been removed. Legacy neighborhood
# arguments are accepted via ... so existing scripts do not hard-error, but they
# are ignored with a warning. This shim is removed once callers are migrated.
.deprecate_neighborhood_args <- function(dots, fn = "summarize_celltype_communication") {
  if (!length(dots)) return(invisible(NULL))
  legacy <- c("knn_mat", "seurat_obj", "graph_name", "remove_self_edges",
              "graph_symmetrize", "edge_weight_mode")
  hit <- intersect(names(dots), legacy)
  m <- dots[["mode"]]
  if (!is.null(m) && !as.character(m)[1] %in% c("global", "celltype", "auto")) {
    hit <- c("mode", hit)
  }
  if (length(hit)) {
    warning(fn, "(): argument(s) ", paste(hit, collapse = ", "),
            " are deprecated and ignored. LogicComm now scores communication at ",
            "the cell-type level / global co-expression (no per-cell neighborhood graph).",
            call. = FALSE)
  }
  invisible(NULL)
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
    lcs = numeric(), lcs_global = numeric(), lcs_unweighted = numeric(),
    n_edges = integer(), edge_weight_sum = numeric(), n_active_edges = numeric(),
    n_global_possible = numeric(), n_active_global = numeric(), sender_cell_count = integer(), receiver_cell_count = integer(), ligand_active_n_sender = numeric(),
    receptor_active_n_receiver = numeric(), ligand_active_frac_sender = numeric(), receptor_active_frac_receiver = numeric(),
    active = logical(), stringsAsFactors = FALSE
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
    sum_lcs = numeric(nrow(group_defs)),
    sum_lcs_all = numeric(nrow(group_defs)),
    mean_lcs_active = rep(NA_real_, nrow(group_defs)),
    mean_lcs_all = rep(NA_real_, nrow(group_defs)),
    sum_active_edges = numeric(nrow(group_defs)),
    mean_edge_support_fraction_active = rep(NA_real_, nrow(group_defs)),
    top_lr_pair = rep(NA_character_, nrow(group_defs)),
    top_pathway = rep(NA_character_, nrow(group_defs)),
    stringsAsFactors = FALSE
  )
}

.empty_pathway_summary <- function() {
  data.frame(
    sender_type = character(), receiver_type = character(), pathway = character(), n_active_lr = integer(),
    sum_lcs = numeric(), sum_lcs_all = numeric(), mean_lcs_active = numeric(), sum_active_edges = numeric(),
    mean_edge_support_fraction_active = numeric(), top_lr_pair = character(), stringsAsFactors = FALSE
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
