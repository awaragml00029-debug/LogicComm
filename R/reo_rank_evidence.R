# R/reo_rank_evidence.R

#' Calculate a REO Rank Score Matrix
#'
#' Returns within-cell rank percentiles for ligand-receptor genes. Unlike the
#' binary REO matrix, this preserves how far a gene sits above or below other
#' genes in the same cell, which is useful for RankComp-inspired weak-signal
#' and threshold-stability analyses.
#'
#' @param expr_mat Expression input accepted by \code{\link{calc_REO_matrix}}.
#' @param genes Optional genes to retain. If \code{NULL}, all ligand-receptor
#'   genes in \code{lr_db} are used.
#' @param lr_db Ligand-receptor database used when \code{genes = NULL}.
#' @param layer Seurat layer or slot.
#' @param chunk_size Number of cells processed per chunk.
#' @param anchor_genes Optional anchor gene universe for rank calculation.
#' @param verbose Whether to print progress messages.
#' @return Dense numeric matrix of rank percentiles, genes x cells.
#' @examples
#' expr <- matrix(
#'   c(5, 1, 4, 2, 1, 5, 3, 4, 4, 2, 5, 1),
#'   nrow = 3,
#'   dimnames = list(c("L1", "R1", "T1"), paste0("cell", 1:4))
#' )
#' reo <- expr >= 3
#' rank_mat <- apply(expr, 2, rank) / nrow(expr)
#' lr_db <- data.frame(
#'   ligand = "L1",
#'   receptor = "R1",
#'   pathway = "toy",
#'   stringsAsFactors = FALSE
#' )
#' lcs <- data.frame(
#'   ligand = "L1",
#'   receptor = "R1",
#'   pathway = "toy",
#'   sample = c("S1", "S2"),
#'   group = c("control", "case"),
#'   sender = "A",
#'   receiver = "B",
#'   celltype_sender = "A",
#'   celltype_receiver = "B",
#'   LCS = c(0.2, 0.5),
#'   lcs = c(0.2, 0.5),
#'   mean_lcs = c(0.2, 0.5),
#'   delta_lcs = c(0.0, 0.3),
#'   p_value = c(0.5, 0.01),
#'   p_adj = c(0.5, 0.02),
#'   fdr = c(0.5, 0.02),
#'   stringsAsFactors = FALSE
#' )
#' sample_ct_list <- list(S1 = lcs, S2 = lcs)
#' group_info <- c(S1 = "control", S2 = "case")
#' knn <- matrix(1, nrow = 4, ncol = 4, dimnames = list(colnames(expr), colnames(expr)))
#' diag(knn) <- 0
#' toy_args <- list(
#'   x = lcs, result = lcs, results = lcs, lcs_df = lcs, ct_comm = lcs,
#'   comm_df = lcs, communication = lcs, celltype_comm = lcs,
#'   celltype_results = lcs, differential_results = lcs, diff_comm = lcs,
#'   glm_result = lcs, role_df = lcs, roles = lcs, specificity = lcs,
#'   null_pair = list(observed = lcs, null = lcs), reo_mat = reo,
#'   rank_mat = rank_mat, expr_mat = expr, expression = expr,
#'   lr_db = lr_db, samples = list(S1 = expr, S2 = expr),
#'   sample_ct_list = sample_ct_list, group_info = group_info,
#'   group_labels = group_info, groups = group_info, knn_mat = knn,
#'   output_dir = tempfile("logiccomm"), file = tempfile(fileext = ".csv"),
#'   path = tempfile(fileext = ".csv")
#' )
#' fun <- get("calc_REO_rank_score_matrix")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
calc_REO_rank_score_matrix <- function(expr_mat,
                                       genes = NULL,
                                       lr_db = lr_pairs_human,
                                       layer = "counts",
                                       chunk_size = 5000,
                                       anchor_genes = NULL,
                                       verbose = TRUE) {
  if (is.null(genes)) {
    if (is.null(lr_db)) stop("Either genes or lr_db must be provided.")
    genes <- all_lr_genes(lr_db)
  }
  genes <- unique(as.character(genes))
  genes <- genes[nzchar(genes) & !is.na(genes)]
  if (length(genes) == 0) stop("genes must contain at least one non-empty gene name.")

  res <- calc_REO_matrix(
    expr_mat,
    lr_genes = genes,
    rank_threshold = 0.5,
    layer = layer,
    chunk_size = chunk_size,
    return_rank = TRUE,
    rank_output = "percentile",
    anchor_genes = anchor_genes,
    verbose = verbose
  )
  rank_mat <- res$rank
  attr(rank_mat, "logic_type") <- "REO_rank_percentile"
  attr(rank_mat, "anchor_n_genes") <- attr(res$logic, "anchor_n_genes")
  attr(rank_mat, "anchor_source") <- attr(res$logic, "anchor_source")
  class(rank_mat) <- c("LogicCommRankMatrix", class(rank_mat))
  rank_mat
}

#' Score Per-Cell Ligand-Receptor Rank Activity
#'
#' Converts a gene-level REO rank matrix into continuous sender and receiver
#' rank scores for each ligand-receptor pair.
#'
#' @param rank_mat Rank-percentile matrix, or a \code{LogicCommREOResult} with
#'   a non-null \code{rank} element.
#' @param lr_db Ligand-receptor database.
#' @param mode One of \code{"sender"}, \code{"receiver"}, or \code{"both"}.
#' @param complex_aggregate How to combine multi-subunit ligand/receptor ranks.
#'   \code{"min"} is conservative and is the default.
#' @param aggregate Whether to return per-cell aggregate scores.
#' @param verbose Whether to print progress messages.
#' @return List with LR-pair x cell sender/receiver rank matrices and optional
#'   aggregate per-cell rank scores.
#' @examples
#' expr <- matrix(
#'   c(5, 1, 4, 2, 1, 5, 3, 4, 4, 2, 5, 1),
#'   nrow = 3,
#'   dimnames = list(c("L1", "R1", "T1"), paste0("cell", 1:4))
#' )
#' reo <- expr >= 3
#' rank_mat <- apply(expr, 2, rank) / nrow(expr)
#' lr_db <- data.frame(
#'   ligand = "L1",
#'   receptor = "R1",
#'   pathway = "toy",
#'   stringsAsFactors = FALSE
#' )
#' lcs <- data.frame(
#'   ligand = "L1",
#'   receptor = "R1",
#'   pathway = "toy",
#'   sample = c("S1", "S2"),
#'   group = c("control", "case"),
#'   sender = "A",
#'   receiver = "B",
#'   celltype_sender = "A",
#'   celltype_receiver = "B",
#'   LCS = c(0.2, 0.5),
#'   lcs = c(0.2, 0.5),
#'   mean_lcs = c(0.2, 0.5),
#'   delta_lcs = c(0.0, 0.3),
#'   p_value = c(0.5, 0.01),
#'   p_adj = c(0.5, 0.02),
#'   fdr = c(0.5, 0.02),
#'   stringsAsFactors = FALSE
#' )
#' sample_ct_list <- list(S1 = lcs, S2 = lcs)
#' group_info <- c(S1 = "control", S2 = "case")
#' knn <- matrix(1, nrow = 4, ncol = 4, dimnames = list(colnames(expr), colnames(expr)))
#' diag(knn) <- 0
#' toy_args <- list(
#'   x = lcs, result = lcs, results = lcs, lcs_df = lcs, ct_comm = lcs,
#'   comm_df = lcs, communication = lcs, celltype_comm = lcs,
#'   celltype_results = lcs, differential_results = lcs, diff_comm = lcs,
#'   glm_result = lcs, role_df = lcs, roles = lcs, specificity = lcs,
#'   null_pair = list(observed = lcs, null = lcs), reo_mat = reo,
#'   rank_mat = rank_mat, expr_mat = expr, expression = expr,
#'   lr_db = lr_db, samples = list(S1 = expr, S2 = expr),
#'   sample_ct_list = sample_ct_list, group_info = group_info,
#'   group_labels = group_info, groups = group_info, knn_mat = knn,
#'   output_dir = tempfile("logiccomm"), file = tempfile(fileext = ".csv"),
#'   path = tempfile(fileext = ".csv")
#' )
#' fun <- get("score_lr_rank_activity")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
score_lr_rank_activity <- function(rank_mat,
                                   lr_db = lr_pairs_human,
                                   mode = c("both", "sender", "receiver"),
                                   complex_aggregate = c("min", "mean", "max"),
                                   aggregate = TRUE,
                                   verbose = TRUE) {
  mode <- match.arg(mode)
  complex_aggregate <- match.arg(complex_aggregate)
  rank_mat <- .as_logiccomm_rank_matrix(rank_mat)

  n_pairs <- nrow(lr_db)
  n_cells <- ncol(rank_mat)
  cell_names <- colnames(rank_mat)
  if (verbose) {
    message(sprintf("[RankScore] Scoring %d L-R pairs for %d cells...", n_pairs, n_cells))
  }

  sender_mat <- receiver_mat <- NULL
  if (mode %in% c("sender", "both")) {
    sender_mat <- matrix(NA_real_, nrow = n_pairs, ncol = n_cells,
                         dimnames = list(lr_db$lr_pair, cell_names))
  }
  if (mode %in% c("receiver", "both")) {
    receiver_mat <- matrix(NA_real_, nrow = n_pairs, ncol = n_cells,
                           dimnames = list(lr_db$lr_pair, cell_names))
  }

  for (i in seq_len(n_pairs)) {
    if (!is.null(sender_mat)) {
      sender_mat[i, ] <- .resolve_complex_rank(
        .lr_list_entry(lr_db, "ligand_genes", i),
        rank_mat,
        aggregate = complex_aggregate
      )
    }
    if (!is.null(receiver_mat)) {
      receiver_mat[i, ] <- .resolve_complex_rank(
        .lr_list_entry(lr_db, "receptor_genes", i),
        rank_mat,
        aggregate = complex_aggregate
      )
    }
  }

  out <- list()
  if (!is.null(sender_mat)) {
    out$sender_rank_mat <- sender_mat
    out$sender_rank_score <- .col_means_or_na(sender_mat)
  }
  if (!is.null(receiver_mat)) {
    out$receiver_rank_mat <- receiver_mat
    out$receiver_rank_score <- .col_means_or_na(receiver_mat)
  }
  if (isTRUE(aggregate) && mode == "both") {
    combined <- pmax(sender_mat, receiver_mat, na.rm = TRUE)
    combined[!is.finite(combined)] <- NA_real_
    out$comm_rank_score <- .col_means_or_na(combined)
  }

  if (verbose) message("[RankScore] Done.")
  out
}

#' Identify Rank-Aware Logic Consensus Evidence
#'
#' Computes ligand-receptor evidence that keeps the original binary LCS but adds
#' continuous rank-dominance, rank-margin, threshold-stability, and optional
#' cell-label specificity scores.
#'
#' @param reo_mat Optional binary REO matrix or \code{LogicCommREOResult}.
#' @param rank_mat Optional rank-percentile matrix. Required unless supplied via
#'   \code{reo_mat} or calculated from \code{expr_mat}.
#' @param expr_mat Optional expression input used to calculate \code{rank_mat}
#'   when no rank matrix is supplied.
#' @param seurat_obj Optional Seurat object for KNN graph extraction.
#' @param knn_mat Optional cells x cells adjacency matrix.
#' @param lr_db Ligand-receptor database.
#' @param graph_name Optional Seurat graph name.
#' @param rank_threshold Binary activity threshold applied to ranks when
#'   \code{reo_mat} is not supplied.
#' @param threshold_grid Rank thresholds used for stability scoring.
#' @param cell_labels Optional named cell-type labels for specificity scoring.
#' @param layer Seurat layer used when \code{expr_mat} is provided.
#' @param chunk_size Chunk size used when \code{expr_mat} is provided.
#' @param anchor_genes Optional anchor genes used when \code{expr_mat} is provided.
#' @param remove_self_edges Whether to remove diagonal graph edges.
#' @param graph_symmetrize Graph symmetrization mode.
#' @param edge_weight_mode \code{"binary"} or \code{"weighted"} graph scoring.
#' @param complex_aggregate How to combine multi-subunit ranks.
#' @param lcs_threshold Activity threshold used only for evidence-tier labels.
#' @param verbose Whether to print progress messages.
#' @return Data frame with one row per LR pair and rank-aware evidence columns.
#' @examples
#' expr <- matrix(
#'   c(5, 1, 4, 2, 1, 5, 3, 4, 4, 2, 5, 1),
#'   nrow = 3,
#'   dimnames = list(c("L1", "R1", "T1"), paste0("cell", 1:4))
#' )
#' reo <- expr >= 3
#' rank_mat <- apply(expr, 2, rank) / nrow(expr)
#' lr_db <- data.frame(
#'   ligand = "L1",
#'   receptor = "R1",
#'   pathway = "toy",
#'   stringsAsFactors = FALSE
#' )
#' lcs <- data.frame(
#'   ligand = "L1",
#'   receptor = "R1",
#'   pathway = "toy",
#'   sample = c("S1", "S2"),
#'   group = c("control", "case"),
#'   sender = "A",
#'   receiver = "B",
#'   celltype_sender = "A",
#'   celltype_receiver = "B",
#'   LCS = c(0.2, 0.5),
#'   lcs = c(0.2, 0.5),
#'   mean_lcs = c(0.2, 0.5),
#'   delta_lcs = c(0.0, 0.3),
#'   p_value = c(0.5, 0.01),
#'   p_adj = c(0.5, 0.02),
#'   fdr = c(0.5, 0.02),
#'   stringsAsFactors = FALSE
#' )
#' sample_ct_list <- list(S1 = lcs, S2 = lcs)
#' group_info <- c(S1 = "control", S2 = "case")
#' knn <- matrix(1, nrow = 4, ncol = 4, dimnames = list(colnames(expr), colnames(expr)))
#' diag(knn) <- 0
#' toy_args <- list(
#'   x = lcs, result = lcs, results = lcs, lcs_df = lcs, ct_comm = lcs,
#'   comm_df = lcs, communication = lcs, celltype_comm = lcs,
#'   celltype_results = lcs, differential_results = lcs, diff_comm = lcs,
#'   glm_result = lcs, role_df = lcs, roles = lcs, specificity = lcs,
#'   null_pair = list(observed = lcs, null = lcs), reo_mat = reo,
#'   rank_mat = rank_mat, expr_mat = expr, expression = expr,
#'   lr_db = lr_db, samples = list(S1 = expr, S2 = expr),
#'   sample_ct_list = sample_ct_list, group_info = group_info,
#'   group_labels = group_info, groups = group_info, knn_mat = knn,
#'   output_dir = tempfile("logiccomm"), file = tempfile(fileext = ".csv"),
#'   path = tempfile(fileext = ".csv")
#' )
#' fun <- get("IdentifyRankLogicConsensus")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
IdentifyRankLogicConsensus <- function(reo_mat = NULL,
                                       rank_mat = NULL,
                                       expr_mat = NULL,
                                       seurat_obj = NULL,
                                       knn_mat = NULL,
                                       lr_db = lr_pairs_human,
                                       graph_name = NULL,
                                       rank_threshold = 0.5,
                                       threshold_grid = c(0.4, 0.5, 0.6, 0.7),
                                       cell_labels = NULL,
                                       layer = "counts",
                                       chunk_size = 5000,
                                       anchor_genes = NULL,
                                       remove_self_edges = TRUE,
                                       graph_symmetrize = c("or", "none", "max"),
                                       edge_weight_mode = c("binary", "weighted"),
                                       complex_aggregate = c("min", "mean", "max"),
                                       lcs_threshold = 0.01,
                                       verbose = TRUE) {
  stopifnot(is.numeric(rank_threshold), length(rank_threshold) == 1,
            rank_threshold > 0, rank_threshold < 1)
  threshold_grid <- sort(unique(as.numeric(threshold_grid)))
  threshold_grid <- threshold_grid[is.finite(threshold_grid) & threshold_grid > 0 & threshold_grid < 1]
  if (length(threshold_grid) == 0) stop("threshold_grid must contain values between 0 and 1.")
  edge_weight_mode <- match.arg(edge_weight_mode)
  complex_aggregate <- match.arg(complex_aggregate)
  if (missing(graph_symmetrize)) {
    graph_symmetrize <- if (edge_weight_mode == "binary" && isTRUE(remove_self_edges)) "or" else "none"
  } else {
    graph_symmetrize <- match.arg(graph_symmetrize)
  }

  parsed <- .resolve_rank_logic_inputs(
    reo_mat = reo_mat,
    rank_mat = rank_mat,
    expr_mat = expr_mat,
    lr_db = lr_db,
    rank_threshold = rank_threshold,
    layer = layer,
    chunk_size = chunk_size,
    anchor_genes = anchor_genes,
    verbose = verbose
  )
  logic_mat <- parsed$logic
  rank_mat <- parsed$rank
  cell_names <- colnames(rank_mat)

  graph <- .resolve_rank_graph(
    seurat_obj = seurat_obj,
    knn_mat = knn_mat,
    graph_name = graph_name,
    cell_names = cell_names,
    remove_self_edges = remove_self_edges,
    graph_symmetrize = graph_symmetrize,
    edge_weight_mode = edge_weight_mode,
    verbose = verbose
  )
  use_knn <- graph$use_knn

  labels <- .align_rank_cell_labels(cell_labels, cell_names)
  if (verbose) {
    mode_msg <- if (use_knn) sprintf("KNN mode with %d edges", length(graph$edges$i)) else "global mode"
    message(sprintf("[RankLCS] Scoring %d L-R pairs (%s).", nrow(lr_db), mode_msg))
  }

  rows <- vector("list", nrow(lr_db))
  for (i in seq_len(nrow(lr_db))) {
    rows[[i]] <- .score_one_rank_lr_pair(
      i = i,
      lr_db = lr_db,
      logic_mat = logic_mat,
      rank_mat = rank_mat,
      use_knn = use_knn,
      edges = graph$edges,
      rank_threshold = rank_threshold,
      threshold_grid = threshold_grid,
      labels = labels,
      complex_aggregate = complex_aggregate,
      lcs_threshold = lcs_threshold
    )
  }

  out <- do.call(rbind, rows)
  out$rank <- rank(-out$specificity_weighted_rank_lcs, ties.method = "min", na.last = "keep")
  out <- out[order(out$rank, out$lr_pair, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("LogicCommRankEvidence", "data.frame")
  attr(out, "rank_threshold") <- rank_threshold
  attr(out, "threshold_grid") <- threshold_grid
  attr(out, "mode") <- if (use_knn) "knn" else "global"
  if (verbose) {
    n_scored <- sum(is.finite(out$rank_dominance_lcs))
    message(sprintf("[RankLCS] Done. %d / %d pairs scored.", n_scored, nrow(out)))
  }
  out
}

#' Compare Rank-Aware Logic Evidence Between Groups
#'
#' Compares continuous rank-aware LR evidence across biological samples. This is
#' intended to complement \code{\link{CompareLogicGroups}} when binary active
#' frequency is saturated but rank strength still differs between conditions.
#'
#' @param rank_result_list Named list of data frames from
#'   \code{\link{IdentifyRankLogicConsensus}}.
#' @param group_info Named vector mapping samples to group labels.
#' @param score_col Rank evidence column to compare.
#' @param case_label Case group label.
#' @param ctrl_label Control group label.
#' @param rank_score_threshold Threshold for calling a sample rank-positive.
#' @param min_samples_per_group Minimum samples required per group.
#' @param lr_db Optional LR database for metadata.
#' @param verbose Whether to print progress messages.
#' @return Data frame sorted by case-control rank evidence difference.
#' @examples
#' expr <- matrix(
#'   c(5, 1, 4, 2, 1, 5, 3, 4, 4, 2, 5, 1),
#'   nrow = 3,
#'   dimnames = list(c("L1", "R1", "T1"), paste0("cell", 1:4))
#' )
#' reo <- expr >= 3
#' rank_mat <- apply(expr, 2, rank) / nrow(expr)
#' lr_db <- data.frame(
#'   ligand = "L1",
#'   receptor = "R1",
#'   pathway = "toy",
#'   stringsAsFactors = FALSE
#' )
#' lcs <- data.frame(
#'   ligand = "L1",
#'   receptor = "R1",
#'   pathway = "toy",
#'   sample = c("S1", "S2"),
#'   group = c("control", "case"),
#'   sender = "A",
#'   receiver = "B",
#'   celltype_sender = "A",
#'   celltype_receiver = "B",
#'   LCS = c(0.2, 0.5),
#'   lcs = c(0.2, 0.5),
#'   mean_lcs = c(0.2, 0.5),
#'   delta_lcs = c(0.0, 0.3),
#'   p_value = c(0.5, 0.01),
#'   p_adj = c(0.5, 0.02),
#'   fdr = c(0.5, 0.02),
#'   stringsAsFactors = FALSE
#' )
#' sample_ct_list <- list(S1 = lcs, S2 = lcs)
#' group_info <- c(S1 = "control", S2 = "case")
#' knn <- matrix(1, nrow = 4, ncol = 4, dimnames = list(colnames(expr), colnames(expr)))
#' diag(knn) <- 0
#' toy_args <- list(
#'   x = lcs, result = lcs, results = lcs, lcs_df = lcs, ct_comm = lcs,
#'   comm_df = lcs, communication = lcs, celltype_comm = lcs,
#'   celltype_results = lcs, differential_results = lcs, diff_comm = lcs,
#'   glm_result = lcs, role_df = lcs, roles = lcs, specificity = lcs,
#'   null_pair = list(observed = lcs, null = lcs), reo_mat = reo,
#'   rank_mat = rank_mat, expr_mat = expr, expression = expr,
#'   lr_db = lr_db, samples = list(S1 = expr, S2 = expr),
#'   sample_ct_list = sample_ct_list, group_info = group_info,
#'   group_labels = group_info, groups = group_info, knn_mat = knn,
#'   output_dir = tempfile("logiccomm"), file = tempfile(fileext = ".csv"),
#'   path = tempfile(fileext = ".csv")
#' )
#' fun <- get("CompareRankLogicGroups")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
CompareRankLogicGroups <- function(rank_result_list,
                                   group_info,
                                   score_col = "specificity_weighted_rank_lcs",
                                   case_label = "Case",
                                   ctrl_label = "Ctrl",
                                   rank_score_threshold = 0.05,
                                   min_samples_per_group = 2,
                                   lr_db = NULL,
                                   verbose = TRUE) {
  if (is.null(names(rank_result_list))) {
    names(rank_result_list) <- paste0("S", seq_along(rank_result_list))
  }
  if (is.null(names(group_info))) names(group_info) <- names(rank_result_list)
  group_info <- group_info[names(rank_result_list)]
  if (anyNA(group_info)) stop("group_info must provide labels for all rank_result_list samples.")

  case_s <- names(group_info)[group_info == case_label]
  ctrl_s <- names(group_info)[group_info == ctrl_label]
  if (length(case_s) < min_samples_per_group || length(ctrl_s) < min_samples_per_group) {
    stop("Not enough samples per group for CompareRankLogicGroups().")
  }

  missing_col <- vapply(rank_result_list, function(x) !score_col %in% names(x), logical(1))
  if (any(missing_col)) {
    stop("score_col '", score_col, "' missing from: ",
         paste(names(rank_result_list)[missing_col], collapse = ", "))
  }

  all_pairs <- sort(unique(unlist(lapply(rank_result_list, function(x) as.character(x$lr_pair)),
                                use.names = FALSE)))
  score_mat <- matrix(NA_real_, nrow = length(all_pairs), ncol = length(rank_result_list),
                      dimnames = list(all_pairs, names(rank_result_list)))
  binary_mat <- score_mat

  for (s in names(rank_result_list)) {
    df <- as.data.frame(rank_result_list[[s]], stringsAsFactors = FALSE)
    idx <- match(df$lr_pair, all_pairs)
    score_mat[idx, s] <- as.numeric(df[[score_col]])
    if ("binary_lcs" %in% names(df)) binary_mat[idx, s] <- as.numeric(df$binary_lcs)
  }

  case_mat <- score_mat[, case_s, drop = FALSE]
  ctrl_mat <- score_mat[, ctrl_s, drop = FALSE]
  case_mean <- rowMeans(case_mat, na.rm = TRUE)
  ctrl_mean <- rowMeans(ctrl_mat, na.rm = TRUE)
  case_mean[!is.finite(case_mean)] <- NA_real_
  ctrl_mean[!is.finite(ctrl_mean)] <- NA_real_

  p_wilcox <- vapply(seq_along(all_pairs), function(i) {
    x <- case_mat[i, ]; y <- ctrl_mat[i, ]
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    if (length(x) < 2 || length(y) < 2) return(1)
    tryCatch(stats::wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) 1)
  }, numeric(1))

  out <- data.frame(
    lr_pair = all_pairs,
    case_mean_rank_score = case_mean,
    ctrl_mean_rank_score = ctrl_mean,
    delta_rank_score = case_mean - ctrl_mean,
    log2fc_rank_score = log2((case_mean + 1e-6) / (ctrl_mean + 1e-6)),
    case_freq_rank_positive = rowMeans(case_mat >= rank_score_threshold, na.rm = TRUE),
    ctrl_freq_rank_positive = rowMeans(ctrl_mat >= rank_score_threshold, na.rm = TRUE),
    p_wilcox = p_wilcox,
    fdr_wilcox = stats::p.adjust(p_wilcox, method = "BH"),
    n_case_avail = rowSums(is.finite(case_mat)),
    n_ctrl_avail = rowSums(is.finite(ctrl_mat)),
    stringsAsFactors = FALSE
  )
  out$rank_positive_asymmetry <- out$case_freq_rank_positive - out$ctrl_freq_rank_positive

  case_bin <- binary_mat[, case_s, drop = FALSE]
  ctrl_bin <- binary_mat[, ctrl_s, drop = FALSE]
  out$case_mean_binary_lcs <- rowMeans(case_bin, na.rm = TRUE)
  out$ctrl_mean_binary_lcs <- rowMeans(ctrl_bin, na.rm = TRUE)
  out$case_mean_binary_lcs[!is.finite(out$case_mean_binary_lcs)] <- NA_real_
  out$ctrl_mean_binary_lcs[!is.finite(out$ctrl_mean_binary_lcs)] <- NA_real_

  out <- .append_rank_lr_metadata(out, rank_result_list, lr_db)
  out <- out[order(-out$delta_rank_score, out$fdr_wilcox, out$lr_pair, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("LogicCommRankGroupComparison", "data.frame")
  attr(out, "score_col") <- score_col
  attr(out, "case_label") <- case_label
  attr(out, "ctrl_label") <- ctrl_label
  attr(out, "score_mat") <- score_mat
  if (verbose) {
    message(sprintf("[RankCompare] Compared %d LR pairs across %d samples.", nrow(out), ncol(score_mat)))
  }
  out
}

.as_logiccomm_rank_matrix <- function(rank_mat) {
  if (is.list(rank_mat) && !is.null(rank_mat$rank)) rank_mat <- rank_mat$rank
  if (is.null(rank_mat) || is.null(dim(rank_mat))) stop("rank_mat must be a genes x cells matrix.")
  if (is.null(rownames(rank_mat)) || is.null(colnames(rank_mat))) {
    stop("rank_mat must have gene rownames and cell colnames.")
  }
  rank_mat <- as.matrix(rank_mat)
  storage.mode(rank_mat) <- "double"
  rank_mat
}

.resolve_rank_logic_inputs <- function(reo_mat,
                                       rank_mat,
                                       expr_mat,
                                       lr_db,
                                       rank_threshold,
                                       layer,
                                       chunk_size,
                                       anchor_genes,
                                       verbose) {
  logic_mat <- NULL
  if (inherits(reo_mat, "LogicCommREOResult")) {
    logic_mat <- reo_mat$logic
    if (is.null(rank_mat)) rank_mat <- reo_mat$rank
  } else if (!is.null(reo_mat)) {
    logic_mat <- reo_mat
  }

  if (is.null(rank_mat)) {
    if (is.null(expr_mat)) {
      stop("rank_mat is required unless reo_mat is a LogicCommREOResult with ranks or expr_mat is supplied.")
    }
    rank_mat <- calc_REO_rank_score_matrix(
      expr_mat,
      genes = all_lr_genes(lr_db),
      lr_db = lr_db,
      layer = layer,
      chunk_size = chunk_size,
      anchor_genes = anchor_genes,
      verbose = verbose
    )
  }
  rank_mat <- .as_logiccomm_rank_matrix(rank_mat)

  if (is.null(logic_mat)) {
    logic_mat <- Matrix::Matrix((rank_mat >= rank_threshold) * 1L, sparse = TRUE)
    rownames(logic_mat) <- rownames(rank_mat)
    colnames(logic_mat) <- colnames(rank_mat)
  } else {
    if (is.null(rownames(logic_mat)) || is.null(colnames(logic_mat))) {
      stop("reo_mat must have gene rownames and cell colnames.")
    }
    if (!all(colnames(logic_mat) %in% colnames(rank_mat))) {
      stop("rank_mat must contain all reo_mat cells.")
    }
    if (!all(rownames(logic_mat) %in% rownames(rank_mat))) {
      missing_genes <- setdiff(rownames(logic_mat), rownames(rank_mat))
      stop("rank_mat must contain all reo_mat genes. Missing: ",
           paste(utils::head(missing_genes, 10), collapse = ", "))
    }
    rank_mat <- rank_mat[rownames(logic_mat), colnames(logic_mat), drop = FALSE]
    if (!inherits(logic_mat, "sparseMatrix")) {
      logic_mat <- Matrix::Matrix(logic_mat, sparse = TRUE)
    }
  }

  list(logic = logic_mat, rank = rank_mat)
}

.resolve_rank_graph <- function(seurat_obj,
                                knn_mat,
                                graph_name,
                                cell_names,
                                remove_self_edges,
                                graph_symmetrize,
                                edge_weight_mode,
                                verbose) {
  use_knn <- FALSE
  edges <- list(i = integer(0), j = integer(0), w = numeric(0))
  if (!is.null(knn_mat)) {
    knn_mat <- .validate_and_align_knn(knn_mat, cell_names)
    knn_mat <- .symmetrize_sparse_graph(knn_mat, graph_symmetrize)
    edges <- .sparse_to_edges(knn_mat, remove_self_edges = remove_self_edges,
                              edge_weight_mode = edge_weight_mode)
    use_knn <- length(edges$i) > 0
  } else if (!is.null(seurat_obj)) {
    tryCatch({
      knn_mat <- .extract_knn(seurat_obj, graph_name)
      knn_mat <- .validate_and_align_knn(knn_mat, cell_names)
      knn_mat <- .symmetrize_sparse_graph(knn_mat, graph_symmetrize)
      edges <<- .sparse_to_edges(knn_mat, remove_self_edges = remove_self_edges,
                                 edge_weight_mode = edge_weight_mode)
      use_knn <<- length(edges$i) > 0
    }, error = function(e) {
      if (verbose) message("[RankLCS] KNN graph not found; using global mode.")
    })
  }
  list(use_knn = use_knn, edges = edges)
}

.score_one_rank_lr_pair <- function(i,
                                    lr_db,
                                    logic_mat,
                                    rank_mat,
                                    use_knn,
                                    edges,
                                    rank_threshold,
                                    threshold_grid,
                                    labels,
                                    complex_aggregate,
                                    lcs_threshold) {
  lig_genes <- .lr_list_entry(lr_db, "ligand_genes", i)
  rec_genes <- .lr_list_entry(lr_db, "receptor_genes", i)
  lig_found <- length(intersect(lig_genes, rownames(rank_mat))) == length(lig_genes) && length(lig_genes) > 0
  rec_found <- length(intersect(rec_genes, rownames(rank_mat))) == length(rec_genes) && length(rec_genes) > 0

  base <- .empty_rank_evidence_row(lr_db, i, lig_found, rec_found, rank_threshold)
  if (!lig_found || !rec_found) return(base)

  lig_rank <- .resolve_complex_rank(lig_genes, rank_mat, aggregate = complex_aggregate)
  rec_rank <- .resolve_complex_rank(rec_genes, rank_mat, aggregate = complex_aggregate)
  lig_logic <- .resolve_complex_logic(lig_genes, logic_mat)
  rec_logic <- .resolve_complex_logic(rec_genes, logic_mat)

  if (use_knn) {
    idx_i <- edges$i
    idx_j <- edges$j
    weights <- edges$w
    pair_rank <- pmin(lig_rank[idx_i], rec_rank[idx_j])
    active <- lig_logic[idx_i] & rec_logic[idx_j]
    threshold_scores <- vapply(threshold_grid, function(thr) {
      .weighted_mean01(lig_rank[idx_i] >= thr & rec_rank[idx_j] >= thr, weights)
    }, numeric(1))
    spec <- .rank_event_specificity(
      active = active,
      rank_value = pair_rank,
      edge_i = idx_i,
      edge_j = idx_j,
      labels = labels,
      weights = weights
    )
  } else {
    weights <- rep(1, length(lig_rank))
    pair_rank <- pmin(lig_rank, rec_rank)
    active <- lig_logic & rec_logic
    threshold_scores <- vapply(threshold_grid, function(thr) {
      .weighted_mean01(lig_rank >= thr & rec_rank >= thr, weights)
    }, numeric(1))
    spec <- .rank_cell_specificity(active, pair_rank, labels, weights)
  }

  binary_lcs <- .weighted_mean01(active, weights)
  rank_dominance_lcs <- .weighted_mean_numeric(pair_rank, weights)
  rank_margin_lcs <- .weighted_mean_numeric(pmax(0, pair_rank - rank_threshold) / (1 - rank_threshold), weights)
  threshold_stability_lcs <- mean(threshold_scores, na.rm = TRUE)
  if (!is.finite(threshold_stability_lcs)) threshold_stability_lcs <- NA_real_
  active_rank_dominance <- .weighted_mean_numeric(pair_rank[active], weights[active])
  n_events <- sum(is.finite(pair_rank))
  n_active_events <- sum(active & is.finite(pair_rank), na.rm = TRUE)
  specificity_score <- spec$specificity_score
  specificity_multiplier <- if (is.finite(specificity_score)) 0.5 + 0.5 * specificity_score else 1
  rank_evidence_score <- 0.5 * rank_dominance_lcs + 0.3 * threshold_stability_lcs + 0.2 * rank_margin_lcs
  specificity_weighted_rank_lcs <- rank_evidence_score * specificity_multiplier

  base$binary_lcs <- binary_lcs
  base$rank_dominance_lcs <- rank_dominance_lcs
  base$active_rank_dominance <- active_rank_dominance
  base$rank_margin_lcs <- rank_margin_lcs
  base$threshold_stability_lcs <- threshold_stability_lcs
  base$rank_evidence_score <- rank_evidence_score
  base$specificity_score <- specificity_score
  base$specificity_weighted_rank_lcs <- specificity_weighted_rank_lcs
  base$n_events <- n_events
  base$n_active_events <- n_active_events
  base$top_specificity_group <- spec$top_group
  base$n_specificity_groups <- spec$n_groups
  base$evidence_tier <- .rank_evidence_tier(
    binary_lcs,
    rank_margin_lcs,
    threshold_stability_lcs,
    specificity_score,
    lcs_threshold
  )
  base
}

.empty_rank_evidence_row <- function(lr_db, i, lig_found, rec_found, rank_threshold) {
  data.frame(
    lr_pair = as.character(lr_db$lr_pair[i]),
    ligand = if ("ligand" %in% names(lr_db)) as.character(lr_db$ligand[i]) else NA_character_,
    receptor = if ("receptor" %in% names(lr_db)) as.character(lr_db$receptor[i]) else NA_character_,
    pathway = if ("pathway" %in% names(lr_db)) as.character(lr_db$pathway[i]) else NA_character_,
    annotation = if ("annotation" %in% names(lr_db)) as.character(lr_db$annotation[i]) else NA_character_,
    ligand_found = lig_found,
    receptor_found = rec_found,
    rank_threshold = rank_threshold,
    binary_lcs = NA_real_,
    rank_dominance_lcs = NA_real_,
    active_rank_dominance = NA_real_,
    rank_margin_lcs = NA_real_,
    threshold_stability_lcs = NA_real_,
    rank_evidence_score = NA_real_,
    specificity_score = NA_real_,
    specificity_weighted_rank_lcs = NA_real_,
    n_events = 0L,
    n_active_events = 0L,
    top_specificity_group = NA_character_,
    n_specificity_groups = 0L,
    evidence_tier = "unscored_missing_genes",
    stringsAsFactors = FALSE
  )
}

.rank_evidence_tier <- function(binary_lcs,
                                rank_margin_lcs,
                                threshold_stability_lcs,
                                specificity_score,
                                lcs_threshold) {
  if (!is.finite(binary_lcs) || !is.finite(rank_margin_lcs) || !is.finite(threshold_stability_lcs)) {
    return("unscored")
  }
  specific_enough <- !is.finite(specificity_score) || specificity_score >= 0.5
  if (binary_lcs >= lcs_threshold && rank_margin_lcs >= 0.15 &&
      threshold_stability_lcs >= lcs_threshold && specific_enough) {
    return("Tier 1: strong, specific, stable")
  }
  if (binary_lcs >= lcs_threshold && rank_margin_lcs >= 0.05 &&
      threshold_stability_lcs >= lcs_threshold) {
    return("Tier 2: active rank-supported axis")
  }
  if (rank_margin_lcs >= 0.05 && specific_enough) {
    return("Tier 3: rank-supported weak candidate")
  }
  "Tier 4: weak or threshold-sensitive"
}

.rank_event_specificity <- function(active, rank_value, edge_i, edge_j, labels, weights) {
  if (is.null(labels)) return(list(specificity_score = NA_real_, top_group = NA_character_, n_groups = 0L))
  keep <- active & is.finite(rank_value) & !is.na(labels[edge_i]) & !is.na(labels[edge_j])
  if (!any(keep)) return(list(specificity_score = NA_real_, top_group = NA_character_, n_groups = 0L))
  group <- paste(labels[edge_i[keep]], labels[edge_j[keep]], sep = "|")
  contrib <- rank_value[keep] * weights[keep]
  .specificity_from_group_weights(group, contrib)
}

.rank_cell_specificity <- function(active, rank_value, labels, weights) {
  if (is.null(labels)) return(list(specificity_score = NA_real_, top_group = NA_character_, n_groups = 0L))
  keep <- active & is.finite(rank_value) & !is.na(labels)
  if (!any(keep)) return(list(specificity_score = NA_real_, top_group = NA_character_, n_groups = 0L))
  .specificity_from_group_weights(labels[keep], rank_value[keep] * weights[keep])
}

.specificity_from_group_weights <- function(group, weight) {
  group <- as.character(group)
  weight <- as.numeric(weight)
  keep <- nzchar(group) & !is.na(group) & is.finite(weight) & weight > 0
  if (!any(keep)) return(list(specificity_score = NA_real_, top_group = NA_character_, n_groups = 0L))
  totals <- tapply(weight[keep], group[keep], sum)
  totals <- totals[is.finite(totals) & totals > 0]
  if (length(totals) == 0) return(list(specificity_score = NA_real_, top_group = NA_character_, n_groups = 0L))
  p <- totals / sum(totals)
  if (length(p) == 1) {
    specificity <- 1
  } else {
    entropy <- -sum(p * log(p))
    specificity <- 1 - entropy / log(length(p))
  }
  list(
    specificity_score = as.numeric(max(0, min(1, specificity))),
    top_group = names(totals)[which.max(totals)],
    n_groups = length(totals)
  )
}

.weighted_mean01 <- function(x, weights = NULL) {
  x <- as.logical(x)
  .weighted_mean_numeric(as.numeric(x), weights)
}

.weighted_mean_numeric <- function(x, weights = NULL) {
  x <- as.numeric(x)
  if (is.null(weights)) weights <- rep(1, length(x))
  weights <- as.numeric(weights)
  keep <- is.finite(x) & is.finite(weights) & weights > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * weights[keep]) / sum(weights[keep])
}

.col_means_or_na <- function(x) {
  out <- colMeans(x, na.rm = TRUE)
  out[!is.finite(out)] <- NA_real_
  out
}

.align_rank_cell_labels <- function(cell_labels, cell_names) {
  if (is.null(cell_labels)) return(NULL)
  label_names <- names(cell_labels)
  labels <- as.character(cell_labels)
  if (!is.null(label_names)) {
    names(labels) <- label_names
    labels <- labels[cell_names]
  } else if (length(labels) != length(cell_names)) {
    stop("cell_labels must be named by cell, or have length equal to ncol(rank_mat).")
  }
  names(labels) <- cell_names
  labels
}

.append_rank_lr_metadata <- function(out, rank_result_list, lr_db = NULL) {
  meta_cols <- c("lr_pair", "ligand", "receptor", "pathway", "annotation")
  meta <- NULL
  if (!is.null(lr_db)) {
    meta <- as.data.frame(lr_db, stringsAsFactors = FALSE)
    meta <- meta[, intersect(meta_cols, names(meta)), drop = FALSE]
  } else {
    for (x in rank_result_list) {
      cols <- intersect(meta_cols, names(x))
      if (length(cols) > 1) {
        meta <- unique(as.data.frame(x[, cols, drop = FALSE], stringsAsFactors = FALSE))
        break
      }
    }
  }
  if (is.null(meta) || !"lr_pair" %in% names(meta)) return(out)
  meta <- meta[!duplicated(meta$lr_pair), , drop = FALSE]
  merge(out, meta, by = "lr_pair", all.x = TRUE, sort = FALSE)
}

