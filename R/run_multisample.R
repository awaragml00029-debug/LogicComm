# R/run_multisample.R

#' Run Full LogicComm Pipeline on Multiple Samples
#'
#' One-stop function for multi-sample cell-cell communication analysis.
#' Runs \code{\link{calc_REO_matrix}} and \code{\link{IdentifyLogicConsensus}}
#' for each sample, optionally in parallel, then calls
#' \code{\link{CompareLogicGroups}}.
#'
#' @details
#' Parallel processing is handled via the \strong{parallel} package (base R).
#' On Unix/macOS, \code{mc_cores > 1} uses \code{mclapply}. On Windows,
#' set \code{mc_cores = 1} (fork-based parallelism is unsupported on Windows).
#'
#' Memory tip: for very large datasets (> 100k cells per sample), pass
#' BPCells \code{IterableMatrix} objects in \code{sample_list}; LogicComm
#' will stream them chunk by chunk.
#'
#' @param sample_list Named list of expression inputs per sample. Each element
#'   can be:
#'   \itemize{
#'     \item A numeric or sparse matrix (genes x cells)
#'     \item A \code{Seurat} object
#'     \item A \code{BPCells} \code{IterableMatrix}
#'   }
#' @param group_info Named character vector: sample name -> group label
#'   (e.g. \code{c(s1="Case", s2="Ctrl")}). Must cover all names in
#'   \code{sample_list}. Alternatively, an unnamed vector in the same order.
#' @param lr_db LR database. Default: \code{lr_pairs_human}.
#' @param rank_threshold Per-cell REO anchor quantile. Default: \code{0.5}.
#' @param anchor_genes Optional gene universe for REO anchors. Default \code{NULL}
#'   uses all genes in each sample.
#' @param lcs_threshold Minimum LCS for "active" call. Default: \code{0.01}.
#' @param case_label Label for case group. Default: \code{"Case"}.
#' @param ctrl_label Label for control group. Default: \code{"Ctrl"}.
#' @param knn_list Optional named list of KNN matrices (cells x cells),
#'   one per sample. If \code{NULL}, Seurat graphs are used when available,
#'   otherwise global LCS mode is used.
#' @param graph_name Optional Seurat graph name passed to
#'   \code{IdentifyLogicConsensus()}, for example \code{"RNA_nn"}.
#' @param layer Seurat assay layer/slot for REO conversion. Default: \code{"counts"}.
#' @param remove_self_edges Remove KNN self-loops before LCS scoring. Default: \code{TRUE}.
#' @param graph_symmetrize How to symmetrize KNN/SNN graphs before scoring.
#' @param edge_weight_mode \code{"binary"} or \code{"weighted"}.
#' @param chunk_size Cells per chunk for REO computation. Default: \code{5000}.
#' @param mc_cores Number of parallel cores. Default: \code{1}.
#'   Set > 1 on Unix/macOS only.
#' @param min_samples_per_group Minimum number of non-missing samples per group
#'   required by \code{CompareLogicGroups}. Default: \code{1}.
#' @param verbose Logical. Default: \code{TRUE}.
#'
#' @return A list of class \code{LogicCommMulti} with:
#'   \describe{
#'     \item{lcs_list}{Named list of LCS vectors (one per sample)}
#'     \item{comparison}{A \code{\link{CompareLogicGroups}} result data frame}
#'     \item{group_info}{Group assignments used}
#'     \item{params}{List of parameters used}
#'   }
#'
#' @seealso \code{\link{CompareLogicGroups}}, \code{\link{plot_lcs_bubble}}
#'
#' @examples
#' make_sample <- function(l_value, r_value) {
#'   matrix(
#'     c(l_value, 0,
#'       0, r_value,
#'       1, 1),
#'     nrow = 3, byrow = TRUE,
#'     dimnames = list(c("L", "R", "BG"), c("C1", "C2"))
#'   )
#' }
#' samples <- list(
#'   Case1 = make_sample(10, 10),
#'   Ctrl1 = make_sample(5, 5)
#' )
#' groups <- c(Case1 = "Case", Ctrl1 = "Ctrl")
#' lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
#'                     stringsAsFactors = FALSE)
#' lr_db$ligand_genes <- list("L")
#' lr_db$receptor_genes <- list("R")
#' result <- run_multisample(samples, group_info = groups, lr_db = lr_db,
#'                           verbose = FALSE)
#' print(result$comparison)
#'
#' @importFrom parallel mclapply
#' @export
run_multisample <- function(sample_list,
                             group_info,
                             lr_db          = lr_pairs_human,
                             rank_threshold = 0.5,
                             anchor_genes   = NULL,
                             lcs_threshold  = 0.01,
                             case_label     = "Case",
                             ctrl_label     = "Ctrl",
                             knn_list       = NULL,
                             graph_name     = NULL,
                             layer          = "counts",
                             remove_self_edges = TRUE,
                             graph_symmetrize = c("none", "or", "max"),
                             edge_weight_mode = c("binary", "weighted"),
                             chunk_size     = 5000,
                             mc_cores       = 1,
                             min_samples_per_group = 1,
                             verbose        = TRUE) {

  graph_symmetrize <- match.arg(graph_symmetrize)
  edge_weight_mode <- match.arg(edge_weight_mode)

  # 1. Validate
  if (is.null(names(sample_list))) {
    names(sample_list) <- paste0("Sample", seq_along(sample_list))
  }
  if (is.null(names(group_info))) {
    if (length(group_info) != length(sample_list))
      stop("Unnamed group_info must have same length as sample_list.")
    names(group_info) <- names(sample_list)
  }
  group_info <- group_info[names(sample_list)]
  if (anyNA(group_info)) {
    missing <- names(sample_list)[is.na(group_info)]
    stop("group_info must provide a label for every sample in sample_list; missing labels for: ",
         paste(utils::head(missing, 5), collapse = ", "), ".")
  }

  lr_genes <- all_lr_genes(lr_db)
  n_samples <- length(sample_list)

  if (verbose) {
    n_case <- sum(group_info == case_label)
    n_ctrl <- sum(group_info == ctrl_label)
    message(sprintf(
      "[MultiSample] %d samples: %d %s + %d %s | %d L-R pairs | %d cores",
      n_samples, n_case, case_label, n_ctrl, ctrl_label,
      nrow(lr_db), mc_cores))
  }

  # 2. Per-sample REO + LCS
  sample_names <- names(sample_list)

  run_one <- function(sname) {
    if (verbose) message(sprintf("  [%s] Building REO matrix...", sname))
    obj <- sample_list[[sname]]

    reo <- tryCatch(
      calc_REO_matrix(obj, lr_genes = lr_genes,
                      rank_threshold = rank_threshold,
                      anchor_genes = anchor_genes,
                      layer = layer,
                      chunk_size = chunk_size, verbose = FALSE),
      error = function(e) {
        warning(sprintf("Sample '%s' REO failed: %s", sname, conditionMessage(e)))
        return(NULL)
      })
    if (is.null(reo)) return(NULL)

    knn <- if (!is.null(knn_list)) knn_list[[sname]] else NULL
    seurat_obj <- if (inherits(obj, "Seurat")) obj else NULL

    if (verbose) message(sprintf("  [%s] Computing LCS...", sname))
    lcs <- tryCatch(
      IdentifyLogicConsensus(reo, seurat_obj = seurat_obj, knn_mat = knn,
                             graph_name = graph_name,
                             lr_db = lr_db, lcs_threshold = lcs_threshold,
                             remove_self_edges = remove_self_edges,
                             graph_symmetrize = graph_symmetrize,
                             edge_weight_mode = edge_weight_mode,
                             verbose = FALSE),
      error = function(e) {
        warning(sprintf("Sample '%s' LCS failed: %s", sname, conditionMessage(e)))
        return(NULL)
      })
    lcs
  }

  # Parallel or sequential
  if (mc_cores > 1 && .Platform$OS.type == "unix") {
    lcs_list <- parallel::mclapply(sample_names, run_one,
                                    mc.cores = mc_cores,
                                    mc.preschedule = FALSE)
  } else {
    lcs_list <- lapply(sample_names, run_one)
  }
  names(lcs_list) <- sample_names

  # Remove failed samples
  failed <- vapply(lcs_list, is.null, logical(1))
  if (any(failed)) {
    warning(sprintf("%d sample(s) failed and were excluded: %s",
                    sum(failed), paste(sample_names[failed], collapse=", ")))
    lcs_list   <- lcs_list[!failed]
    group_info <- group_info[names(lcs_list)]
  }
  if (length(lcs_list) == 0) stop("All samples failed. Check input data.")

  # 3. Group comparison
  if (verbose) message("[MultiSample] Running group comparison...")
  comparison <- CompareLogicGroups(
    lcs_list, group_info = group_info,
    case_label = case_label, ctrl_label = ctrl_label,
    lcs_threshold = lcs_threshold, lr_db = lr_db,
    min_samples_per_group = min_samples_per_group,
    verbose = verbose)

  # 4. Assemble output
  result <- structure(
    list(
      lcs_list   = lcs_list,
      comparison = comparison,
      group_info = group_info,
      params     = list(
        rank_threshold = rank_threshold,
        anchor_n_genes = if (is.null(anchor_genes)) NA_integer_ else length(anchor_genes),
        lcs_threshold  = lcs_threshold,
        case_label     = case_label,
        ctrl_label     = ctrl_label,
        graph_name     = graph_name,
        layer          = layer,
        remove_self_edges = remove_self_edges,
        graph_symmetrize = graph_symmetrize,
        edge_weight_mode = edge_weight_mode,
        n_lr_pairs     = nrow(lr_db),
        n_samples      = length(lcs_list)
      )
    ),
    class = "LogicCommMulti"
  )

  if (verbose) {
    n_sig <- sum(comparison$fdr_fisher < 0.05 & comparison$asymmetry > 0,
                 na.rm = TRUE)
    message(sprintf(
      "[MultiSample] Done. %d significant %s-enriched pairs (FDR<0.05).",
      n_sig, case_label))
  }
  result
}

#' Print method for LogicCommMulti
#' @param x LogicCommMulti object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}, called for side effects.
#' @examples
#' comparison <- structure(
#'   data.frame(lr_pair = "L_R", case_freq = 0.8, ctrl_freq = 0.2,
#'              asymmetry = 0.6, log2fc_lcs = 1, fdr_fisher = 0.01),
#'   class = c("LogicCommResult", "data.frame"),
#'   case_label = "Case", ctrl_label = "Ctrl"
#' )
#' multi <- structure(
#'   list(comparison = comparison,
#'        params = list(n_samples = 2, n_lr_pairs = 1, case_label = "Case",
#'                      ctrl_label = "Ctrl", rank_threshold = 0.5,
#'                      lcs_threshold = 0.01)),
#'   class = "LogicCommMulti"
#' )
#' print(multi)
#' @export
print.LogicCommMulti <- function(x, ...) {
  p <- x$params
  cat(sprintf("LogicCommMulti | %d samples | %d L-R pairs\n",
              p$n_samples, p$n_lr_pairs))
  cat(sprintf("Groups: %s vs %s | rank_threshold=%.2f | lcs_threshold=%.3f\n",
              p$case_label, p$ctrl_label, p$rank_threshold, p$lcs_threshold))
  cat("--- Top enriched pairs --------------------------\n")
  print.LogicCommResult(x$comparison, n = 5)
  invisible(x)
}

#' Extract the comparison result from a LogicCommMulti object
#' @param x LogicCommMulti object
#' @return LogicCommResult data frame
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
#' fun <- get("get_comparison")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
get_comparison <- function(x) {
  stopifnot(inherits(x, "LogicCommMulti"))
  x$comparison
}

#' Export LogicComm results to a CSV file
#'
#' @param result A \code{LogicCommResult} or \code{LogicCommMulti} object.
#' @param file Output file path (e.g. \code{"results.csv"}).
#' @param sig_only If \code{TRUE}, only export significant pairs. Default: \code{FALSE}.
#' @param fdr_cutoff FDR cutoff for \code{sig_only}. Default: \code{0.05}.
#' @return Invisibly returns the data frame written.
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
#' fun <- get("export_results")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
export_results <- function(result, file, sig_only = FALSE, fdr_cutoff = 0.05) {
  if (inherits(result, "LogicCommMulti")) result <- result$comparison
  stopifnot(inherits(result, "LogicCommResult"))

  df <- result
  # Drop list columns if any
  df <- df[, !vapply(df, is.list, logical(1)), drop = FALSE]

  if (sig_only) {
    df <- df[!is.na(df$fdr_fisher) & df$fdr_fisher <= fdr_cutoff, ]
  }
  utils::write.csv(df, file = file, row.names = FALSE)
  if (nrow(df) > 0) {
    message(sprintf("Exported %d L-R pairs to: %s", nrow(df), file))
  }
  invisible(df)
}
