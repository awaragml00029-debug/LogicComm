# R/CompareLogicGroups.R

#' Compare Logic Consensus Scores Between Two Sample Groups
#'
#' Takes a named list of per-sample LCS vectors and compares the frequency and
#' magnitude of ligand-receptor activity between Case and Control groups.
#'
#' @details
#' A sample is called positive for a ligand-receptor pair when its LCS is greater
#' than or equal to \code{lcs_threshold}. Missing LCS values are treated as
#' unavailable rather than negative: frequencies and Fisher tables use the number
#' of non-missing samples for that pair in each group.
#'
#' @param lcs_list Named list where each element is a numeric LCS vector returned
#'   by \code{IdentifyLogicConsensus()}.
#' @param group_info Named character vector mapping sample names to group labels.
#'   If unnamed, it must have the same order and length as \code{lcs_list}.
#' @param case_label Character. Label for the case group. Default: \code{"Case"}.
#' @param ctrl_label Character. Label for the control group. Default: \code{"Ctrl"}.
#' @param lcs_threshold LCS cutoff for calling a sample positive. Default:
#'   \code{0.01}.
#' @param min_case_freq Minimum Case frequency for inclusion in output. Default:
#'   \code{0}.
#' @param lr_db Optional LR database used to append metadata columns such as
#'   \code{ligand}, \code{receptor}, \code{pathway}, and \code{annotation}. If
#'   \code{NULL}, the built-in \code{lr_pairs_human} metadata is used when
#'   available.
#' @param min_samples_per_group Minimum number of non-missing samples required in
#'   each group for Fisher/Wilcoxon tests and frequencies. Pairs below this
#'   threshold receive \code{NA} frequencies and p-values of 1. Default:
#'   \code{1}.
#' @param verbose Print progress messages. Default: \code{TRUE}.
#'
#' @return A data frame (class \code{LogicCommResult}) with one row per L-R pair.
#'   Important columns include \code{case_freq}, \code{ctrl_freq},
#'   \code{asymmetry}, mean LCS values, Fisher/Wilcoxon p-values and FDRs, total
#'   group sizes (\code{n_case}, \code{n_ctrl}), and per-pair non-missing sample
#'   counts (\code{n_case_avail}, \code{n_ctrl_avail}).
#'
#' @seealso \code{\link{IdentifyLogicConsensus}}, \code{\link{plot_lcs_bubble}},
#'   \code{\link{plot_lcs_heatmap}}
#'
#' @examples
#' lcs_list <- list(
#'   Case1 = c(L_R = 0.4),
#'   Case2 = c(L_R = 0.3),
#'   Ctrl1 = c(L_R = 0.0),
#'   Ctrl2 = c(L_R = 0.1)
#' )
#' groups <- c(Case1 = "Case", Case2 = "Case", Ctrl1 = "Ctrl", Ctrl2 = "Ctrl")
#' result <- CompareLogicGroups(lcs_list, group_info = groups, verbose = FALSE)
#' subset(result, asymmetry >= 0.1)
#'
#' @export
CompareLogicGroups <- function(lcs_list,
                               group_info,
                               case_label    = "Case",
                               ctrl_label    = "Ctrl",
                               lcs_threshold = 0.01,
                               min_case_freq = 0,
                               lr_db         = NULL,
                               min_samples_per_group = 1,
                               verbose       = TRUE) {

  stopifnot(is.list(lcs_list), length(lcs_list) >= 2)
  stopifnot(is.numeric(lcs_threshold), length(lcs_threshold) == 1)
  stopifnot(is.numeric(min_samples_per_group), length(min_samples_per_group) == 1,
            min_samples_per_group >= 1)
  min_samples_per_group <- as.integer(min_samples_per_group)

  if (is.null(names(lcs_list))) names(lcs_list) <- paste0("S", seq_along(lcs_list))

  if (is.null(names(group_info))) {
    if (length(group_info) != length(lcs_list)) {
      stop("group_info length must match lcs_list length when unnamed.")
    }
    names(group_info) <- names(lcs_list)
  }
  group_info <- group_info[names(lcs_list)]
  if (anyNA(group_info)) {
    stop("group_info must contain labels for all samples in lcs_list.")
  }

  case_samples <- names(group_info)[group_info == case_label]
  ctrl_samples <- names(group_info)[group_info == ctrl_label]
  if (length(case_samples) == 0) stop("No samples found with case_label='", case_label, "'")
  if (length(ctrl_samples) == 0) stop("No samples found with ctrl_label='", ctrl_label, "'")

  n_case <- length(case_samples)
  n_ctrl <- length(ctrl_samples)
  if (verbose) message(sprintf("[Compare] %d %s samples vs %d %s samples.",
                               n_case, case_label, n_ctrl, ctrl_label))

  pair_lists <- lapply(lcs_list, function(v) {
    nm <- names(v)
    if (is.null(nm)) character(0) else nm
  })
  all_pairs <- Reduce(union, pair_lists)
  if (length(all_pairs) == 0) stop("No named LCS pairs found in lcs_list.")

  lcs_mat <- do.call(cbind, lapply(lcs_list, function(v) {
    out <- setNames(rep(NA_real_, length(all_pairs)), all_pairs)
    nm <- names(v)
    if (!is.null(nm) && length(nm) > 0) out[nm] <- as.numeric(v)
    out
  }))
  rownames(lcs_mat) <- all_pairs
  colnames(lcs_mat) <- names(lcs_list)

  if (verbose) message(sprintf("[Compare] Computing stats for %d L-R pairs...",
                               length(all_pairs)))

  n_pairs <- length(all_pairs)
  case_mat <- lcs_mat[, case_samples, drop = FALSE]
  ctrl_mat <- lcs_mat[, ctrl_samples, drop = FALSE]

  case_n_avail <- rowSums(!is.na(case_mat))
  ctrl_n_avail <- rowSums(!is.na(ctrl_mat))
  enough <- case_n_avail >= min_samples_per_group &
            ctrl_n_avail >= min_samples_per_group

  case_pos <- rowSums(case_mat >= lcs_threshold, na.rm = TRUE)
  ctrl_pos <- rowSums(ctrl_mat >= lcs_threshold, na.rm = TRUE)
  case_freq <- rep(NA_real_, n_pairs)
  ctrl_freq <- rep(NA_real_, n_pairs)
  case_freq[enough] <- case_pos[enough] / case_n_avail[enough]
  ctrl_freq[enough] <- ctrl_pos[enough] / ctrl_n_avail[enough]
  asymmetry <- case_freq - ctrl_freq

  case_mean <- rowMeans(case_mat, na.rm = TRUE)
  ctrl_mean <- rowMeans(ctrl_mat, na.rm = TRUE)
  case_mean[is.nan(case_mean)] <- NA_real_
  ctrl_mean[is.nan(ctrl_mean)] <- NA_real_
  log2fc <- log2((case_mean + 1e-6) / (ctrl_mean + 1e-6))

  p_fisher <- vapply(seq_len(n_pairs), function(i) {
    if (!enough[i]) return(1.0)
    tab <- matrix(c(case_pos[i], case_n_avail[i] - case_pos[i],
                    ctrl_pos[i], ctrl_n_avail[i] - ctrl_pos[i]),
                  nrow = 2, byrow = TRUE)
    tryCatch(stats::fisher.test(tab)$p.value, error = function(e) 1.0)
  }, numeric(1))

  p_wilcox <- vapply(seq_len(n_pairs), function(i) {
    if (!enough[i]) return(1.0)
    cv <- case_mat[i, ]; cv <- cv[!is.na(cv)]
    rv <- ctrl_mat[i, ]; rv <- rv[!is.na(rv)]
    if (length(cv) < 2 || length(rv) < 2) return(1.0)
    tryCatch(stats::wilcox.test(cv, rv, exact = FALSE)$p.value,
             error = function(e) 1.0)
  }, numeric(1))

  fdr_fisher <- .bh_adjust(p_fisher)
  fdr_wilcox <- .bh_adjust(p_wilcox)
  sig_stars <- .pval_stars(fdr_fisher)

  result <- data.frame(
    lr_pair       = all_pairs,
    case_freq     = round(case_freq, 4),
    ctrl_freq     = round(ctrl_freq, 4),
    asymmetry     = round(asymmetry, 4),
    case_mean_lcs = round(case_mean, 5),
    ctrl_mean_lcs = round(ctrl_mean, 5),
    log2fc_lcs    = round(log2fc, 4),
    p_fisher      = signif(p_fisher, 4),
    fdr_fisher    = signif(fdr_fisher, 4),
    p_wilcox      = signif(p_wilcox, 4),
    fdr_wilcox    = signif(fdr_wilcox, 4),
    sig_stars     = sig_stars,
    n_case        = n_case,
    n_ctrl        = n_ctrl,
    n_case_avail  = case_n_avail,
    n_ctrl_avail  = ctrl_n_avail,
    stringsAsFactors = FALSE
  )

  meta_source <- lr_db
  if (is.null(meta_source) && exists("lr_pairs_human", inherits = TRUE)) {
    meta_source <- lr_pairs_human
  }
  result <- .append_lr_metadata(result, meta_source)

  keep <- !is.na(result$case_freq) & result$case_freq >= min_case_freq
  result <- result[keep, , drop = FALSE]
  result <- result[order(-result$asymmetry, result$fdr_fisher, na.last = TRUE), ]
  rownames(result) <- NULL

  class(result) <- c("LogicCommResult", "data.frame")
  attr(result, "case_label") <- case_label
  attr(result, "ctrl_label") <- ctrl_label
  attr(result, "lcs_threshold") <- lcs_threshold
  attr(result, "min_samples_per_group") <- min_samples_per_group
  attr(result, "lcs_mat") <- lcs_mat

  n_sig <- sum(result$fdr_fisher < 0.05 & result$asymmetry > 0, na.rm = TRUE)
  if (verbose) message(sprintf(
    "[Compare] Done. %d pairs enriched in %s (FDR<0.05, asymmetry>0).",
    n_sig, case_label))

  result
}

#' Append LR metadata while preserving result row order
#' @keywords internal
.append_lr_metadata <- function(result, lr_db = NULL) {
  if (is.null(lr_db) || !is.data.frame(lr_db) || !"lr_pair" %in% names(lr_db)) {
    return(result)
  }
  meta_cols <- intersect(c("lr_pair", "ligand", "receptor", "pathway", "annotation"),
                         names(lr_db))
  if (length(meta_cols) <= 1) return(result)
  meta <- lr_db[, meta_cols, drop = FALSE]
  meta <- meta[!duplicated(meta$lr_pair), , drop = FALSE]
  m <- match(result$lr_pair, meta$lr_pair)
  for (col in setdiff(meta_cols, "lr_pair")) {
    result[[col]] <- meta[[col]][m]
  }
  result
}

#' Print method for LogicCommResult
#' @param x LogicCommResult data frame.
#' @param n Number of rows to print. Default: 10.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}, called for side effects.
#' @examples
#' result <- structure(
#'   data.frame(lr_pair = "L_R", case_freq = 0.8, ctrl_freq = 0.2,
#'              asymmetry = 0.6, log2fc_lcs = 1, fdr_fisher = 0.01),
#'   class = c("LogicCommResult", "data.frame"),
#'   case_label = "Case", ctrl_label = "Ctrl"
#' )
#' print(result)
#' @export
print.LogicCommResult <- function(x, n = 10, ...) {
  case_l <- attr(x, "case_label")
  ctrl_l <- attr(x, "ctrl_label")
  cat(sprintf("LogicCommResult: %d L-R pairs | %s vs %s\n",
              nrow(x), case_l, ctrl_l))
  sig <- sum(x$fdr_fisher < 0.05 & x$asymmetry > 0, na.rm = TRUE)
  cat(sprintf("Significant (FDR<0.05, asymmetry>0): %d\n\n", sig))

  show_cols <- intersect(
    c("lr_pair", "pathway", "case_freq", "ctrl_freq", "asymmetry", "log2fc_lcs",
      "fdr_fisher", "sig_stars", "n_case_avail", "n_ctrl_avail"),
    names(x))
  print(utils::head(as.data.frame(x)[, show_cols, drop = FALSE], n), row.names = FALSE)
  invisible(x)
}

#' Filter LogicCommResult by significance or asymmetry
#'
#' @param result A \code{LogicCommResult} object.
#' @param min_asymmetry Minimum absolute asymmetry. Default: \code{0.3}.
#' @param max_fdr Maximum BH-adjusted FDR (Fisher). Default: \code{0.05}.
#' @param direction \code{"up"} (Case > Ctrl), \code{"down"} (Ctrl > Case),
#'   or \code{"both"}. Default: \code{"up"}.
#' @return Filtered \code{LogicCommResult}.
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
#' fun <- get("filter_lcs")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
filter_lcs <- function(result,
                       min_asymmetry = 0.3,
                       max_fdr       = 0.05,
                       direction     = "up") {
  stopifnot(inherits(result, "LogicCommResult"))
  direction <- match.arg(direction, c("up", "down", "both"))
  idx <- result$fdr_fisher <= max_fdr & abs(result$asymmetry) >= min_asymmetry
  if (direction == "up") idx <- idx & result$asymmetry > 0
  if (direction == "down") idx <- idx & result$asymmetry < 0
  idx[is.na(idx)] <- FALSE

  out <- result[idx, , drop = FALSE]
  keep_attrs <- attributes(out)[c("names", "row.names")]
  attributes(out) <- keep_attrs
  class(out) <- c("LogicCommResult", "data.frame")
  for (nm in c("case_label", "ctrl_label", "lcs_threshold",
               "min_samples_per_group", "lcs_mat")) {
    attr(out, nm) <- attr(result, nm)
  }
  rownames(out) <- NULL
  out
}
