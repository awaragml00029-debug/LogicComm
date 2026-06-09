# R/discover.R

#' One-call cell-type communication discovery pipeline
#'
#' Wires the LogicComm discovery workflow end to end with footgun-free defaults:
#' cell-type co-expression scoring, specificity + proliferation-confound
#' annotation, an \strong{axis-level} permutation null (so every sender ->
#' receiver -> L-R axis gets its own empirical p and BH FDR, rather than a single
#' cell-type-pair p-value broadcast across all of its L-R pairs), evidence
#' ranking, and a confound-filtered discovery view plus an FDR-passing shortlist.
#'
#' @details
#' This is the recommended entry point. It exists mainly to make the
#' axis-level permutation the default and to fold in the proliferation /
#' transcriptional-breadth confound, which are easy to miss when calling the
#' steps by hand (e.g. passing \code{metric = "sum_lcs"} silently gives a
#' pair-level null). Equivalent to:
#' \preformatted{
#'   ct   <- summarize_celltype_communication(reo, cell_labels, lr_db)
#'   ct   <- score_communication_specificity(ct)
#'   ct   <- diagnose_proliferation_confound(ct, expr = counts)   # if expr given
#'   null <- permute_celltype_communication(ct, reo_mat = reo, metric = "lcs",
#'                                          n_perm = 1000, n_cores = 8)
#'   ranked <- rank_communication_axes(ct, null_pair = null)
#'   view   <- communication_discovery_view(ranked, drop_proliferation = TRUE)
#' }
#'
#' @param reo_mat Binary REO matrix (genes x cells) or \code{LogicCommREOResult},
#'   from \code{\link{calc_REO_matrix}}.
#' @param cell_labels Named cell-type labels (or \code{NULL} to take them from a
#'   \code{LogicCommREOResult}/Seurat-derived workflow upstream).
#' @param lr_db Ligand-receptor database.
#' @param expr Optional gene x cell expression matrix (raw counts or
#'   log-normalized). When supplied, the proliferation / transcriptional-breadth
#'   confound is diagnosed and (by default) dropped from the discovery view.
#'   Strongly recommended whenever cycling cells are present (e.g. a
#'   \code{Cycling_Lymphocytes} cluster), which otherwise dominate the ranking.
#' @param n_perm,n_cores,seed Permutation settings. Co-expression nulls need no
#'   graph rescan, so a publication-grade \code{n_perm} (e.g. 1000) is affordable;
#'   \code{n_cores > 1} forks the loop on Unix/macOS.
#' @param fdr_cutoff Axes with \code{permutation_fdr <= fdr_cutoff} are returned in
#'   \code{shortlist}. Default \code{0.1}.
#' @param drop_broad,drop_proliferation,drop_identity Confound filters passed to
#'   \code{\link{communication_discovery_view}}. \code{drop_proliferation} is
#'   honoured only when \code{expr} was supplied.
#' @param top_n Optional cap on the discovery view.
#' @param verbose Print progress.
#' @param ... Passed to \code{\link{summarize_celltype_communication}}.
#' @return A list with \code{ct_comm}, \code{null} (axis-level permutation),
#'   \code{ranked}, \code{view} (confound-filtered discovery view), and
#'   \code{shortlist} (axes passing \code{permutation_fdr <= fdr_cutoff}, ordered
#'   by FDR then discovery score).
#' @seealso \code{\link{summarize_celltype_communication}},
#'   \code{\link{permute_celltype_communication}},
#'   \code{\link{rank_communication_axes}},
#'   \code{\link{communication_discovery_view}},
#'   \code{\link{diagnose_proliferation_confound}}
#' @export
discover_celltype_communication <- function(reo_mat,
                                            cell_labels = NULL,
                                            lr_db = lr_pairs_human,
                                            expr = NULL,
                                            n_perm = 1000,
                                            n_cores = 1L,
                                            seed = NULL,
                                            fdr_cutoff = 0.1,
                                            drop_broad = TRUE,
                                            drop_proliferation = TRUE,
                                            drop_identity = FALSE,
                                            top_n = NULL,
                                            verbose = TRUE,
                                            ...) {
  ct <- summarize_celltype_communication(reo_mat, cell_labels = cell_labels,
                                         lr_db = lr_db, verbose = verbose, ...)
  ct <- score_communication_specificity(ct, verbose = verbose)

  did_prolif <- FALSE
  if (!is.null(expr)) {
    ct <- diagnose_proliferation_confound(ct, expr = expr, verbose = verbose)
    did_prolif <- TRUE
  } else if (isTRUE(drop_proliferation) && isTRUE(verbose)) {
    message("discover_celltype_communication: no 'expr' supplied, so the ",
            "proliferation / transcriptional-breadth confound cannot be diagnosed. ",
            "Pass the expression matrix as 'expr' when cycling cells are present ",
            "(e.g. a Cycling_Lymphocytes cluster), or they will dominate the ranking.")
  }

  # Axis-level null (metric = 'lcs'): each L-R axis gets its own empirical p / FDR.
  null <- permute_celltype_communication(ct, reo_mat = reo_mat, lr_db = lr_db,
                                         n_perm = n_perm, n_cores = n_cores,
                                         seed = seed, metric = "lcs", verbose = verbose)
  ranked <- rank_communication_axes(ct, null_pair = null)
  view <- communication_discovery_view(ranked,
                                       drop_broad = drop_broad,
                                       drop_proliferation = isTRUE(drop_proliferation) && did_prolif,
                                       drop_identity = drop_identity,
                                       top_n = top_n, verbose = verbose)

  if ("permutation_fdr" %in% names(view)) {
    keep <- is.finite(view$permutation_fdr) & view$permutation_fdr <= fdr_cutoff
    shortlist <- view[keep, , drop = FALSE]
    shortlist <- shortlist[order(shortlist$permutation_fdr, -shortlist$discovery_score), , drop = FALSE]
  } else {
    shortlist <- view[0, , drop = FALSE]
  }
  if (isTRUE(verbose)) {
    message(sprintf("[discover] %d axis(es) pass permutation_fdr <= %.3g (of %d in the discovery view).",
                    nrow(shortlist), fdr_cutoff, nrow(view)))
  }

  list(ct_comm = ct, null = null, ranked = ranked, view = view, shortlist = shortlist)
}
