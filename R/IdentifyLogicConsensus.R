# R/IdentifyLogicConsensus.R

#' Identify Logic Consensus Scores (LCS) for all L-R pairs
#'
#' For each ligand-receptor pair in the database, computes a Logic Consensus
#' Score (LCS): the fraction of cells in which the complete ligand logic and the
#' complete receptor logic are simultaneously realized (co-expression).
#'
#' @details
#' \deqn{LCS = \frac{\text{cells with both L=TRUE AND R=TRUE}}{n_{cells}}}
#'
#' LogicComm scores communication from REO co-expression and does not aggregate
#' over a per-cell KNN/SNN neighborhood graph: for dissociated scRNA-seq that
#' graph lives in expression space, not physical space, so it cannot license
#' spatial juxtacrine/paracrine claims.
#'
#' **Multimer complexes:** A complex like "TGFBR1+TGFBR2" requires ALL subunit
#' genes to be active (\code{TRUE}) for the receptor to be counted as expressed.
#'
#' @param reo_mat Binary REO matrix (genes x cells), output of
#'   \code{\link{calc_REO_matrix}}.
#' @param lr_db Ligand-receptor database. Data frame with columns
#'   \code{lr_pair}, \code{ligand_genes}, \code{receptor_genes} (list columns).
#'   Defaults to the built-in \code{lr_pairs_human}.
#' @param lcs_threshold Minimum LCS for a pair to be considered "active" in
#'   downstream comparisons. Default: \code{0.01}.
#' @param verbose Print progress messages. Default: \code{TRUE}.
#' @param ... Deprecated neighborhood arguments (\code{knn_mat}, \code{seurat_obj},
#'   \code{graph_name}, \code{remove_self_edges}, \code{graph_symmetrize},
#'   \code{edge_weight_mode}), accepted for backward compatibility and ignored
#'   with a warning.
#'
#' @return A named numeric vector of LCS values, one per L-R pair in
#'   \code{lr_db}. Names are \code{lr_pair} identifiers. Pairs where neither
#'   the ligand nor receptor is found in \code{reo_mat} get \code{NA}.
#'
#' @seealso \code{\link{calc_REO_matrix}}, \code{\link{CompareLogicGroups}}
#'
#' @examples
#' \dontrun{
#' data(lr_pairs_human)
#' lr_genes <- all_lr_genes(lr_pairs_human)
#' reo <- calc_REO_matrix(my_seurat, lr_genes = lr_genes)
#'
#' lcs <- IdentifyLogicConsensus(reo)
#'
#' # View top pairs
#' sort(lcs, decreasing = TRUE)[1:20]
#' }
#'
#' @importFrom Matrix summary
#' @export
IdentifyLogicConsensus <- function(reo_mat,
                                    lr_db         = lr_pairs_human,
                                    lcs_threshold = 0.01,
                                    verbose       = TRUE,
                                    ...) {

  stopifnot(inherits(reo_mat, c("dgCMatrix","REOMatrix","sparseMatrix","matrix")))
  .deprecate_neighborhood_args(list(...), "IdentifyLogicConsensus")

  if (verbose) message("[LCS] Global co-expression mode.")

  # ── 2. Pre-calculate complex logic ────────────────────────────────────────
  if (verbose) message("[LCS] Resolving complex logic for unique gene sets...")
  
  unique_complexes <- unique(c(lr_db$ligand_genes, lr_db$receptor_genes))
  complex_keys <- vapply(unique_complexes, function(gs) paste(sort(gs), collapse = "|"), character(1))
  
  # Map each unique complex to its logic vector
  complex_logic_map <- lapply(unique_complexes, function(gs) {
    .resolve_complex_logic(gs, reo_mat)
  })
  names(complex_logic_map) <- complex_keys
  
  # Helper to get pre-calculated logic
  get_logic <- function(gs) {
    key <- paste(sort(gs), collapse = "|")
    complex_logic_map[[key]]
  }

  # ── 3. Compute LCS per L-R pair ───────────────────────────────────────────
  n_pairs <- nrow(lr_db)
  lcs_vec <- setNames(rep(NA_real_, n_pairs), lr_db$lr_pair)

  if (verbose) message(sprintf("[LCS] Scoring %d L-R pairs...", n_pairs))

  for (i in seq_len(n_pairs)) {
    lig_genes <- lr_db$ligand_genes[[i]]
    rec_genes <- lr_db$receptor_genes[[i]]

    # Skip if no subunit found in REO matrix
    if (!any(lig_genes %in% rownames(reo_mat)) || !any(rec_genes %in% rownames(reo_mat))) {
      next
    }

    lig_logic <- get_logic(lig_genes)
    rec_logic <- get_logic(rec_genes)

    lcs_vec[i] <- compute_lcs_global(lig_logic, rec_logic)
  }

  n_active <- sum(!is.na(lcs_vec) & lcs_vec >= lcs_threshold, na.rm = TRUE)
  if (verbose) message(sprintf(
    "[LCS] Done. %d / %d pairs active (LCS >= %.3f).",
    n_active, n_pairs, lcs_threshold))

  class(lcs_vec) <- c("LCSVector", class(lcs_vec))
  lcs_vec
}

#' Print method for LCSVector
#' @param x LCSVector
#' @param ... ignored
#' @export
print.LCSVector <- function(x, ...) {
  n_total  <- length(x)
  n_scored <- sum(!is.na(x))
  n_active <- sum(x > 0.01, na.rm = TRUE)
  cat(sprintf(
    "LCSVector: %d pairs | scored: %d | active (>0.01): %d\n",
    n_total, n_scored, n_active))
  cat("Top 5:\n")
  top <- head(sort(x[!is.na(x)], decreasing = TRUE), 5)
  print(round(top, 4))
  invisible(x)
}
