# R/modulators.R — CellChatDB modulator helpers for LogicComm

#' Resolve rank percentiles for a gene or complex
#'
#' @param genes Character vector of genes/subunits.
#' @param rank_mat Numeric rank-percentile matrix, genes x cells.
#' @param require_all_subunits Logical; if TRUE, missing subunits return NA.
#' @param aggregate Aggregation for multi-subunit ranks. "min" is conservative.
#' @return Numeric vector of per-cell rank percentiles.
#' @keywords internal
.resolve_complex_rank <- function(genes, rank_mat, require_all_subunits = TRUE,
                                  aggregate = c("min", "mean", "max")) {
  aggregate <- match.arg(aggregate)
  genes <- unique(genes[nzchar(genes) & !is.na(genes)])
  if (length(genes) == 0 || is.null(rank_mat)) return(rep(NA_real_, ncol(rank_mat)))
  available <- intersect(genes, rownames(rank_mat))
  if (length(available) == 0) return(rep(NA_real_, ncol(rank_mat)))
  if (require_all_subunits && length(available) < length(genes)) {
    return(rep(NA_real_, ncol(rank_mat)))
  }
  mat_sub <- rank_mat[available, , drop = FALSE]
  if (nrow(mat_sub) == 1) return(as.numeric(mat_sub[1, ]))
  fun <- switch(aggregate, min = min, mean = mean, max = max)
  res <- apply(mat_sub, 2, fun, na.rm = TRUE)
  # A fully-NA column yields Inf/-Inf from min/max with na.rm; treat as missing
  # so it cannot contaminate downstream rank comparisons.
  res[!is.finite(res)] <- NA_real_
  res
}

#' Get a list-column entry safely
#' @keywords internal
.lr_list_entry <- function(lr_db, col, i) {
  if (!col %in% names(lr_db)) return(character(0))
  x <- lr_db[[col]][[i]]
  x <- as.character(x)
  x[nzchar(x) & !is.na(x)]
}

#' Check whether an LR row has any genes in a modulator list column
#' @keywords internal
.has_lr_genes <- function(lr_db, col, i) {
  length(.lr_list_entry(lr_db, col, i)) > 0
}

#' Summarize CellChatDB modulator activity in a LogicComm dataset
#'
#' Reports whether converted CellChatDB modulators are present, active, and (when
#' rank percentiles are supplied) likely to block ligand-receptor logic by rank
#' order. This function does not change LCS values; it is intended for diagnostic
#' interpretation before running \code{IdentifyLogicGateConsensus()}.
#'
#' @param reo_mat Binary REO matrix, genes x cells.
#' @param lr_db LogicComm LR database, optionally with CellChatDB modulator list
#'   columns produced by \code{as_logiccomm_lr_db_from_cellchat()}.
#' @param rank_mat Optional rank-percentile matrix aligned to \code{reo_mat}.
#' @return Data frame with one row per LR pair.
#' @export
summarize_lr_modulators <- function(reo_mat, lr_db, rank_mat = NULL) {
  stopifnot(inherits(reo_mat, c("dgCMatrix", "sparseMatrix", "matrix")))
  if (!is.null(rank_mat)) {
    if (!identical(rownames(rank_mat), rownames(reo_mat)) ||
        !identical(colnames(rank_mat), colnames(reo_mat))) {
      stop("rank_mat must have the same rownames and colnames as reo_mat.")
    }
  }

  n_pairs <- nrow(lr_db)
  out <- vector("list", n_pairs)
  mod_cols <- c("agonist_genes", "antagonist_genes", "co_A_receptor_genes", "co_I_receptor_genes")

  for (i in seq_len(n_pairs)) {
    lig_genes <- .lr_list_entry(lr_db, "ligand_genes", i)
    rec_genes <- .lr_list_entry(lr_db, "receptor_genes", i)
    lig_rank <- if (!is.null(rank_mat)) .resolve_complex_rank(lig_genes, rank_mat) else NULL
    rec_rank <- if (!is.null(rank_mat)) .resolve_complex_rank(rec_genes, rank_mat) else NULL

    row <- data.frame(
      lr_pair = lr_db$lr_pair[i],
      has_agonist = .has_lr_genes(lr_db, "agonist_genes", i),
      has_antagonist = .has_lr_genes(lr_db, "antagonist_genes", i),
      has_co_A_receptor = .has_lr_genes(lr_db, "co_A_receptor_genes", i),
      has_co_I_receptor = .has_lr_genes(lr_db, "co_I_receptor_genes", i),
      agonist_active_rate = NA_real_,
      antagonist_active_rate = NA_real_,
      co_A_active_rate = NA_real_,
      co_I_active_rate = NA_real_,
      antagonist_outranks_ligand_rate = NA_real_,
      antagonist_outranks_receptor_rate = NA_real_,
      coI_outranks_receptor_rate = NA_real_,
      missing_modulator_genes = NA_character_,
      stringsAsFactors = FALSE
    )

    missing <- character(0)
    for (col in mod_cols) {
      genes <- .lr_list_entry(lr_db, col, i)
      missing <- c(missing, setdiff(genes, rownames(reo_mat)))
    }
    row$missing_modulator_genes <- paste(unique(missing), collapse = ";")

    active_specs <- list(
      agonist_active_rate = "agonist_genes",
      antagonist_active_rate = "antagonist_genes",
      co_A_active_rate = "co_A_receptor_genes",
      co_I_active_rate = "co_I_receptor_genes"
    )
    for (nm in names(active_specs)) {
      genes <- .lr_list_entry(lr_db, active_specs[[nm]], i)
      if (length(genes) > 0) {
        vec <- .resolve_complex_logic(genes, reo_mat, require_all_subunits = FALSE, mode = "any")
        row[[nm]] <- mean(vec, na.rm = TRUE)
      }
    }

    if (!is.null(rank_mat)) {
      ant_genes <- .lr_list_entry(lr_db, "antagonist_genes", i)
      if (length(ant_genes) > 0) {
        ant_rank <- .resolve_complex_rank(ant_genes, rank_mat, require_all_subunits = FALSE, aggregate = "max")
        row$antagonist_outranks_ligand_rate <- mean(!is.na(ant_rank) & !is.na(lig_rank) & ant_rank >= lig_rank)
        row$antagonist_outranks_receptor_rate <- mean(!is.na(ant_rank) & !is.na(rec_rank) & ant_rank >= rec_rank)
      }
      coi_genes <- .lr_list_entry(lr_db, "co_I_receptor_genes", i)
      if (length(coi_genes) > 0) {
        coi_rank <- .resolve_complex_rank(coi_genes, rank_mat, require_all_subunits = FALSE, aggregate = "max")
        row$coI_outranks_receptor_rate <- mean(!is.na(coi_rank) & !is.na(rec_rank) & coi_rank >= rec_rank)
      }
    }
    out[[i]] <- row
  }
  do.call(rbind, out)
}
