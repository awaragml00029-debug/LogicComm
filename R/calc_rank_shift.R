# R/calc_rank_shift.R

#' Compute Gene Expression Rank Shift Between Groups
#'
#' Detects genes whose **within-cell expression rank** changes dramatically
#' between Case and Control groups — even if the absolute fold-change is modest.
#'
#' @details
#' **Motivation:** A gene whose raw expression only doubles between conditions
#' may still represent a major signaling shift if its rank within the cell's
#' transcriptome jumps from position 20,000 to position 500. Standard
#' differential expression tests miss this because they operate on absolute
#' values.
#'
#' **Algorithm:**
#' \enumerate{
#'   \item For each cell, rank all expressed genes (1 = highest expressed).
#'   \item Normalize rank to [0,1]: \eqn{r_{ij} = \text{rank}_{ij} / n_{\text{expressed},j}}
#'   \item For each gene and each sample, compute the median normalized rank
#'         across cells that express the gene (excludes dropout zeros).
#'   \item Compare Case vs Control sample medians via Wilcoxon test.
#'   \item Rank shift = median(Ctrl ranks) - median(Case ranks).
#'         A positive shift means the gene moved to a **higher** (lower
#'         number = more dominant) rank in Case — i.e., it became more prominent.
#' }
#'
#' @param sample_list Named list of expression matrices or Seurat objects.
#' @param group_info Named character vector mapping sample names to group labels.
#' @param genes Character vector of genes to analyze. Defaults to all L-R genes
#'   in \code{lr_db}.
#' @param lr_db LR database for default gene set. Default: \code{lr_pairs_human}.
#' @param case_label Case group label. Default: \code{"Case"}.
#' @param ctrl_label Control group label. Default: \code{"Ctrl"}.
#' @param min_detection_frac Minimum fraction of cells expressing a gene
#'   (non-zero) for it to be included. Default: \code{0.05}.
#' @param layer Seurat layer. Default: \code{"counts"}.
#' @param verbose Logical. Default: \code{TRUE}.
#'
#' @return A data frame (class \code{RankShiftResult}) sorted by \code{shift_score},
#'   with columns:
#'   \describe{
#'     \item{gene}{Gene symbol}
#'     \item{role}{Whether gene is "ligand", "receptor", "both", or "other" in lr_db}
#'     \item{median_rank_case}{Median normalized rank across Case samples (0=top, 1=bottom)}
#'     \item{median_rank_ctrl}{Median normalized rank across Control samples}
#'     \item{shift_score}{ctrl_rank - case_rank. Positive = gene moved UP in Case}
#'     \item{p_wilcox, fdr_wilcox}{Wilcoxon test p-value and BH-adjusted FDR}
#'     \item{mean_expr_case, mean_expr_ctrl}{Mean expression per group}
#'     \item{log2fc_expr}{Traditional log2 fold-change in expression}
#'     \item{sig_stars}{Significance annotation}
#'   }
#'
#' @seealso \code{\link{plot_rank_shift}}, \code{\link{CompareLogicGroups}}
#'
#' @examples
#' \dontrun{
#' groups <- c(rep("Case", 10), rep("Ctrl", 10))
#' names(groups) <- names(sample_list)
#'
#' rs <- calc_rank_shift(sample_list, group_info = groups)
#' head(rs[rs$shift_score > 0.2, ])  # genes that rose in rank in Case
#' plot_rank_shift(rs)
#' }
#'
#' @export
calc_rank_shift <- function(sample_list,
                             group_info,
                             genes              = NULL,
                             lr_db              = lr_pairs_human,
                             case_label         = "Case",
                             ctrl_label         = "Ctrl",
                             min_detection_frac = 0.05,
                             layer              = "counts",
                             verbose            = TRUE) {

  # ── 1. Setup ──────────────────────────────────────────────────────────────
  if (is.null(names(sample_list))) {
    names(sample_list) <- paste0("S", seq_along(sample_list))
  }
  if (is.null(names(group_info))) {
    names(group_info) <- names(sample_list)
  }
  group_info <- group_info[names(sample_list)]

  case_s <- names(group_info)[group_info == case_label]
  ctrl_s <- names(group_info)[group_info == ctrl_label]
  if (length(case_s) == 0 || length(ctrl_s) == 0)
    stop("No samples found for one or both group labels.")

  if (is.null(genes)) {
    genes <- all_lr_genes(lr_db)
    if (verbose) message(sprintf("[RankShift] Using %d L-R genes.", length(genes)))
  }

  # ── 2. Per-sample median normalized rank ─────────────────────────────────
  if (verbose) message(sprintf("[RankShift] Processing %d samples...", length(sample_list)))

  sample_ranks <- lapply(names(sample_list), function(sname) {
    mat <- .extract_matrix(sample_list[[sname]], layer = layer)
    if (length(intersect(genes, rownames(mat))) == 0) return(NULL)
    # Rank each gene within the cell's FULL transcriptome, then keep the genes of
    # interest. Subsetting to the L-R panel before ranking would rank each gene
    # only among the panel (so a transcriptome-dominant ligand among ~120 panel
    # genes would look mid-ranked), which contradicts the documented method and
    # the within-cell REO philosophy used throughout LogicComm.
    .median_norm_rank_per_gene(mat, min_detection_frac, target_genes = genes)
  })
  names(sample_ranks) <- names(sample_list)

  # Remove NULLs
  sample_ranks <- sample_ranks[!vapply(sample_ranks, is.null, logical(1))]
  if (length(sample_ranks) == 0) {
    stop("No analyzable genes found in any sample. Check genes, rownames, and min_detection_frac.")
  }
  if (length(intersect(case_s, names(sample_ranks))) == 0 ||
      length(intersect(ctrl_s, names(sample_ranks))) == 0) {
    stop("No analyzable samples remain in one or both groups after gene filtering.")
  }

  # ── 3. Build rank matrix (genes x samples) ───────────────────────────────
  all_genes_found <- Reduce(union, lapply(sample_ranks, names))
  if (length(all_genes_found) == 0) {
    stop("No genes passed min_detection_frac in any analyzable sample.")
  }
  rank_mat <- matrix(NA_real_, nrow = length(all_genes_found),
                     ncol = length(sample_ranks),
                     dimnames = list(all_genes_found, names(sample_ranks)))
  for (sname in names(sample_ranks)) {
    rv <- sample_ranks[[sname]]
    rank_mat[names(rv), sname] <- rv
  }

  expr_means <- .compute_expr_means(sample_list, all_genes_found, group_info,
                                    case_label, ctrl_label, layer)

  # ── 4. Statistics ─────────────────────────────────────────────────────────
  if (verbose) message(sprintf("[RankShift] Computing statistics for %d genes...",
                               length(all_genes_found)))

  case_mat  <- rank_mat[, intersect(case_s, colnames(rank_mat)), drop = FALSE]
  ctrl_mat  <- rank_mat[, intersect(ctrl_s, colnames(rank_mat)), drop = FALSE]

  median_case <- apply(case_mat, 1, stats::median, na.rm = TRUE)
  median_ctrl <- apply(ctrl_mat, 1, stats::median, na.rm = TRUE)
  shift_score <- median_ctrl - median_case   # positive = moved UP in Case

  p_wilcox <- vapply(seq_len(nrow(rank_mat)), function(i) {
    cv <- case_mat[i, ]; cv <- cv[!is.na(cv)]
    rv <- ctrl_mat[i, ]; rv <- rv[!is.na(rv)]
    if (length(cv) < 2 || length(rv) < 2) return(1.0)
    tryCatch(stats::wilcox.test(cv, rv, exact = FALSE)$p.value,
             error = function(e) 1.0)
  }, numeric(1))
  fdr_wilcox <- .bh_adjust(p_wilcox)

  # Gene role annotation
  lig_genes <- unique(unlist(lr_db$ligand_genes))
  rec_genes <- unique(unlist(lr_db$receptor_genes))
  role <- vapply(all_genes_found, function(g) {
    is_l <- g %in% lig_genes
    is_r <- g %in% rec_genes
    if (is_l && is_r) "both" else if (is_l) "ligand" else if (is_r) "receptor" else "other"
  }, character(1))

  result <- data.frame(
    gene              = all_genes_found,
    role              = role,
    median_rank_case  = round(median_case, 4),
    median_rank_ctrl  = round(median_ctrl, 4),
    shift_score       = round(shift_score, 4),
    p_wilcox          = signif(p_wilcox, 4),
    fdr_wilcox        = signif(fdr_wilcox, 4),
    mean_expr_case    = round(expr_means$case[all_genes_found], 4),
    mean_expr_ctrl    = round(expr_means$ctrl[all_genes_found], 4),
    stringsAsFactors  = FALSE
  )
  result$log2fc_expr <- round(
    log2((result$mean_expr_case + 0.1) / (result$mean_expr_ctrl + 0.1)), 3)
  result$sig_stars <- .pval_stars(result$fdr_wilcox)

  result <- result[order(-result$shift_score), ]
  rownames(result) <- NULL

  class(result) <- c("RankShiftResult", "data.frame")
  attr(result, "case_label") <- case_label
  attr(result, "ctrl_label") <- ctrl_label

  n_up <- sum(result$shift_score > 0.1 & result$fdr_wilcox < 0.05, na.rm = TRUE)
  if (verbose) message(sprintf(
    "[RankShift] Done. %d genes significantly up-ranked in %s (|shift|>0.1, FDR<0.05).",
    n_up, case_label))

  result
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' Compute median normalized rank per gene across cells in one sample
#'
#' For each cell, all expressed genes are ranked within the cell's full
#' transcriptome (rank 1 = highest expressed) and normalized by the number of
#' expressed genes, so a smaller value means a more transcriptome-dominant gene.
#' Only the requested \code{target_genes} are returned, but the ranking is always
#' computed over the whole transcriptome so the score reflects each gene's
#' position among all expressed genes, not within a panel. Memory stays bounded
#' by \code{length(target_genes) x ncol(mat)}.
#' @keywords internal
.median_norm_rank_per_gene <- function(mat, min_detection_frac, target_genes = NULL) {
  keep <- if (is.null(target_genes)) rownames(mat) else intersect(target_genes, rownames(mat))
  if (length(keep) == 0) return(stats::setNames(numeric(0), character(0)))
  keep_idx <- match(keep, rownames(mat))
  n_genes <- nrow(mat)
  n_cells <- ncol(mat)

  # Store only the target genes' ranks; rank over the full transcriptome per cell.
  tgt_ranks <- matrix(NA_real_, length(keep), n_cells, dimnames = list(keep, colnames(mat)))
  for (j in seq_len(n_cells)) {
    col <- as.numeric(mat[, j])
    expressed <- col > 0
    n_exp <- sum(expressed)
    if (n_exp == 0) next
    # ties.method = "average"; rank(-x) makes rank 1 = highest expression.
    norm_r <- rep(NA_real_, n_genes)
    norm_r[expressed] <- rank(-col[expressed], ties.method = "average") / n_exp
    tgt_ranks[, j] <- norm_r[keep_idx]
  }

  # Per-gene median normalized rank over cells where it is expressed, then drop
  # genes detected in fewer than min_detection_frac of cells.
  detection <- rowMeans(mat[keep, , drop = FALSE] > 0, na.rm = TRUE)
  med_ranks <- apply(tgt_ranks, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    stats::median(x)
  })
  med_ranks[detection < min_detection_frac] <- NA_real_
  med_ranks[!is.na(med_ranks)]
}

#' Compute mean expression per gene per group
#' @keywords internal
.compute_expr_means <- function(sample_list, genes, group_info,
                                case_label, ctrl_label, layer) {
  case_s <- names(group_info)[group_info == case_label]
  ctrl_s <- names(group_info)[group_info == ctrl_label]

  accum <- function(sample_names) {
    mats <- lapply(sample_names, function(s) {
      m <- .extract_matrix(sample_list[[s]], layer = layer)
      keep <- intersect(genes, rownames(m))
      rowMeans(m[keep, , drop = FALSE], na.rm = TRUE)
    })
    # Average across the samples that actually contain each gene (dividing by the
    # total sample count would under-estimate genes absent from some samples).
    all_g <- Reduce(union, lapply(mats, names))
    sum_v <- setNames(rep(0.0, length(all_g)), all_g)
    cnt_v <- setNames(rep(0L, length(all_g)), all_g)
    for (mv in mats) {
      sum_v[names(mv)] <- sum_v[names(mv)] + mv
      cnt_v[names(mv)] <- cnt_v[names(mv)] + 1L
    }
    sum_v / pmax(cnt_v, 1L)
  }

  list(
    case = accum(intersect(case_s, names(sample_list))),
    ctrl = accum(intersect(ctrl_s, names(sample_list)))
  )
}

#' Print method for RankShiftResult
#' @param x RankShiftResult
#' @param n Rows to display. Default: 10.
#' @param ... ignored
#' @export
print.RankShiftResult <- function(x, n = 10, ...) {
  case_l <- attr(x, "case_label")
  ctrl_l <- attr(x, "ctrl_label")
  n_up   <- sum(x$shift_score > 0.1 & x$fdr_wilcox < 0.05, na.rm = TRUE)
  n_dn   <- sum(x$shift_score < -0.1 & x$fdr_wilcox < 0.05, na.rm = TRUE)
  cat(sprintf("RankShiftResult: %d genes | Up in %s: %d | Down: %d\n",
              nrow(x), case_l, n_up, n_dn))
  cat("Top up-ranked genes in", case_l, ":\n")
  show_cols <- c("gene","role","shift_score","log2fc_expr","fdr_wilcox","sig_stars")
  show_cols <- intersect(show_cols, names(x))
  tmp <- as.data.frame(x)[x$shift_score > 0, show_cols, drop = FALSE]
  print(utils::head(tmp, n), row.names = FALSE)
  invisible(x)
}
