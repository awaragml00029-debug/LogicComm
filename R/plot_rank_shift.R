# R/plot_rank_shift.R

#' Rank-Shift Volcano Plot
#'
#' A modified volcano plot where the x-axis shows gene expression rank shift
#' (how much a gene's within-cell dominance changed in Case vs Control) and
#' the y-axis shows statistical significance. Unlike a standard DE volcano,
#' this captures genes that became **transcriptionally dominant** in Case cells
#' without necessarily having a large absolute fold-change.
#'
#' @param rs_result A \code{RankShiftResult} from \code{\link{calc_rank_shift}}.
#' @param top_n_label Number of top genes (by |shift_score|) to label.
#'   Default: \code{20}.
#' @param fdr_cutoff FDR significance line. Default: \code{0.05}.
#' @param min_shift_label Minimum |shift_score| to label a gene. Default: \code{0.1}.
#' @param color_by_role Logical. Color points by gene role (ligand/receptor/both).
#'   Default: \code{TRUE}.
#' @param highlight_genes Optional character vector of specific genes to label
#'   regardless of significance thresholds.
#' @param point_size Numeric. Base point size. Default: \code{2}.
#' @param label_size Numeric. ggrepel label font size. Default: \code{3}.
#' @param title Plot title. Default: auto.
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline scale_color_manual labs theme_bw theme element_text annotate
#' @importFrom ggrepel geom_text_repel
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
#' fun <- get("plot_rank_shift")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
plot_rank_shift <- function(rs_result,
                             top_n_label      = 20,
                             fdr_cutoff       = 0.05,
                             min_shift_label  = 0.1,
                             color_by_role    = TRUE,
                             highlight_genes  = NULL,
                             point_size       = 2,
                             label_size       = 3,
                             title            = NULL) {

  stopifnot(inherits(rs_result, "RankShiftResult"))

  case_l <- attr(rs_result, "case_label") %||% "Case"
  ctrl_l <- attr(rs_result, "ctrl_label") %||% "Ctrl"

  df <- rs_result
  df <- df[!is.na(df$p_wilcox), ]
  df$neg_log10_fdr <- -log10(pmax(df$fdr_wilcox, 1e-10))

  # Significance category
  df$sig_cat <- dplyr::case_when(
    df$fdr_wilcox < fdr_cutoff & df$shift_score >  min_shift_label ~ paste0("Up in ", case_l),
    df$fdr_wilcox < fdr_cutoff & df$shift_score < -min_shift_label ~ paste0("Up in ", ctrl_l),
    TRUE                                                             ~ "ns"
  )

  # Role colors
  if (color_by_role) {
    df$point_color <- interaction(df$sig_cat, df$role, sep = "|")
    # Simplify: use sig_cat for color (role shown by shape/size if needed)
    df$color_var <- df$sig_cat
    simple_colors <- c(
      setNames("#D1495B", paste0("Up in ", case_l)),
      setNames("#2E86AB", paste0("Up in ", ctrl_l)),
      "ns" = "grey78"
    )
    df$point_size_var <- ifelse(df$sig_cat == "ns", point_size * 0.7, point_size * 1.2)
    df$alpha_var <- ifelse(df$sig_cat == "ns", 0.4, 0.85)
  } else {
    df$color_var <- df$sig_cat
    simple_colors <- c(
      setNames("#D1495B", paste0("Up in ", case_l)),
      setNames("#2E86AB", paste0("Up in ", ctrl_l)),
      "ns" = "grey78"
    )
    df$point_size_var <- point_size
    df$alpha_var <- 0.7
  }

  # Label selection
  sig_df <- df[df$sig_cat != "ns" & abs(df$shift_score) >= min_shift_label, ]
  sig_df <- sig_df[order(-abs(sig_df$shift_score)), ]
  if (nrow(sig_df) > top_n_label) sig_df <- sig_df[seq_len(top_n_label), ]

  if (!is.null(highlight_genes)) {
    extra <- df[df$gene %in% highlight_genes & !df$gene %in% sig_df$gene, ]
    sig_df <- rbind(sig_df, extra)
  }

  if (is.null(title)) {
    title <- sprintf("Gene Rank Shift: %s vs %s", case_l, ctrl_l)
  }

  fdr_line <- -log10(fdr_cutoff)

  # Count up/down for annotation
  n_up <- sum(df$sig_cat == paste0("Up in ", case_l))
  n_dn <- sum(df$sig_cat == paste0("Up in ", ctrl_l))

  p <- ggplot2::ggplot(df, ggplot2::aes(
      x = shift_score, y = neg_log10_fdr, color = color_var)) +
    ggplot2::geom_point(size = df$point_size_var, alpha = df$alpha_var) +
    ggplot2::geom_hline(yintercept = fdr_line,
                        linetype = "dashed", color = "grey40", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = c(-min_shift_label, min_shift_label),
                        linetype = "dotted", color = "grey60", linewidth = 0.4) +
    ggrepel::geom_text_repel(
      data = sig_df,
      ggplot2::aes(label = gene),
      size = label_size, max.overlaps = 50,
      segment.color = "grey40", segment.size = 0.3,
      box.padding = 0.3) +
    ggplot2::scale_color_manual(values = simple_colors, name = "Direction") +
    ggplot2::annotate("text", x = max(df$shift_score, na.rm=TRUE) * 0.85,
                      y = max(df$neg_log10_fdr, na.rm=TRUE) * 0.95,
                      label = sprintf("n=%d", n_up), color = "#E53935",
                      size = 4, fontface = "bold") +
    ggplot2::annotate("text", x = min(df$shift_score, na.rm=TRUE) * 0.85,
                      y = max(df$neg_log10_fdr, na.rm=TRUE) * 0.95,
                      label = sprintf("n=%d", n_dn), color = "#1E88E5",
                      size = 4, fontface = "bold") +
    ggplot2::labs(
      title    = title,
      subtitle = sprintf("Rank shift > 0: gene became more dominant in %s cells", case_l),
      x        = sprintf("Rank Shift Score\n(<- dominant in %s   dominant in %s ->)",
                         ctrl_l, case_l),
      y        = "-log10(FDR)",
      caption  = sprintf("FDR cutoff = %.2f | |shift| threshold = %.2f",
                         fdr_cutoff, min_shift_label)) +
    theme_logiccomm() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey50", size = 10),
      legend.position = "top")

  p
}

#' Dot Plot of Top Rank-Shifted L-R Pairs
#'
#' Combines rank shift and LCS comparison to show the top differentially active
#' L-R pairs alongside the rank shift of their component genes. Useful for
#' identifying pairs where both the expression dominance AND the communication
#' logic have shifted.
#'
#' @param rs_result A \code{RankShiftResult}.
#' @param lcs_result A \code{LogicCommResult} from \code{\link{CompareLogicGroups}}.
#' @param top_n Top N pairs (by |asymmetry|) to display. Default: \code{20}.
#' @param fdr_cutoff FDR cutoff for LCS significance. Default: \code{0.1}.
#'
#' @return A ggplot2 object.
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
#' fun <- get("plot_shift_lcs_combined")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
plot_shift_lcs_combined <- function(rs_result,
                                     lcs_result,
                                     top_n     = 20,
                                     fdr_cutoff = 0.1) {
  stopifnot(inherits(rs_result, "RankShiftResult"))
  stopifnot(inherits(lcs_result, "LogicCommResult"))

  case_l <- attr(lcs_result, "case_label") %||% "Case"

  # Select top LCS pairs
  sig_lcs <- lcs_result[!is.na(lcs_result$fdr_fisher) &
                         lcs_result$fdr_fisher <= fdr_cutoff, ]
  sig_lcs <- utils::head(sig_lcs[order(-abs(sig_lcs$asymmetry)), ], top_n)

  if (nrow(sig_lcs) == 0) {
    stop("No significant LCS pairs found. Relax fdr_cutoff.")
  }

  # For each pair, get ligand/receptor rank shift
  rs_lookup <- setNames(rs_result$shift_score, rs_result$gene)

  rows <- lapply(seq_len(nrow(sig_lcs)), function(i) {
    row    <- sig_lcs[i, ]
    lig    <- gsub("\\+.*", "", row$ligand)
    rec    <- gsub("\\+.*", "", row$receptor)
    data.frame(
      lr_pair       = row$lr_pair,
      pathway       = if ("pathway" %in% names(row)) row$pathway else NA,
      asymmetry     = row$asymmetry,
      ligand        = lig,
      receptor      = rec,
      lig_shift     = rs_lookup[lig] %||% NA_real_,
      rec_shift     = rs_lookup[rec] %||% NA_real_,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  df$combined_shift <- rowMeans(
    cbind(df$lig_shift, df$rec_shift), na.rm = TRUE)
  df$lr_pair <- factor(df$lr_pair, levels = rev(df$lr_pair))

  ggplot2::ggplot(df, ggplot2::aes(y = lr_pair)) +
    ggplot2::geom_point(ggplot2::aes(x = lig_shift, color = "Ligand shift"),
                        size = 4, shape = 16) +
    ggplot2::geom_point(ggplot2::aes(x = rec_shift, color = "Receptor shift"),
                        size = 4, shape = 17) +
    ggplot2::geom_segment(ggplot2::aes(x = lig_shift, xend = rec_shift,
                                       yend = lr_pair),
                          linewidth = 0.6, color = "grey60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::scale_color_manual(
      values = c("Ligand shift" = "#D1495B", "Receptor shift" = "#2E86AB"),
      name = "") +
    ggplot2::labs(
      title = sprintf("Rank Shift + LCS: Top %d Pairs in %s", top_n, case_l),
      x     = "Rank Shift Score (-> more dominant in Case)",
      y     = NULL) +
    theme_logiccomm(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "top",
      axis.text.y = ggplot2::element_text(size = 9))
}
