# R/plot_functions.R

#' Bubble Plot of Logic Consensus Scores (Case vs Control)
#'
#' Visualizes L-R pair enrichment between two groups as a bubble plot, where
#' x-axis = Ctrl frequency, y-axis = Case frequency, bubble size = |asymmetry|,
#' and color = log2 fold-change in mean LCS.
#'
#' @param result A \code{LogicCommResult} from \code{\link{CompareLogicGroups}}.
#' @param top_n Number of top pairs to label. Default: \code{20}.
#' @param fdr_cutoff FDR cutoff for labeling. Default: \code{0.05}.
#' @param min_asymmetry Minimum asymmetry for labeling. Default: \code{0.2}.
#' @param color_by \code{"log2fc"} or \code{"asymmetry"}. Default: \code{"log2fc"}.
#' @param pathway_filter Optional character vector of pathways to highlight.
#' @param point_alpha Point transparency. Default: \code{0.75}.
#' @param label_size ggrepel label size. Default: \code{3}.
#' @param title Plot title. Default: auto-generated.
#'
#' @return A ggplot2 object.
#' @importFrom ggplot2 ggplot aes geom_point geom_abline scale_color_gradient2 scale_size_continuous labs theme_bw theme element_text
#' @importFrom ggrepel geom_text_repel
#' @export
plot_lcs_bubble <- function(result,
                             top_n          = 20,
                             fdr_cutoff     = 0.05,
                             min_asymmetry  = 0.2,
                             color_by       = "log2fc",
                             pathway_filter = NULL,
                             point_alpha    = 0.75,
                             label_size     = 3,
                             title          = NULL) {

  stopifnot(inherits(result, "LogicCommResult"))
  color_by <- match.arg(color_by, c("log2fc", "asymmetry"))
  case_l <- attr(result, "case_label") %||% "Case"
  ctrl_l <- attr(result, "ctrl_label") %||% "Ctrl"

  df <- result
  df <- df[!is.na(df$case_freq) & !is.na(df$ctrl_freq), ]

  # Color variable
  df$color_val <- if (color_by == "log2fc") df$log2fc_lcs else df$asymmetry
  col_label    <- if (color_by == "log2fc") "log2(LCS FC)" else "Asymmetry"

  # Label candidates
  df$do_label <- df$fdr_fisher <= fdr_cutoff & abs(df$asymmetry) >= min_asymmetry
  if (!is.null(pathway_filter)) {
    df$do_label <- df$do_label & (df$pathway %in% pathway_filter)
  }
  # Limit to top_n
  df_label <- df[df$do_label, , drop = FALSE]
  if (nrow(df_label) > top_n) {
    df_label <- df_label[order(-abs(df_label$asymmetry)), ][seq_len(top_n), ]
  }

  # Point size = |asymmetry|, clamped
  df$size_val <- pmin(abs(df$asymmetry) + 0.02, 1)

  if (is.null(title)) {
    title <- sprintf("L-R Logic Consensus: %s vs %s", case_l, ctrl_l)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(
      x = ctrl_freq, y = case_freq,
      size = size_val, color = color_val)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "grey70") +
    ggplot2::geom_point(alpha = point_alpha) +
    ggrepel::geom_text_repel(
      data = df_label,
      ggplot2::aes(label = lr_pair),
      size = label_size, max.overlaps = 40,
      fontface = "bold",
      segment.color = "grey50", segment.size = 0.3) +
    scale_color_logiccomm_diverging(midpoint = 0, name = col_label) +
    ggplot2::scale_size_continuous(range = c(1, 8), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-0.02, 1.05),
                                labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_y_continuous(limits = c(-0.02, 1.05),
                                labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      title = title,
      x     = paste0(ctrl_l, " frequency"),
      y     = paste0(case_l, " frequency"),
      caption = sprintf("Bubble size = |asymmetry|; labeled: FDR<%.2f & |asymm|>%.1f",
                        fdr_cutoff, min_asymmetry)) +
    theme_logiccomm() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank())

  p
}

#' Heatmap of LCS Across All Samples
#'
#' Displays a genes x samples heatmap of raw LCS values for the top
#' differentially active L-R pairs.
#'
#' @param result A \code{LogicCommResult} from \code{\link{CompareLogicGroups}}.
#' @param group_info Named character vector mapping sample names to group labels.
#' @param top_n Number of top L-R pairs (by |asymmetry|) to show. Default: \code{30}.
#' @param fdr_cutoff FDR threshold to filter pairs. Default: \code{0.1}.
#' @param scale_rows Logical: row-scale the LCS values (Z-score). Default: \code{TRUE}.
#' @param color_palette Color scale. Default: blue-white-red.
#' @param ... Additional arguments passed to \code{pheatmap::pheatmap}.
#'
#' @return Invisibly returns the pheatmap object.
#' @export
plot_lcs_heatmap <- function(result,
                              group_info,
                              top_n        = 30,
                              fdr_cutoff   = 0.1,
                              scale_rows   = TRUE,
                              color_palette = NULL,
                              ...) {

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' required for heatmap. Install with: install.packages('pheatmap')")
  }

  stopifnot(inherits(result, "LogicCommResult"))
  lcs_mat <- attr(result, "lcs_mat")
  if (is.null(lcs_mat)) stop("LCS matrix not found in result attributes.")

  case_l <- attr(result, "case_label") %||% "Case"
  ctrl_l <- attr(result, "ctrl_label") %||% "Ctrl"

  # Select top pairs
  sig <- result[!is.na(result$fdr_fisher) & result$fdr_fisher <= fdr_cutoff, ]
  sig <- sig[order(-abs(sig$asymmetry)), ]
  sig <- utils::head(sig, top_n)
  top_pairs <- sig$lr_pair

  if (length(top_pairs) == 0) {
    stop("No pairs pass fdr_cutoff=", fdr_cutoff, ". Try relaxing the threshold.")
  }

  mat <- lcs_mat[top_pairs, , drop = FALSE]
  mat[is.na(mat)] <- 0

  # Row scale
  if (scale_rows && nrow(mat) > 1) {
    mat <- t(scale(t(mat)))
    mat[is.nan(mat)] <- 0
  }

  # Annotation
  if (is.null(names(group_info))) names(group_info) <- colnames(mat)
  col_ann <- data.frame(Group = group_info[colnames(mat)], row.names = colnames(mat))

  ann_colors <- list(Group = c(
    setNames("firebrick", case_l),
    setNames("steelblue", ctrl_l)
  ))

  if (is.null(color_palette)) {
    color_palette <- grDevices::colorRampPalette(
      c("steelblue","white","firebrick"))(100)
  }

  pheatmap::pheatmap(mat,
    annotation_col = col_ann,
    annotation_colors = ann_colors,
    color = color_palette,
    cluster_cols = FALSE,
    fontsize_row = 8,
    fontsize_col = 7,
    main = sprintf("LCS Heatmap: %s vs %s (top %d pairs)", case_l, ctrl_l, length(top_pairs)),
    ...)
}

#' UMAP Plot Colored by Single L-R Logic Activity
#'
#' Colors each cell on the UMAP by whether it is an "active sender" (ligand > anchor),
#' "active receiver" (receptor > anchor), "both", or "neither" for a given L-R pair.
#'
#' @param reo_mat REO binary matrix (genes x cells).
#' @param seurat_obj Seurat object with UMAP coordinates in \code{reductions$umap}.
#' @param lr_pair Character. L-R pair name (must match \code{lr_db$lr_pair}).
#' @param lr_db LR database. Default: \code{lr_pairs_human}.
#' @param pt_size Point size. Default: \code{0.8}.
#' @param colors Named vector with colors for categories. Default: traffic-light scheme.
#'
#' @return A ggplot2 object.
#' @export
plot_umap_logic <- function(reo_mat,
                             seurat_obj,
                             lr_pair,
                             lr_db    = lr_pairs_human,
                             pt_size  = 0.8,
                             colors   = NULL) {

  if (!requireNamespace("Seurat", quietly = TRUE) &&
      !requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Seurat package required.")
  }

  # Get UMAP coords
  umap_emb <- tryCatch(
    seurat_obj@reductions$umap@cell.embeddings[, 1:2],
    error = function(e) stop("No UMAP found. Run RunUMAP() first."))
  colnames(umap_emb) <- c("UMAP1","UMAP2")

  # Find LR entry
  lr_row <- lr_db[lr_db$lr_pair == lr_pair, ]
  if (nrow(lr_row) == 0) stop("lr_pair '", lr_pair, "' not found in lr_db.")

  lig_genes <- lr_row$ligand_genes[[1]]
  rec_genes <- lr_row$receptor_genes[[1]]

  common_cells <- intersect(rownames(umap_emb), colnames(reo_mat))
  reo_sub <- reo_mat[, common_cells, drop = FALSE]

  lig_logic <- .resolve_complex_logic(lig_genes, reo_sub)
  rec_logic <- .resolve_complex_logic(rec_genes, reo_sub)

  plot_df <- data.frame(
    UMAP1 = umap_emb[common_cells, 1],
    UMAP2 = umap_emb[common_cells, 2],
    status = dplyr::case_when(
      lig_logic & rec_logic  ~ "Sender & Receiver",
      lig_logic & !rec_logic ~ "Sender (Ligand)",
      !lig_logic & rec_logic ~ "Receiver (Receptor)",
      TRUE                   ~ "Inactive"
    ),
    stringsAsFactors = FALSE
  )
  plot_df$status <- factor(plot_df$status,
    levels = c("Inactive","Sender (Ligand)","Receiver (Receptor)","Sender & Receiver"))

  if (is.null(colors)) {
    colors <- c("Inactive"             = "grey88",
                "Sender (Ligand)"      = "#1F6F78",
                "Receiver (Receptor)"  = "#E8A33D",
                "Sender & Receiver"    = "#D1495B")
  }

  ggplot2::ggplot(plot_df[order(plot_df$status), ],
    ggplot2::aes(x = UMAP1, y = UMAP2, color = status)) +
    ggplot2::geom_point(size = pt_size, alpha = 0.6) +
    ggplot2::scale_color_manual(values = colors, name = "Logic Status") +
    ggplot2::labs(title = lr_pair,
                  subtitle = paste0("Ligand: ", lr_row$ligand[1],
                                    "  |  Receptor: ", lr_row$receptor[1])) +
    theme_logiccomm() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "grey80")
    )
}
