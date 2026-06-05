# R/celltype_plots.R

# Filter an lr_table to requested sender/receiver cell types, erroring with the
# list of valid values when a requested type is not present at all. This turns a
# silent empty plot (a common cause of "the plot is broken" reports when a
# cell-type name is mistyped) into an actionable message.
.filter_celltype_axis <- function(df, senders = NULL, receivers = NULL) {
  if (!is.null(senders)) {
    bad <- setdiff(senders, unique(df$sender_type))
    if (length(bad) == length(senders)) {
      stop("None of the requested sender(s) are cell types in this object: ",
           paste(bad, collapse = ", "), ".\nAvailable sender types: ",
           paste(sort(unique(df$sender_type)), collapse = ", "), call. = FALSE)
    }
    df <- df[df$sender_type %in% senders, , drop = FALSE]
  }
  if (!is.null(receivers)) {
    bad <- setdiff(receivers, unique(df$receiver_type))
    if (length(bad) == length(receivers)) {
      stop("None of the requested receiver(s) are cell types in this object: ",
           paste(bad, collapse = ", "), ".\nAvailable receiver types: ",
           paste(sort(unique(df$receiver_type)), collapse = ", "), call. = FALSE)
    }
    df <- df[df$receiver_type %in% receivers, , drop = FALSE]
  }
  df
}

# Build an informative message when a cell-type axis exists but carries no active
# events, so users can tell "no signal here" apart from "you mistyped a name".
.no_active_events_message <- function(ct_comm, senders = NULL, receivers = NULL) {
  scope <- c(
    if (!is.null(senders)) paste0("sender(s) ", paste(senders, collapse = ", ")),
    if (!is.null(receivers)) paste0("receiver(s) ", paste(receivers, collapse = ", "))
  )
  scope_txt <- if (length(scope)) paste0(" for ", paste(scope, collapse = " and ")) else ""
  paste0(
    "No active L-R events to plot", scope_txt, ".\n",
    "The requested cell type(s) exist but have no events passing the active-call ",
    "gates (lcs_threshold, min_active_edges, min_expr_frac). Try a different ",
    "sender/receiver, lower min_expr_frac/lcs_threshold in ",
    "summarize_celltype_communication(), or inspect ct_comm$lr_table directly."
  )
}

#' Plot Cell-Type-Level Communication Network
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param metric Metric for edge weights.
#' @param layout Node layout: "auto" (force-directed via igraph when available,
#'   otherwise a circle), "fr", or "circle".
#' @param min_weight Filter edges below this threshold.
#' @param arrow_size Size of edge arrows.
#' @param show_labels Show cell-type names.
#' @param color_edges_by How to color edges: \code{"top_pathway"} or \code{"range"}.
#' @return A ggplot2 object.
#' @export
plot_celltype_network <- function(ct_comm,
                                  metric = c("sum_lcs", "n_active_lr"),
                                  layout = c("auto", "fr", "circle"),
                                  min_weight = 0.05,
                                  arrow_size = 0.2,
                                  show_labels = TRUE,
                                  color_edges_by = c("top_pathway", "range")) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  metric <- match.arg(metric)
  layout <- match.arg(layout)
  color_edges_by <- match.arg(color_edges_by)
  pair_sum <- ct_comm$pair_summary
  if (!metric %in% names(pair_sum)) stop("metric not found in pair_summary: ", metric)
  pair_sum <- pair_sum[is.finite(pair_sum[[metric]]) & pair_sum[[metric]] >= min_weight, , drop = FALSE]
  if (nrow(pair_sum) == 0) stop("No edges above min_weight.")

  nodes <- ct_comm$role_summary
  coords <- .celltype_graph_layout(nodes$cell_type, pair_sum, layout = layout)
  nodes$x <- coords$x[match(nodes$cell_type, coords$cell_type)]
  nodes$y <- coords$y[match(nodes$cell_type, coords$cell_type)]
  nodes$node_label <- .short_label(as.character(nodes$cell_type), 18L)
  edges <- merge(pair_sum, nodes[, c("cell_type", "x", "y", "hub_score", "sender_receiver_balance")],
                 by.x = "sender_type", by.y = "cell_type", all.x = TRUE)
  edges <- merge(edges, nodes[, c("cell_type", "x", "y")],
                 by.x = "receiver_type", by.y = "cell_type", all.x = TRUE, suffixes = c("", "end"))
  if (!"communication_support_label" %in% names(edges)) {
    edges$communication_support_label <- ifelse(edges$n_active_lr > 0 & edges$n_local_active == 0 & edges$n_distal_candidate > 0,
                                                "global_only_candidate",
                                                ifelse(edges$n_distal_candidate > 0, "mixed_local_global", "local_graph_supported"))
  }
  edges$edge_color <- if (color_edges_by == "range") edges$dominant_communication_range else edges$top_pathway
  edges$edge_linetype <- ifelse(edges$dominant_communication_range == "distal/endocrine", "distal/endocrine", "local/mixed")
  edges$edge_support_alpha <- ifelse(edges$communication_support_label == "global_only_candidate", 0.35,
                                     ifelse(edges$communication_support_label == "mixed_local_global", 0.6, 0.8))
  edges$value <- edges[[metric]]
  self_edge <- is.finite(edges$x) & is.finite(edges$y) & is.finite(edges$xend) & is.finite(edges$yend) &
    abs(edges$x - edges$xend) < .Machine$double.eps & abs(edges$y - edges$yend) < .Machine$double.eps
  if (any(self_edge)) {
    tx <- -edges$y[self_edge]
    ty <- edges$x[self_edge]
    norm <- sqrt(tx^2 + ty^2)
    norm[norm == 0] <- 1
    tx <- tx / norm
    ty <- ty / norm
    rx <- edges$x[self_edge]
    ry <- edges$y[self_edge]
    rnorm <- sqrt(rx^2 + ry^2)
    rnorm[rnorm == 0] <- 1
    rx <- rx / rnorm
    ry <- ry / rnorm
    loop_radius <- 0.14
    inward_shift <- 0.05
    edges$x[self_edge] <- edges$x[self_edge] + loop_radius * tx - inward_shift * rx
    edges$y[self_edge] <- edges$y[self_edge] + loop_radius * ty - inward_shift * ry
    edges$xend[self_edge] <- edges$xend[self_edge] - loop_radius * tx - inward_shift * rx
    edges$yend[self_edge] <- edges$yend[self_edge] - loop_radius * ty - inward_shift * ry
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_curve(data = edges,
                        ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                                     linewidth = value, color = edge_color, linetype = edge_linetype,
                                     alpha = edge_support_alpha),
                        curvature = 0.18,
                        arrow = grid::arrow(length = grid::unit(arrow_size, "cm"), type = "closed"),
                        lineend = "round") +
    ggplot2::geom_point(data = nodes,
                        ggplot2::aes(x = x, y = y, fill = sender_receiver_balance, size = hub_score),
                        shape = 21, color = "white", stroke = 0.6) +
    ggplot2::scale_linewidth_continuous(range = c(0.3, 2.5), name = metric) +
    scale_fill_logiccomm_diverging(midpoint = 0, name = "S-R Balance") +
    ggplot2::scale_size_continuous(range = c(3, 12), name = "Hub score") +
    ggplot2::scale_linetype_manual(values = c("local/mixed" = "solid", "distal/endocrine" = "dashed"), name = "Range style") +
    ggplot2::scale_alpha_identity() +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "Cell-type communication network", color = if (color_edges_by == "range") "Range" else "Top pathway",
                  caption = "Faded dashed edges are global-only distal candidates without local active edge support.") +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.25), colour = logiccomm_brand$ink),
      plot.caption = ggplot2::element_text(colour = "grey45", size = ggplot2::rel(0.78), hjust = 0),
      legend.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.9)),
      legend.position = "right")

  if (color_edges_by == "range") {
    p <- p + ggplot2::scale_color_manual(values = .logiccomm_palettes$range, na.value = "grey70")
  } else {
    p <- p + scale_color_logiccomm_d(na.value = "grey75")
  }
  if (show_labels) {
    p <- p + ggrepel::geom_text_repel(data = nodes, ggplot2::aes(x = x, y = y, label = node_label),
                                      size = 3.5, fontface = "bold", max.overlaps = 20)
  }
  p
}

#' Plot a Cell-Type Pair Heatmap
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param metric Pair-summary metric to display.
#' @param title Optional title.
#' @param rotate_x Rotation angle for x-axis labels.
#' @return A ggplot2 object.
#' @export
plot_celltype_heatmap <- function(ct_comm,
                                  metric = c("n_active_lr", "active_lr_event_count", "sum_lcs",
                                             "mean_lcs_active", "sum_lcs_all", "n_edges", "edge_weight_sum",
                                             "sum_active_edges", "sum_active_edge_weight"),
                                  title = NULL,
                                  rotate_x = 45) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  metric <- match.arg(metric)
  df <- ct_comm$pair_summary
  if (!metric %in% names(df)) stop("metric not found in pair_summary: ", metric)
  if (is.null(title)) title <- paste("Cell-type communication:", metric)
  ggplot2::ggplot(df, ggplot2::aes(x = sender_type, y = receiver_type, fill = .data[[metric]])) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_logiccomm_c(name = metric, na.value = "grey90") +
    ggplot2::labs(title = title, x = "Sender type", y = "Receiver type") +
    theme_logiccomm() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = rotate_x, hjust = 1))
}

#' Plot Communication Range Summary
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @return A ggplot2 object.
#' @export
plot_communication_range_summary <- function(ct_comm) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  df <- ct_comm$lr_table
  df <- df[df$active %in% TRUE, , drop = FALSE]
  if (nrow(df) == 0) stop("No active events to plot.")
  pw_df <- as.data.frame(table(df$pathway, df$communication_range), stringsAsFactors = FALSE)
  colnames(pw_df) <- c("Pathway", "Range", "Count")
  pw_totals <- tapply(pw_df$Count, pw_df$Pathway, sum)
  pw_df$Pathway <- factor(pw_df$Pathway, levels = names(sort(pw_totals)))
  ggplot2::ggplot(pw_df, ggplot2::aes(x = Pathway, y = Count, fill = Range)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = .logiccomm_palettes$range, na.value = "grey80") +
    ggplot2::coord_flip() +
    theme_logiccomm() +
    ggplot2::labs(title = "Signaling range by pathway", y = "Number of active L-R pairs", x = "Pathway")
}

#' Advanced Bubble Plot of L-R Logic Events
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param senders Optional sender filter.
#' @param receivers Optional receiver filter.
#' @param top_n_per_pair Number of rows per sender-receiver pair.
#' @return A ggplot2 object.
#' @export
plot_lr_bubble_advanced <- function(ct_comm, senders = NULL, receivers = NULL, top_n_per_pair = 5) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  df <- ct_comm$lr_table
  df <- .filter_celltype_axis(df, senders, receivers)
  df <- df[df$active %in% TRUE & is.finite(df$lcs), , drop = FALSE]
  if (nrow(df) == 0) {
    stop(.no_active_events_message(ct_comm, senders, receivers), call. = FALSE)
  }
  keys <- paste(df$sender_type, df$receiver_type, sep = "|||")
  keep <- unlist(lapply(split(seq_len(nrow(df)), keys), function(ii) ii[order(df$lcs[ii], decreasing = TRUE)][seq_len(min(top_n_per_pair, length(ii)))]), use.names = FALSE)
  df <- df[keep, , drop = FALSE]
  df$point_alpha <- ifelse(df$communication_range == "distal/endocrine" & df$n_active_neighborhood == 0, 0.45, 0.8)
  ggplot2::ggplot(df, ggplot2::aes(x = sender_type, y = receiver_type, size = lcs, color = communication_range, alpha = point_alpha)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~pathway, scales = "free") +
    ggplot2::scale_color_manual(values = .logiccomm_palettes$range, na.value = "grey70") +
    ggplot2::scale_size_continuous(range = c(2, 8)) +
    ggplot2::scale_alpha_identity() +
    theme_logiccomm() +
    ggplot2::labs(title = "Cell-type LogicComm hotspots", x = "Sender type", y = "Receiver type",
                  caption = "Faded distal bubbles are global-only candidates without local active edge support.") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot L-R Bubbles by Cell Type
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param sender Optional sender filter.
#' @param receiver Optional receiver filter.
#' @param top_n Number of rows to show.
#' @param active_only Whether to show active candidate rows only.
#' @param color_by Field used for color.
#' @param title Optional title.
#' @return A ggplot2 object.
#' @export
plot_lr_bubble_by_celltype <- function(ct_comm,
                                       sender = NULL,
                                       receiver = NULL,
                                       top_n = 30,
                                       active_only = TRUE,
                                       color_by = c("pathway", "lcs", "n_active_edges", "ligand_active_frac_sender", "receptor_active_frac_receiver", "communication_range"),
                                       title = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  color_by <- match.arg(color_by)
  df <- ct_comm$lr_table
  df <- .filter_celltype_axis(df, sender, receiver)
  if (isTRUE(active_only)) df <- df[df$active %in% TRUE, , drop = FALSE]
  df <- df[is.finite(df$lcs), , drop = FALSE]
  if (nrow(df) == 0) stop(.no_active_events_message(ct_comm, sender, receiver), call. = FALSE)
  df <- df[order(df$lcs, decreasing = TRUE), , drop = FALSE]
  df <- utils::head(df, top_n)
  df$axis_label <- paste(df$sender_type, df$receiver_type, sep = " -> ")
  df$lr_pair <- factor(df$lr_pair, levels = rev(unique(df$lr_pair)))
  if (is.null(title)) title <- "Top cell-type L-R communication events"
  p <- ggplot2::ggplot(df, ggplot2::aes(x = axis_label, y = lr_pair, size = lcs, color = .data[[color_by]])) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_size_continuous(range = c(2, 8), name = "LCS") +
    ggplot2::labs(title = title, x = "Cell-type pair", y = "L-R pair", color = color_by) +
    theme_logiccomm() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (is.numeric(df[[color_by]])) p + scale_color_logiccomm_c(name = color_by) else p
}

#' Plot Pathway Heatmap
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param metric Pathway-summary metric.
#' @param sender Optional sender filter.
#' @param receiver Optional receiver filter.
#' @param top_n_pathways Number of pathways to display.
#' @param title Optional title.
#' @param cluster_rows Cluster pathway rows using the selected metric across cell-type pairs.
#' @param cluster_cols Cluster sender-receiver pair columns using the selected metric across pathways.
#' @return A ggplot2 object.
#' @export
plot_pathway_heatmap <- function(ct_comm,
                                 metric = c("sum_lcs", "n_active_lr", "mean_lcs_active"),
                                 sender = NULL,
                                 receiver = NULL,
                                 top_n_pathways = 30,
                                 title = NULL,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  metric <- match.arg(metric)
  df <- ct_comm$pathway_summary
  if (!metric %in% names(df)) stop("metric not found in pathway_summary: ", metric)
  if (!is.null(sender)) df <- df[df$sender_type %in% sender, , drop = FALSE]
  if (!is.null(receiver)) df <- df[df$receiver_type %in% receiver, , drop = FALSE]
  if (nrow(df) == 0) stop("No pathway rows to plot.")
  totals <- tapply(df[[metric]], df$pathway, sum, na.rm = TRUE)
  keep <- names(sort(totals, decreasing = TRUE))[seq_len(min(top_n_pathways, length(totals)))]
  df <- df[df$pathway %in% keep, , drop = FALSE]
  df$celltype_pair <- paste(df$sender_type, df$receiver_type, sep = " -> ")
  if (isTRUE(cluster_rows) || isTRUE(cluster_cols)) {
    cluster_values <- df[[metric]]
    cluster_values[!is.finite(cluster_values) | is.na(cluster_values)] <- 0
    cluster_df <- stats::aggregate(
      cluster_values,
      by = list(pathway = df$pathway, celltype_pair = df$celltype_pair),
      FUN = sum,
      na.rm = TRUE
    )
    mat <- stats::xtabs(x ~ pathway + celltype_pair, data = cluster_df)
    if (isTRUE(cluster_rows) && nrow(mat) > 1) {
      ord <- tryCatch(stats::hclust(stats::dist(mat))$order, error = function(e) seq_len(nrow(mat)))
      df$pathway <- factor(df$pathway, levels = rownames(mat)[ord])
    }
    if (isTRUE(cluster_cols) && ncol(mat) > 1) {
      ord <- tryCatch(stats::hclust(stats::dist(t(mat)))$order, error = function(e) seq_len(ncol(mat)))
      df$celltype_pair <- factor(df$celltype_pair, levels = colnames(mat)[ord])
    }
  }
  if (is.null(title)) title <- paste("Pathway communication:", metric)
  ggplot2::ggplot(df, ggplot2::aes(x = celltype_pair, y = pathway, fill = .data[[metric]])) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_logiccomm_c(name = metric, na.value = "grey90") +
    ggplot2::labs(title = title, x = "Cell-type pair", y = "Pathway") +
    theme_logiccomm() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Explain a Cell-Type L-R Interaction
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param sender Sender cell type.
#' @param receiver Receiver cell type.
#' @param lr_pair L-R pair.
#' @param reo_mat Optional REO matrix.
#' @param seurat_obj Optional Seurat object.
#' @param lr_db Optional LR database.
#' @param cell_labels Optional cell labels.
#' @param pt_size Point size for optional plots.
#' @return A list with query, evidence, interpretation, and plot.
#' @export
explain_celltype_interaction <- function(ct_comm, sender, receiver, lr_pair,
                                         reo_mat = NULL, seurat_obj = NULL,
                                         lr_db = NULL, cell_labels = NULL, pt_size = 0.7) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  ev <- ct_comm$lr_table[ct_comm$lr_table$sender_type == sender &
                           ct_comm$lr_table$receiver_type == receiver &
                           ct_comm$lr_table$lr_pair == lr_pair, , drop = FALSE]
  if (nrow(ev) == 0) stop("Selected sender/receiver/lr_pair not found.")
  ev <- ev[1, , drop = FALSE]
  interpretation <- paste0(
    sender, " -> ", receiver, " ", lr_pair, " is classified as ", ev$communication_range,
    " with local LCS=", signif(ev$lcs_neighborhood, 3),
    " and global potential=", signif(ev$lcs_global, 3), "."
  )
  plt <- tryCatch(plot_lr_evidence(ct_comm, sender, receiver, lr_pair), error = function(e) NULL)
  list(query = list(sender = sender, receiver = receiver, lr_pair = lr_pair), evidence = ev,
       interpretation = interpretation, plot = plt)
}

#' Plot Ligand/Receptor Activity Balance
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param sender Optional sender filter.
#' @param receiver Optional receiver filter.
#' @param active_only Whether to keep active rows.
#' @param top_n Number of strongest rows to display.
#' @param title Optional title.
#' @return A ggplot2 object.
#' @export
plot_lr_activity_balance <- function(ct_comm, sender = NULL, receiver = NULL, active_only = TRUE, top_n = 100, title = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  df <- ct_comm$lr_table
  df <- .filter_celltype_axis(df, sender, receiver)
  if (isTRUE(active_only)) df <- df[df$active %in% TRUE, , drop = FALSE]
  df <- df[is.finite(df$lcs), , drop = FALSE]
  if (nrow(df) == 0) stop(.no_active_events_message(ct_comm, sender, receiver), call. = FALSE)
  df <- utils::head(df[order(df$lcs, decreasing = TRUE), , drop = FALSE], top_n)
  if (is.null(title)) title <- "Ligand and receptor logic activity balance"
  ggplot2::ggplot(df, ggplot2::aes(x = ligand_active_frac_sender, y = receptor_active_frac_receiver)) +
    ggplot2::geom_point(ggplot2::aes(size = lcs, color = communication_range), alpha = 0.75) +
    ggrepel::geom_text_repel(ggplot2::aes(label = lr_pair), size = 3, max.overlaps = 15) +
    ggplot2::scale_color_manual(values = .logiccomm_palettes$range, na.value = "grey70") +
    ggplot2::scale_size_continuous(range = c(2, 8), name = "LCS") +
    ggplot2::labs(title = title, x = "Ligand-active fraction in sender", y = "Receptor-active fraction in receiver") +
    theme_logiccomm()
}

#' Plot Evidence for One Cell-Type L-R Event
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param sender Sender cell type.
#' @param receiver Receiver cell type.
#' @param lr_pair L-R pair.
#' @param title Optional title.
#' @return A ggplot2 object.
#' @export
plot_lr_evidence <- function(ct_comm, sender, receiver, lr_pair, title = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  row <- ct_comm$lr_table[ct_comm$lr_table$sender_type == sender &
                            ct_comm$lr_table$receiver_type == receiver &
                            ct_comm$lr_table$lr_pair == lr_pair, , drop = FALSE]
  if (nrow(row) == 0) stop("Selected sender/receiver/lr_pair not found.")
  row <- row[1, , drop = FALSE]
  support_fraction <- row$n_active_edges / pmax(row$n_edges, 1)
  plot_df <- data.frame(
    metric = c("Primary LCS", "Local LCS", "Global potential", "Ligand sender fraction", "Receptor receiver fraction", "Support fraction"),
    value = c(row$lcs, row$lcs_neighborhood, row$lcs_global, row$ligand_active_frac_sender, row$receptor_active_frac_receiver, support_fraction),
    stringsAsFactors = FALSE
  )
  plot_df$value[!is.finite(plot_df$value)] <- NA_real_
  if (is.null(title)) title <- paste(sender, "->", receiver, lr_pair)
  ggplot2::ggplot(plot_df, ggplot2::aes(x = metric, y = value, fill = metric)) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = title, subtitle = paste("Range:", row$communication_range), x = NULL, y = "Score / fraction") +
    theme_logiccomm() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}

#' Plot Differential Cell-Type Communication Heatmap
#'
#' @param result Output from \code{CompareLogicGroups()}.
#' @param metric Numeric metric.
#' @param top_n Number of features.
#' @param title Optional title.
#' @return A ggplot2 object.
#' @export
plot_differential_celltype_heatmap <- function(result, metric = "asymmetry", top_n = 30, title = NULL) {
  if (inherits(result, "LogicCommDifferential")) {
    df <- result$lr
    if (!metric %in% names(df)) stop("metric not found: ", metric)
    df$sender_receiver <- paste(df$sender_type, "\u2192", df$receiver_type)
    df$lr_label <- df$lr_pair
  } else {
    df <- as.data.frame(result)
    if (!metric %in% names(df)) stop("metric not found: ", metric)
    if (!"lr_pair" %in% names(df)) stop("result must contain lr_pair.")
    parts <- strsplit(as.character(df$lr_pair), "|", fixed = TRUE)
    sender <- vapply(parts, function(x) if (length(x) >= 1) x[1] else "feature", character(1))
    receiver <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
    df$sender_receiver <- ifelse(is.na(receiver), sender, paste(sender, "\u2192", receiver))
    df$lr_label <- vapply(parts, function(x) if (length(x) >= 3) paste(x[3:length(x)], collapse = "|") else x[length(x)], character(1))
  }
  df <- df[order(abs(df[[metric]]), decreasing = TRUE), , drop = FALSE]
  df <- utils::head(df, top_n)
  if (is.null(title)) title <- paste("Differential cell-type communication:", metric)
  ggplot2::ggplot(df, ggplot2::aes(x = sender_receiver, y = lr_label, fill = .data[[metric]])) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_logiccomm_diverging(midpoint = 0, name = metric) +
    ggplot2::labs(title = title, x = "Sender \u2192 receiver", y = "L-R pair") +
    theme_logiccomm() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot Differential Cell-Type Communication Volcano
#'
#' @param result Output from \code{CompareLogicGroups()}.
#' @param x X-axis metric.
#' @param fdr_col FDR column.
#' @param top_n_label Number of labels.
#' @param title Optional title.
#' @return A ggplot2 object.
#' @export
plot_differential_celltype_volcano <- function(result, x = "log2fc_lcs", fdr_col = "fdr_fisher",
                                                fdr_cutoff = 0.05, lfc_threshold = 0,
                                                top_n_label = 15,
                                                title = "Differential cell-type communication",
                                                x_lab = NULL) {
  df <- as.data.frame(result)
  if (!x %in% names(df)) x <- if ("asymmetry" %in% names(df)) "asymmetry" else stop("x metric not found: ", x)
  if (!fdr_col %in% names(df)) fdr_col <- if ("fdr" %in% names(df)) "fdr" else stop("FDR column not found: ", fdr_col)
  df$neglog10 <- -log10(pmax(df[[fdr_col]], .Machine$double.xmin))
  xv <- df[[x]]
  df$significance <- factor(
    ifelse(df[[fdr_col]] <= fdr_cutoff & xv >  lfc_threshold, "Up",
    ifelse(df[[fdr_col]] <= fdr_cutoff & xv < -lfc_threshold, "Down", "n.s.")),
    levels = c("Up", "Down", "n.s."))
  raw_label <- if ("lr_pair" %in% names(df)) as.character(df$lr_pair) else rownames(df)
  df$label <- ifelse(grepl("|", raw_label, fixed = TRUE),
                     .compact_feature_label(raw_label), .short_label(raw_label, 28L))
  sig <- df[df$significance != "n.s.", , drop = FALSE]
  ord <- order(sig[[fdr_col]], -abs(sig[[x]]), na.last = NA)
  lab_df <- sig[utils::head(ord, top_n_label), , drop = FALSE]
  if (is.null(x_lab)) {
    x_lab <- switch(x, log2fc_lcs = expression(log[2]~"fold-change (LCS)"),
                    asymmetry = "Case - Ctrl asymmetry", x)
  }
  subtitle <- sprintf("Coloured: FDR < %.2g%s", fdr_cutoff,
                      if (lfc_threshold > 0) sprintf(" & |effect| > %.2g", lfc_threshold) else "")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = neglog10)) +
    ggplot2::geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed",
                        colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.4)
  if (lfc_threshold > 0) {
    p <- p + ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                                 linetype = "dashed", colour = "grey70", linewidth = 0.35)
  }
  p +
    ggplot2::geom_point(ggplot2::aes(colour = significance), alpha = 0.8, size = 1.9) +
    ggrepel::geom_text_repel(data = lab_df, ggplot2::aes(label = label), size = 2.8,
                             max.overlaps = 16, box.padding = 0.4, min.segment.length = 0,
                             segment.colour = "grey70", segment.size = 0.3,
                             colour = logiccomm_brand$ink, seed = 1) +
    ggplot2::scale_colour_manual(values = .logiccomm_palettes$significance, name = NULL,
                                 drop = FALSE) +
    ggplot2::labs(title = title, subtitle = subtitle, x = x_lab,
                  y = expression(-log[10]~"FDR")) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_logiccomm()
}
