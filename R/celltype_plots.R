# R/celltype_plots.R

#' Plot Cell-Type-Level Communication Network
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param metric Metric for edge weights (\code{"sum_lcs"} or \code{"n_active_lr"}).
#' @param layout Network layout (e.g., \code{"circle"}, \code{"nicely"}).
#' @param min_weight Filter edges below this threshold.
#' @param arrow_size Size of edge arrows.
#' @param show_labels Show cell-type names.
#' @param color_edges_by How to color edges: \code{"range"} or \code{"top_pathway"}.
#' @return A ggplot2 object.
#' @export
plot_celltype_network <- function(ct_comm,
                                  metric = c("sum_lcs", "n_active_lr"),
                                  layout = "circle",
                                  min_weight = 0.05,
                                  arrow_size = 0.2,
                                  show_labels = TRUE,
                                  color_edges_by = c("top_pathway", "range")) {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("igraph required.")
  if (!requireNamespace("ggnetwork", quietly = TRUE)) stop("ggnetwork required.")
  
  metric <- match.arg(metric)
  color_edges_by <- match.arg(color_edges_by)
  pair_sum <- ct_comm$pair_summary
  pair_sum <- pair_sum[pair_sum[[metric]] >= min_weight, , drop = FALSE]
  
  if (nrow(pair_sum) == 0) stop("No edges above min_weight.")
  
  g <- igraph::graph_from_data_frame(pair_sum, directed = TRUE)
  net <- ggnetwork::ggnetwork(g, layout = layout)
  
  # Map role summary stats to nodes
  nodes <- ct_comm$role_summary
  net$hub_score <- nodes$hub_score[match(net$vertex.names, nodes$cell_type)]
  net$balance <- nodes$sender_receiver_balance[match(net$vertex.names, nodes$cell_type)]
  
  # For edge coloring, if range is selected, we need to decide which range to show.
  # Since pair_summary aggregates across L-R pairs, we show the dominant range.
  if (color_edges_by == "range") {
    # Get dominant range per CT pair from lr_table
    df_active <- ct_comm$lr_table[ct_comm$lr_table$active, ]
    range_map <- tapply(df_active$communication_range, paste(df_active$sender_type, df_active$receiver_type, sep = "->"), function(x) names(sort(table(x), decreasing = TRUE))[1])
    net$range <- range_map[paste(net$sender_type, net$receiver_type, sep = "->")]
  }

  p <- ggplot2::ggplot(net, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
    ggnetwork::geom_edges(
      ggplot2::aes(size = !!as.name(metric), 
                   color = if(color_edges_by == "range") range else top_pathway,
                   linetype = if(color_edges_by == "range" && "range" %in% names(net)) range else NULL),
      alpha = 0.6, curvature = 0.15,
      arrow = ggplot2::arrow(length = ggplot2::unit(arrow_size, "cm"), type = "closed")
    ) +
    ggnetwork::geom_nodes(ggplot2::aes(fill = balance, size = hub_score), shape = 21, color = "black") +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                                   midpoint = 0, name = "S-R Balance") +
    ggplot2::scale_size_continuous(range = c(3, 12), name = "Strength") +
    theme_logiccomm()
    
  if (color_edges_by == "range") {
    p <- p + ggplot2::scale_color_manual(values = .logiccomm_palettes$range, name = "Communication Range") +
             ggplot2::scale_linetype_manual(values = c("juxtacrine" = "solid", "paracrine" = "solid", "distal/endocrine" = "dashed"), name = "Range Style")
  } else {
    p <- p + ggplot2::scale_color_viridis_d(option = "turbo", name = "Top Pathway")
  }
    
  if (show_labels) {
    p <- p + ggnetwork::geom_nodetext_repel(ggplot2::aes(label = vertex.names),
                                            size = 3.5, fontface = "bold", box.padding = 0.5)
  }
  p
}

#' Plot Communication Role Landscape
#' @export
plot_celltype_roles <- function(ct_comm) {
  df <- ct_comm$role_summary
  ggplot2::ggplot(df, ggplot2::aes(x = outgoing_strength, y = incoming_strength)) +
    ggplot2::geom_point(ggplot2::aes(size = hub_score, color = sender_receiver_balance)) +
    ggrepel::geom_text_repel(ggplot2::aes(label = cell_type), size = 3) +
    ggplot2::scale_color_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0) +
    ggplot2::scale_size_continuous(range = c(2, 10)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    theme_logiccomm() +
    ggplot2::labs(title = "Signaling Role Landscape", 
                  subtitle = "Outgoing vs Incoming Communication Strength",
                  x = "Outgoing Strength (Source)", y = "Incoming Strength (Target)")
}

#' Plot Logic Communication Range Summary
#' @export
plot_communication_range_summary <- function(ct_comm) {
  df <- ct_comm$lr_table
  df <- df[df$active, , drop = FALSE]
  if (nrow(df) == 0) stop("No active events to plot.")
  
  # Summary by Pathway
  pw_df <- as.data.frame(table(df$pathway, df$communication_range))
  colnames(pw_df) <- c("Pathway", "Range", "Count")
  
  # Sort pathways by total count
  pw_totals <- tapply(pw_df$Count, pw_df$Pathway, sum)
  pw_df$Pathway <- factor(pw_df$Pathway, levels = names(sort(pw_totals)))
  
  ggplot2::ggplot(pw_df, ggplot2::aes(x = Pathway, y = Count, fill = Range)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::scale_fill_manual(values = .logiccomm_palettes$range) +
    ggplot2::coord_flip() +
    theme_logiccomm() +
    ggplot2::labs(title = "Signaling Range by Pathway", 
                  subtitle = "Count of active L-R pairs per range category",
                  y = "Number of active L-R pairs")
}

#' Advanced Bubble Plot of L-R Logic Events
#' @export
plot_lr_bubble_advanced <- function(ct_comm, senders = NULL, receivers = NULL, top_n_per_pair = 5) {
  df <- ct_comm$lr_table
  if (!is.null(senders)) df <- df[df$sender_type %in% senders, ]
  if (!is.null(receivers)) df <- df[df$receiver_type %in% receivers, ]
  df <- df[df$active, ]
  
  # Rank per CT pair
  df <- df |> 
    dplyr::group_by(sender_type, receiver_type) |> 
    dplyr::slice_max(order_by = lcs, n = top_n_per_pair, with_ties = FALSE) |> 
    dplyr::ungroup()
    
  ggplot2::ggplot(df, ggplot2::aes(x = sender_type, y = receiver_type, size = lcs, color = communication_range)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::facet_wrap(~pathway, scales = "free") +
    ggplot2::scale_color_manual(values = .logiccomm_palettes$range) +
    ggplot2::scale_size_continuous(range = c(2, 8)) +
    theme_logiccomm() +
    ggplot2::labs(title = "Cell-Type Logic Communication Hotspots", 
                  subtitle = "Top active L-R pairs facetted by Pathway",
                  x = "Sender Type", y = "Receiver Type") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
