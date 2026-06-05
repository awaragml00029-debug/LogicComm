# R/specificity_plots.R

# Compute 2D node coordinates for a cell-type communication network. Uses an
# igraph force-directed (Fruchterman-Reingold) layout when igraph is available
# (better node separation for many cell types) and otherwise falls back to an
# evenly spaced circle. The user's RNG state is preserved. Returns a
# data.frame(cell_type, x, y) with coordinates normalised to [-1, 1].
.celltype_graph_layout <- function(cell_types, edges, layout = c("auto", "fr", "circle")) {
  layout <- match.arg(layout)
  cell_types <- unique(as.character(cell_types))
  n <- length(cell_types)
  circle <- function() {
    theta <- seq(0, 2 * pi, length.out = n + 1L)[seq_len(n)]
    data.frame(cell_type = cell_types, x = cos(theta), y = sin(theta), stringsAsFactors = FALSE)
  }
  use_fr <- (layout %in% c("auto", "fr")) && requireNamespace("igraph", quietly = TRUE) &&
            n > 2 && !is.null(edges) && nrow(edges) > 0
  if (layout == "circle" || !use_fr) return(circle())
  g <- igraph::graph_from_data_frame(
    data.frame(from = as.character(edges$sender_type),
               to = as.character(edges$receiver_type), stringsAsFactors = FALSE),
    directed = TRUE, vertices = data.frame(name = cell_types, stringsAsFactors = FALSE))
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(1)
  lay <- igraph::layout_with_fr(g)
  if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  norm <- function(v) { r <- range(v); if (diff(r) == 0) return(rep(0, length(v))); 2 * (v - r[1]) / diff(r) - 1 }
  data.frame(cell_type = igraph::V(g)$name, x = norm(lay[, 1]), y = norm(lay[, 2]),
             stringsAsFactors = FALSE)
}

#' Plot Pathway Dominance
#'
#' @param ct_comm Cell-type communication object.
#' @param metric Pathway metric, usually \code{"sum_lcs"}.
#' @param top_n Number of pathways to display.
#' @param title Optional title.
#' @return ggplot2 object.
#' @export
plot_pathway_dominance <- function(ct_comm,
                                   metric = "sum_lcs",
                                   top_n = 15,
                                   title = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  df <- ct_comm$pathway_summary
  if (!metric %in% names(df)) stop("metric not found in pathway_summary: ", metric)
  if (!nrow(df)) stop("No pathway communication rows available.")
  
  agg <- aggregate(df[[metric]], list(pathway = df$pathway), sum, na.rm = TRUE)
  names(agg)[2] <- "value"
  agg <- agg[order(agg$value, decreasing = TRUE), , drop = FALSE]
  agg$fraction <- agg$value / sum(agg$value, na.rm = TRUE)
  agg <- utils::head(agg, top_n)
  agg$pathway <- factor(agg$pathway, levels = rev(agg$pathway))
  
  if (is.null(title)) title <- paste("Pathway dominance:", metric)
  
  ggplot2::ggplot(agg, ggplot2::aes(x = fraction, y = pathway, fill = value)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(), name = "Fraction of total strength") +
    scale_fill_logiccomm_c(name = metric) +
    ggplot2::labs(title = title, y = "Pathway") +
    theme_logiccomm(grid = "x")
}

#' Plot Cell-Type Role Dotplot
#'
#' @param ct_comm Cell-type communication object.
#' @param scaled If TRUE, z-score role scores before plotting.
#' @param title Optional title.
#' @return ggplot2 object.
#' @export
plot_celltype_role_dotplot <- function(ct_comm,
                                       scaled = TRUE,
                                       title = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  df <- ct_comm$role_summary
  role_cols <- c(sender_role_score = "Sender", receiver_role_score = "Receiver",
                 mediator_role_score = "Mediator", influencer_role_score = "Influencer")
  mat <- as.matrix(df[, names(role_cols), drop = FALSE])
  z <- apply(mat, 2, .zscore_safe)
  val <- if (scaled) z else mat
  
  plot_df <- data.frame(
    cell_type = rep(df$cell_type, times = length(role_cols)),
    role = rep(unname(role_cols), each = nrow(df)),
    value = as.numeric(val),
    raw_score = as.numeric(mat),
    evidence = rep(df$communication_evidence_score %||% df$role_confidence %||% 0.5, times = length(role_cols)),
    dominant_role = rep(df$dominant_role, times = length(role_cols)),
    stringsAsFactors = FALSE
  )
  plot_df$role <- factor(plot_df$role, levels = unname(role_cols))
  
  if (is.null(title)) title <- "Cell-type signaling-role profile"
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = role, y = cell_type)) +
    ggplot2::geom_point(ggplot2::aes(size = raw_score, color = value, alpha = evidence)) +
    scale_color_logiccomm_diverging(midpoint = 0, name = if (scaled) "z-score" else "score") +
    ggplot2::scale_alpha_continuous(range = c(0.35, 1), name = "Evidence") +
    ggplot2::scale_size_continuous(range = c(2, 9), name = "Raw score") +
    ggplot2::labs(title = title, x = "Signaling role", y = "Cell type") +
    theme_logiccomm()
}

#' Plot Specificity-Stability Landscape
#'
#' @param sens Output from \code{sensitivity_REO_threshold()}.
#' @param ct_comm Optional cell-type communication object.
#' @param specificity Optional specificity summary.
#' @param top_n_label Number of top points to label.
#' @param label_field Label content: \code{"lr_pair"} (default, short L-R names),
#'   \code{"feature"} (compact sender/receiver/L-R key), or \code{"none"}.
#' @param title Optional title.
#' @return ggplot2 object.
#' @export
plot_specificity_stability <- function(sens,
                                       ct_comm = NULL,
                                       specificity = NULL,
                                       top_n_label = 12,
                                       label_field = c("lr_pair", "feature", "none"),
                                       title = NULL) {
  if (is.null(sens$stability)) stop("sens must be output from sensitivity_REO_threshold().")
  label_field <- match.arg(label_field)
  st <- sens$stability

  # Parse feature names (sender|receiver|lr_pair)
  parts <- .parse_celltype_feature_names(st$feature)
  st$lr_pair <- parts$feature

  if (is.null(specificity)) {
    if (is.null(ct_comm)) stop("Provide ct_comm or specificity.")
    if (is.null(ct_comm$specificity_summary) || !nrow(ct_comm$specificity_summary)) {
      ct_comm <- score_communication_specificity(ct_comm, verbose = FALSE)
    }
    specificity <- ct_comm$specificity_summary
  }

  mi <- match(st$lr_pair, specificity$lr_pair)
  st$pair_specificity <- specificity$pair_specificity[mi]
  st$pathway <- specificity$pathway[mi]
  st$specificity_class <- specificity$specificity_class[mi]
  st$total_lcs <- specificity$total_lcs[mi]
  # Short labels keep the figure readable when cell-type names are long; the
  # full composite key is available via label_field = "feature".
  st$label <- switch(label_field,
                     lr_pair = .short_label(st$lr_pair, 28L),
                     feature = .compact_feature_label(st$feature),
                     none = NA_character_)

  # Label the strongest, most stable, most specific axes; ties (a crowded (1,1)
  # corner) are broken by total LCS so labels do not pile up.
  ord <- order(st$active_fraction, st$pair_specificity, st$total_lcs,
               decreasing = TRUE, na.last = NA)
  st$to_label <- FALSE
  if (label_field != "none" && length(ord) && top_n_label > 0) {
    st$to_label[ord[seq_len(min(top_n_label, length(ord)))]] <- TRUE
  }

  if (is.null(title)) title <- "Communication stability vs specificity"

  p <- ggplot2::ggplot(st, ggplot2::aes(x = active_fraction, y = pair_specificity)) +
    ggplot2::geom_point(ggplot2::aes(size = total_lcs, color = specificity_class), alpha = 0.7) +
    ggplot2::scale_size_continuous(range = c(2, 8), name = "Total LCS") +
    scale_color_logiccomm_d(name = "Specificity class", na.value = "grey70") +
    ggplot2::labs(title = title, x = "REO-threshold stability", y = "Cell-type-pair specificity") +
    theme_logiccomm()
  if (any(st$to_label)) {
    p <- p + ggrepel::geom_text_repel(
      data = st[st$to_label, , drop = FALSE], ggplot2::aes(label = label),
      size = 3, max.overlaps = 12, min.segment.length = 0, box.padding = 0.5,
      segment.color = "grey60", segment.size = 0.3)
  }
  p
}

#' Plot Publication-Style Cell-Type Network
#'
#' @param ct_comm Cell-type communication object.
#' @param metric Pair-level metric.
#' @param top_n_edges Number of strongest non-self edges to show.
#' @param min_value Minimum edge value.
#' @param edge_label "none", "top_pathway", or "n_active_lr".
#' @param layout Node layout: "auto" (force-directed via igraph when available,
#'   otherwise a circle), "fr", or "circle".
#' @param show_self Whether to show self-loops.
#' @param title Optional title.
#' @return ggplot2 object.
#' @export
plot_celltype_network_publication <- function(ct_comm,
                                              metric = "sum_lcs",
                                              top_n_edges = 20,
                                              min_value = 0,
                                              edge_label = c("top_pathway", "none", "n_active_lr"),
                                              layout = c("auto", "fr", "circle"),
                                              show_self = FALSE,
                                              title = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  edge_label <- match.arg(edge_label)
  layout <- match.arg(layout)
  edges <- ct_comm$pair_summary
  if (!metric %in% names(edges)) stop("metric not found in pair_summary: ", metric)
  edges$value <- edges[[metric]]
  edges <- edges[!is.na(edges$value) & edges$value > min_value, , drop = FALSE]
  if (!isTRUE(show_self)) edges <- edges[edges$sender_type != edges$receiver_type, , drop = FALSE]
  edges <- edges[order(edges$value, decreasing = TRUE), , drop = FALSE]
  edges <- utils::head(edges, top_n_edges)

  nodes <- ct_comm$role_summary
  nodes$node_label <- .short_label(as.character(nodes$cell_type), 18L)
  coords <- .celltype_graph_layout(nodes$cell_type, edges, layout = layout)
  nodes$x <- coords$x[match(nodes$cell_type, coords$cell_type)]
  nodes$y <- coords$y[match(nodes$cell_type, coords$cell_type)]

  edge_df <- merge(edges, nodes[, c("cell_type", "x", "y")], by.x = "sender_type", by.y = "cell_type", all.x = TRUE)
  edge_df <- merge(edge_df, nodes[, c("cell_type", "x", "y")], by.x = "receiver_type", by.y = "cell_type", all.x = TRUE, suffixes = c("", "end"))
  edge_df$edge_text <- if (edge_label == "top_pathway") edge_df$top_pathway else if (edge_label == "n_active_lr") as.character(edge_df$n_active_lr) else ""

  if (is.null(title)) title <- paste("Cell-type communication network:", metric)

  p <- ggplot2::ggplot() +
    ggplot2::labs(title = title, subtitle = "Strongest sender -> receiver relationships")

  if (nrow(edge_df)) {
    p <- p + ggplot2::geom_curve(data = edge_df,
                                 ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                                              linewidth = value, colour = top_pathway),
                                 curvature = 0.16,
                                 arrow = grid::arrow(length = grid::unit(0.1, "inches"), type = "closed"),
                                 lineend = "round", alpha = 0.75) +
      ggplot2::scale_linewidth_continuous(range = c(0.4, 3), name = metric) +
      scale_color_logiccomm_d(name = "Top pathway", na.value = "grey75")

    if (edge_label != "none") {
      edge_df$xm <- (edge_df$x + edge_df$xend) / 2
      edge_df$ym <- (edge_df$y + edge_df$yend) / 2
      p <- p + ggrepel::geom_text_repel(data = edge_df, ggplot2::aes(x = xm, y = ym, label = edge_text),
                                        size = 3, max.overlaps = 15, fontface = "italic", colour = "grey30")
    }
  }

  p + ggplot2::geom_point(data = nodes,
                          ggplot2::aes(x = x, y = y, size = hub_score, fill = sender_receiver_balance),
                          shape = 21, colour = "white", stroke = 0.6) +
    ggrepel::geom_text_repel(data = nodes, ggplot2::aes(x = x, y = y, label = node_label),
                             size = 3.4, fontface = "bold", max.overlaps = 30,
                             colour = logiccomm_brand$ink, seed = 1) +
    ggplot2::scale_size_continuous(range = c(4, 13), name = "Hub score") +
    scale_fill_logiccomm_diverging(midpoint = 0, name = "S-R balance") +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.25), colour = logiccomm_brand$ink),
      plot.subtitle = ggplot2::element_text(colour = "grey35", margin = ggplot2::margin(b = 8)),
      legend.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.9)),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}

#' Plot a Pathway-Filtered Cell-Type Network
#'
#' @param ct_comm Cell-type communication object.
#' @param pathway Pathway name to show.
#' @param metric Pathway metric.
#' @param top_n_edges Number of edges to show.
#' @param title Optional title.
#' @return ggplot2 object.
#' @export
plot_pathway_network <- function(ct_comm,
                                 pathway,
                                 metric = "sum_lcs",
                                 top_n_edges = 20,
                                 title = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  df <- ct_comm$pathway_summary
  df <- df[df$pathway %in% pathway, , drop = FALSE]
  if (!nrow(df)) stop("No pathway_summary rows found for pathway: ", paste(pathway, collapse = ", "))
  
  ps <- data.frame(sender_type = df$sender_type, receiver_type = df$receiver_type,
                   n_active_lr = df$n_active_lr, sum_lcs = df$sum_lcs,
                   mean_lcs_active = df$mean_lcs_active,
                   top_pathway = df$pathway, stringsAsFactors = FALSE)
  tmp <- ct_comm
  tmp$pair_summary <- ps
  
  if (is.null(title)) title <- paste("Pathway network:", paste(pathway, collapse = ", "))
  
  plot_celltype_network_publication(tmp, metric = metric, top_n_edges = top_n_edges,
                                    edge_label = "n_active_lr", show_self = FALSE, title = title)
}
