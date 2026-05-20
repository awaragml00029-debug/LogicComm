# R/celltype_roles.R

#' @keywords internal
.communication_role_summary <- function(adj_strength, adj_count,
                                        active_edge_mat = NULL,
                                        lr_table = NULL,
                                        min_role_hub_quantile = 0.2,
                                        min_role_event_count = 5) {
  cell_types <- rownames(adj_strength)
  n <- length(cell_types)
  if (n == 0) return(data.frame())

  # Centrality-inspired scores
  out_strength <- rowSums(adj_strength, na.rm = TRUE)
  in_strength  <- colSums(adj_strength, na.rm = TRUE)
  out_count    <- rowSums(adj_count, na.rm = TRUE)
  in_count     <- colSums(adj_count, na.rm = TRUE)

  # Normalization
  hub_score <- (out_strength + in_strength)
  hub_score <- hub_score / max(hub_score, 1e-9)

  balance <- (out_strength - in_strength) / pmax(out_strength + in_strength, 1e-9)

  # Role labels
  dominant_role <- ifelse(out_strength > in_strength, "Sender", "Receiver")
  
  # More complex centrality metrics usually require igraph
  # Placeholder for more complex metrics if igraph is available
  betweenness <- rep(0, n)
  information <- rep(0, n)
  
  if (requireNamespace("igraph", quietly = TRUE)) {
    g <- igraph::graph_from_adjacency_matrix(adj_strength, mode = "directed", weighted = TRUE)
    betweenness <- igraph::betweenness(g, normalized = TRUE)
    information <- tryCatch(igraph::eigen_centrality(g)$vector, error = function(e) rep(0, n))
  }

  res <- data.frame(
    cell_type = cell_types,
    sender_role_score = as.numeric(out_strength),
    receiver_role_score = as.numeric(in_strength),
    mediator_role_score = as.numeric(betweenness),
    influencer_role_score = as.numeric(information),
    hub_score = as.numeric(hub_score),
    sender_receiver_balance = as.numeric(balance),
    dominant_role = dominant_role,
    stringsAsFactors = FALSE
  )
  
  # Add low-communication flag
  res$low_communication <- res$hub_score < stats::quantile(res$hub_score, min_role_hub_quantile) |
    (out_count + in_count) < min_role_event_count
  
  res$dominant_role[res$low_communication] <- "Low-communication"
  
  res
}

#' Plot Sender-Receiver Role Positioning
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param label Logical; label cell types. Default: \code{TRUE}.
#' @param title Optional plot title.
#' @return A ggplot2 object.
#' @export
plot_celltype_roles <- function(ct_comm, label = TRUE, title = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  df <- ct_comm$role_summary
  if (is.null(title)) title <- "Cell-type communication role positioning"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = sender_role_score, y = receiver_role_score)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey65") +
    ggplot2::geom_point(ggplot2::aes(size = hub_score, color = dominant_role), alpha = 0.85) +
    ggplot2::labs(title = title,
                  x = "Outgoing strength score",
                  y = "Incoming strength score",
                  color = "Dominant role",
                  size = "Hub score") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  if (isTRUE(label)) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = cell_type), size = 3.4,
                                      max.overlaps = Inf)
  }
  p
}

#' Plot Cell-Type Role Scores as a Heatmap
#'
#' Displays the four centrality-inspired role scores (sender, receiver, mediator,
#' influencer) for each cell type. Values can be z-scored by role to emphasize
#' relative positioning.
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param role_cols Role-score columns to display.
#' @param scaled If TRUE, z-score each role column before plotting.
#' @param title Optional title.
#' @return A ggplot2 object.
#' @export
plot_celltype_role_heatmap <- function(ct_comm,
                                       role_cols = c("sender_role_score", "receiver_role_score",
                                                     "mediator_role_score", "influencer_role_score"),
                                       scaled = TRUE,
                                       title = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  df <- ct_comm$role_summary
  missing <- setdiff(role_cols, names(df))
  if (length(missing) > 0) stop("Missing role columns: ", paste(missing, collapse = ", "))
  mat <- as.matrix(df[, role_cols, drop = FALSE])
  if (isTRUE(scaled)) mat <- apply(mat, 2, .zscore_safe)
  plot_df <- data.frame(
    cell_type = rep(df$cell_type, times = length(role_cols)),
    role = rep(role_cols, each = nrow(df)),
    value = as.numeric(mat),
    stringsAsFactors = FALSE
  )
  plot_df$role <- factor(plot_df$role, levels = role_cols)
  if (is.null(title)) title <- if (scaled) "Cell-type role scores (z-scored)" else "Cell-type role scores"
  ggplot2::ggplot(plot_df, ggplot2::aes(x = role, y = cell_type, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                                  midpoint = 0, name = if (scaled) "z-score" else "score") +
    ggplot2::labs(title = title, x = "Role score", y = "Cell type") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}
