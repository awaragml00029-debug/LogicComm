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

  out_strength <- rowSums(adj_strength, na.rm = TRUE)
  in_strength <- colSums(adj_strength, na.rm = TRUE)
  out_count <- rowSums(adj_count, na.rm = TRUE)
  in_count <- colSums(adj_count, na.rm = TRUE)
  total_strength <- out_strength + in_strength
  hub_score <- total_strength / max(total_strength, 1e-9)
  balance <- (out_strength - in_strength) / pmax(total_strength, 1e-9)

  betweenness <- rep(0, n)
  information <- rep(0, n)
  if (n > 1 && any(adj_strength > 0, na.rm = TRUE)) {
    # adj_strength stores communication strength (sum_lcs): larger = stronger.
    # igraph interprets edge weights as distances/costs (smaller = closer), so
    # mediator betweenness must use the INVERSE strength as a distance. Passing
    # strength directly would route shortest paths through the weakest edges and
    # invert the mediator role. Self-loops (diag) never lie on a shortest path
    # and distort centrality, so they are excluded from the graph.
    if (requireNamespace("igraph", quietly = TRUE)) {
      g <- igraph::graph_from_adjacency_matrix(adj_strength, mode = "directed",
                                               weighted = TRUE, diag = FALSE)
      if (igraph::ecount(g) > 0) {
        w <- as.numeric(igraph::E(g)$weight)
        pos_w <- w[is.finite(w) & w > 0]
        dist_w <- 1 / w
        dist_w[!is.finite(dist_w)] <- if (length(pos_w)) 1e3 / min(pos_w) else 1
        betweenness <- tryCatch(
          as.numeric(igraph::betweenness(g, weights = dist_w, normalized = TRUE)),
          error = function(e) .base_weighted_betweenness(adj_strength))
        # PageRank gives the influencer / information-flow score. It is the robust,
        # well-conditioned analogue of incoming eigenvector centrality (a cell type
        # is influential when influential cell types communicate with it) and,
        # unlike eigen_centrality(directed = TRUE), it does not collapse to zeros or
        # to sink nodes on the acyclic / weakly connected communication graphs that
        # are common here. Higher edge weight = stronger random-walk flow.
        information <- tryCatch(
          as.numeric(igraph::page_rank(g, directed = TRUE, weights = w)$vector),
          error = function(e) .base_weighted_pagerank(adj_strength))
      }
    } else {
      betweenness <- .base_weighted_betweenness(adj_strength)
      information <- .base_weighted_pagerank(adj_strength)
    }
  }
  # igraph::betweenness(normalized = TRUE) returns NaN for n == 1 (division by
  # (n-1)(n-2)); centrality scores can also be non-finite on degenerate graphs.
  betweenness[!is.finite(betweenness)] <- 0
  information[!is.finite(information)] <- 0
  names(betweenness) <- names(information) <- cell_types

  active_lr <- if (!is.null(lr_table) && nrow(lr_table) > 0 && "active" %in% names(lr_table)) {
    lr_table[lr_table$active %in% TRUE, , drop = FALSE]
  } else {
    data.frame()
  }

  count_for <- function(type, side, field = NULL) {
    if (nrow(active_lr) == 0) return(0L)
    if (side == "sender") sub <- active_lr[active_lr$sender_type == type, , drop = FALSE] else sub <- active_lr[active_lr$receiver_type == type, , drop = FALSE]
    if (is.null(field)) return(nrow(sub))
    length(unique(sub[[field]]))
  }

  outgoing_lr_event_count <- vapply(cell_types, count_for, integer(1), side = "sender")
  incoming_lr_event_count <- vapply(cell_types, count_for, integer(1), side = "receiver")
  outgoing_target_type_count <- vapply(cell_types, count_for, integer(1), side = "sender", field = "receiver_type")
  incoming_source_type_count <- vapply(cell_types, count_for, integer(1), side = "receiver", field = "sender_type")
  outgoing_unique_lr_count <- vapply(cell_types, count_for, integer(1), side = "sender", field = "lr_pair")
  incoming_unique_lr_count <- vapply(cell_types, count_for, integer(1), side = "receiver", field = "lr_pair")
  total_event_count <- outgoing_lr_event_count + incoming_lr_event_count

  evidence_score <- total_event_count / max(total_event_count, 1L)
  role_scores <- cbind(Sender = out_strength, Receiver = in_strength, Mediator = betweenness, Influencer = information)
  dominant_idx <- max.col(role_scores, ties.method = "first")
  dominant_role_strict <- colnames(role_scores)[dominant_idx]
  secondary_role <- apply(role_scores, 1, function(x) {
    ord <- order(x, decreasing = TRUE)
    if (length(ord) < 2 || x[ord[2]] <= 0) NA_character_ else names(x)[ord[2]]
  })
  top_score <- apply(role_scores, 1, max, na.rm = TRUE)
  second_score <- apply(role_scores, 1, function(x) {
    ord <- order(x, decreasing = TRUE)
    if (length(ord) < 2) 0 else x[ord[2]]
  })
  separation <- ifelse(top_score > 0, (top_score - second_score) / pmax(top_score, 1e-9), 0)

  low_communication <- hub_score < stats::quantile(hub_score, min_role_hub_quantile, na.rm = TRUE) |
    total_event_count < min_role_event_count
  dominant_role <- dominant_role_strict
  dominant_role[low_communication] <- "Low-communication"

  role_confidence <- 0.5 * evidence_score + 0.5 * separation
  communication_evidence_label <- .label_cut(evidence_score, c(0.33, 0.66), c("Low", "Moderate", "High"))
  role_confidence_label <- .label_cut(role_confidence, c(0.33, 0.66), c("Low", "Moderate", "High"))
  role_separation_label <- .label_cut(separation, c(0.25, 0.6), c("Ambiguous", "Mixed", "Clear"))
  role_reliability_label <- ifelse(low_communication, "Low-communication", role_confidence_label)
  role_biological_interpretation <- ifelse(
    dominant_role == "Low-communication", "Insufficient active communication evidence for confident role assignment.",
    ifelse(dominant_role_strict == "Sender", "Candidate source cell type with stronger outgoing than incoming communication.",
           ifelse(dominant_role_strict == "Receiver", "Candidate target cell type with stronger incoming communication.",
                  ifelse(dominant_role_strict == "Mediator", "Candidate bridge cell type connecting communication modules.",
                         "Candidate influencer cell type with broad network centrality.")))
  )

  data.frame(
    cell_type = cell_types,
    sender_role_score = as.numeric(out_strength),
    receiver_role_score = as.numeric(in_strength),
    mediator_role_score = as.numeric(betweenness),
    influencer_role_score = as.numeric(information),
    outdegree_score = as.numeric(out_strength),
    indegree_score = as.numeric(in_strength),
    betweenness_score = as.numeric(betweenness),
    information_score = as.numeric(information),
    hub_score = as.numeric(hub_score),
    sender_receiver_balance = as.numeric(balance),
    outgoing_lr_event_count = as.integer(outgoing_lr_event_count),
    incoming_lr_event_count = as.integer(incoming_lr_event_count),
    outgoing_target_type_count = as.integer(outgoing_target_type_count),
    incoming_source_type_count = as.integer(incoming_source_type_count),
    outgoing_unique_lr_count = as.integer(outgoing_unique_lr_count),
    incoming_unique_lr_count = as.integer(incoming_unique_lr_count),
    communication_evidence_score = as.numeric(evidence_score),
    communication_evidence_label = communication_evidence_label,
    role_confidence = as.numeric(role_confidence),
    role_confidence_label = role_confidence_label,
    role_separation_label = role_separation_label,
    role_reliability_label = role_reliability_label,
    dominant_role_strict = dominant_role_strict,
    secondary_role = secondary_role,
    dominant_role = dominant_role,
    low_communication = as.logical(low_communication),
    role_biological_interpretation = role_biological_interpretation,
    stringsAsFactors = FALSE
  )
}

.base_weighted_betweenness <- function(adj_strength) {
  w <- as.matrix(adj_strength)
  n <- nrow(w)
  if (n < 3) return(rep(0, n))
  w[!is.finite(w) | w <= 0] <- 0
  diag(w) <- 0
  cb <- numeric(n)
  tol <- 1e-12

  for (s in seq_len(n)) {
    pred <- vector("list", n)
    sigma <- numeric(n)
    sigma[s] <- 1
    dist <- rep(Inf, n)
    dist[s] <- 0
    seen <- rep(FALSE, n)
    stack <- integer(0)

    repeat {
      candidates <- which(!seen & is.finite(dist))
      if (!length(candidates)) break
      v <- candidates[which.min(dist[candidates])]
      seen[v] <- TRUE
      stack <- c(stack, v)
      for (to in which(w[v, ] > 0)) {
        alt <- dist[v] + 1 / w[v, to]
        if (alt < dist[to] - tol) {
          dist[to] <- alt
          sigma[to] <- sigma[v]
          pred[[to]] <- v
        } else if (abs(alt - dist[to]) <= tol) {
          sigma[to] <- sigma[to] + sigma[v]
          pred[[to]] <- unique(c(pred[[to]], v))
        }
      }
    }

    delta <- numeric(n)
    for (node in rev(stack)) {
      if (length(pred[[node]]) && sigma[node] > 0) {
        for (p in pred[[node]]) {
          delta[p] <- delta[p] + (sigma[p] / sigma[node]) * (1 + delta[node])
        }
      }
      if (node != s) cb[node] <- cb[node] + delta[node]
    }
  }

  cb / ((n - 1) * (n - 2))
}

.base_weighted_pagerank <- function(adj_strength, damping = 0.85,
                                    max_iter = 200, tol = 1e-10) {
  w <- as.matrix(adj_strength)
  n <- nrow(w)
  if (n == 0) return(numeric(0))
  w[!is.finite(w) | w <= 0] <- 0
  diag(w) <- 0
  row_total <- rowSums(w)
  transition <- matrix(1 / n, nrow = n, ncol = n)
  non_dangling <- row_total > 0
  transition[non_dangling, ] <- w[non_dangling, , drop = FALSE] / row_total[non_dangling]
  pr <- rep(1 / n, n)

  for (i in seq_len(max_iter)) {
    next_pr <- (1 - damping) / n + damping * as.numeric(crossprod(transition, pr))
    if (sum(abs(next_pr - pr)) < tol) {
      pr <- next_pr
      break
    }
    pr <- next_pr
  }
  pr
}

#' Plot Sender-Receiver Role Positioning
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param label Logical; label cell types. Default: \code{TRUE}.
#' @param title Optional plot title.
#' @param subtitle Optional plot subtitle.
#' @return A ggplot2 object.
#' @export
plot_celltype_roles <- function(ct_comm, label = TRUE,
                                title = "Cell-type communication role positioning",
                                subtitle = "Above the diagonal = receiver-leaning; below = sender-leaning") {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  df <- ct_comm$role_summary
  df$node_label <- .short_label(as.character(df$cell_type), 20L)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = sender_role_score, y = receiver_role_score)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey65") +
    ggplot2::geom_point(ggplot2::aes(size = hub_score, fill = dominant_role),
                        shape = 21, colour = "white", stroke = 0.4, alpha = 0.95) +
    ggplot2::scale_fill_manual(values = .logiccomm_palettes$roles, na.value = "grey70",
                               name = "Dominant role") +
    ggplot2::scale_size_continuous(range = c(2.5, 9), name = "Hub score") +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = "Outgoing strength score", y = "Incoming strength score") +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 4))) +
    theme_logiccomm()

  if (isTRUE(label)) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = node_label), size = 3.2, max.overlaps = 20,
                                      colour = logiccomm_brand$ink, box.padding = 0.4,
                                      segment.colour = "grey70", segment.size = 0.3, seed = 1)
  }
  p
}

#' Plot Cell-Type Role Scores as a Heatmap
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
    scale_fill_logiccomm_diverging(midpoint = 0,
                                   name = if (scaled) "z-score" else "score") +
    ggplot2::labs(title = title, x = "Role score", y = "Cell type") +
    theme_logiccomm() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}

#' Plot Cell-Type Role Scores in Polar Coordinates
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param cell_types Optional cell types to display.
#' @param scaled Whether to z-score role scores.
#' @param title Optional title.
#' @return A ggplot2 object.
#' @export
plot_celltype_role_radar <- function(ct_comm, cell_types = NULL, scaled = TRUE, title = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  df <- ct_comm$role_summary
  if (!is.null(cell_types)) df <- df[df$cell_type %in% cell_types, , drop = FALSE]
  if (nrow(df) == 0) stop("No role rows to plot.")
  role_cols <- c("sender_role_score", "receiver_role_score", "mediator_role_score", "influencer_role_score")
  mat <- as.matrix(df[, role_cols, drop = FALSE])
  val <- if (isTRUE(scaled)) apply(mat, 2, .rescale01) else mat
  plot_df <- data.frame(
    cell_type = rep(df$cell_type, times = length(role_cols)),
    role = rep(c("Sender", "Receiver", "Mediator", "Influencer"), each = nrow(df)),
    value = as.numeric(val),
    stringsAsFactors = FALSE
  )
  if (is.null(title)) title <- "Cell-type role radar profile"
  plot_df$role <- factor(plot_df$role, levels = c("Sender", "Receiver", "Mediator", "Influencer"))
  ggplot2::ggplot(plot_df, ggplot2::aes(x = role, y = value, fill = role)) +
    ggplot2::geom_col(width = 1, alpha = 0.9, colour = "white", linewidth = 0.25) +
    ggplot2::coord_polar() +
    ggplot2::facet_wrap(~cell_type) +
    ggplot2::scale_fill_manual(values = .logiccomm_palettes$roles, name = "Role") +
    ggplot2::labs(title = title, x = NULL, y = if (scaled) "Rescaled score" else "Score") +
    theme_logiccomm() +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(size = ggplot2::rel(0.7)))
}

.zscore_safe <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

.rescale01 <- function(x) {
  x <- as.numeric(x)
  if (!any(is.finite(x))) return(rep(0, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || rng[1] == rng[2]) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

.label_cut <- function(x, cuts, labels) {
  out <- ifelse(x < cuts[1], labels[1], ifelse(x < cuts[2], labels[2], labels[3]))
  out[!is.finite(x)] <- labels[1]
  as.character(out)
}
