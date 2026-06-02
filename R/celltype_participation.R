# R/celltype_participation.R — Cell-type communication participation profiles

#' Summarize Per-Cell-Type Communication Participation
#'
#' Quantifies, for each cell type, how broadly its cells participate in
#' communication and what that communication is composed of. Participation is
#' defined on the within-cell REO activity of ligand/receptor genes (it does not
#' require a neighborhood graph): a cell is a candidate sender when at least one
#' ligand complex is active, a candidate receiver when at least one receptor
#' complex is active, and "communicating" when either holds. Pathway and
#' ligand-receptor composition are derived from the active cell-type-resolved
#' events already scored in \code{ct_comm}.
#'
#' @details
#' This complements \code{\link{summarize_celltype_communication}} by answering
#' practical questions about a candidate signaling subgroup:
#' \itemize{
#'   \item What fraction of the cells in a cell type actually express the
#'     ligand/receptor machinery (\code{frac_communicating})?
#'   \item Which signaling pathways and which ligand-receptor pairs make up that
#'     cell type's outgoing and incoming communication
#'     (\code{pathway_composition}, \code{lr_composition})?
#'   \item Is the cell type a major communication hub, combining communication
#'     strength, participation breadth, and active-event support
#'     (\code{importance_label})?
#' }
#' The importance label is a transparent prioritization heuristic for hypothesis
#' generation, not a statistical test.
#'
#' @param ct_comm A \code{LogicCommCellTypeComm} object from
#'   \code{\link{summarize_celltype_communication}}.
#' @param reo_mat The binary REO matrix (or \code{LogicCommREOResult}) used to
#'   build \code{ct_comm}. Needed for per-cell ligand/receptor activity.
#' @param min_active_events Minimum active sender/receiver L-R events required
#'   before a cell type can be labelled a major hub. Default: \code{5}.
#' @param verbose Print progress messages. Default: \code{TRUE}.
#'
#' @return A list of class \code{LogicCommParticipation} with \code{cell_participation},
#'   \code{group_participation}, \code{pathway_composition}, and \code{lr_composition}.
#' @seealso \code{\link{plot_celltype_participation}},
#'   \code{\link{plot_celltype_pathway_composition}},
#'   \code{\link{plot_celltype_communication_profile}}
#' @export
summarize_celltype_participation <- function(ct_comm,
                                             reo_mat,
                                             min_active_events = 5,
                                             verbose = TRUE) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  if (is.list(reo_mat) && !is.null(reo_mat$logic)) reo_mat <- reo_mat$logic
  labels <- ct_comm$cell_labels
  if (is.null(labels) || !length(labels)) stop("ct_comm$cell_labels is empty.")
  lr_db <- ct_comm$lr_db
  cells <- names(labels)
  if (is.null(cells)) stop("ct_comm$cell_labels must be named by cell.")
  if (!all(cells %in% colnames(reo_mat))) {
    stop("reo_mat must contain all cells in ct_comm$cell_labels; recompute on the same cells.")
  }
  reo_sub <- reo_mat[, cells, drop = FALSE]
  labels <- as.character(labels)

  if (isTRUE(verbose)) message("[Participation] Scoring per-cell ligand/receptor activity...")
  act <- score_lr_activity(reo_sub, lr_db = lr_db, mode = "both", aggregate = TRUE, verbose = FALSE)
  sender_n <- as.numeric(Matrix::colSums(act$sender_mat))
  receiver_n <- as.numeric(Matrix::colSums(act$receiver_mat))
  sender_active <- sender_n > 0
  receiver_active <- receiver_n > 0
  communicating <- sender_active | receiver_active
  comm_score <- as.numeric(act$comm_score[cells])

  cell_participation <- data.frame(
    cell = cells,
    cell_type = labels,
    n_active_ligand_pairs = as.integer(sender_n),
    n_active_receptor_pairs = as.integer(receiver_n),
    sender_active = sender_active,
    receiver_active = receiver_active,
    communicating = communicating,
    comm_score = comm_score,
    stringsAsFactors = FALSE
  )

  cell_types <- sort(unique(labels))
  group_participation <- do.call(rbind, lapply(cell_types, function(g) {
    idx <- labels == g
    data.frame(
      cell_type = g,
      n_cells = sum(idx),
      n_communicating = sum(communicating[idx]),
      frac_communicating = mean(communicating[idx]),
      n_sender_active = sum(sender_active[idx]),
      frac_sender_active = mean(sender_active[idx]),
      n_receiver_active = sum(receiver_active[idx]),
      frac_receiver_active = mean(receiver_active[idx]),
      mean_comm_score = mean(comm_score[idx]),
      stringsAsFactors = FALSE
    )
  }))

  # Join communication-strength evidence from the role summary.
  roles <- ct_comm$role_summary
  role_cols <- c("hub_score", "sender_role_score", "receiver_role_score",
                 "outgoing_lr_event_count", "incoming_lr_event_count")
  for (col in role_cols) {
    group_participation[[col]] <- if (!is.null(roles) && col %in% names(roles)) {
      roles[[col]][match(group_participation$cell_type, roles$cell_type)]
    } else {
      rep(NA_real_, nrow(group_participation))
    }
  }
  hub <- ifelse(is.finite(group_participation$hub_score), group_participation$hub_score, 0)
  total_events <- rowSums(cbind(
    ifelse(is.finite(group_participation$outgoing_lr_event_count), group_participation$outgoing_lr_event_count, 0),
    ifelse(is.finite(group_participation$incoming_lr_event_count), group_participation$incoming_lr_event_count, 0)
  ))
  group_participation$total_active_lr_events <- as.integer(total_events)
  importance <- 0.5 * hub + 0.3 * group_participation$frac_communicating + 0.2 * .rescale01(total_events)
  group_participation$importance_score <- importance
  group_participation$importance_label <- ifelse(
    importance >= 0.66 & total_events >= min_active_events, "Major hub",
    ifelse(importance >= 0.33 & total_events >= min_active_events, "Moderate", "Minor / peripheral"))
  group_participation$is_major_hub <- group_participation$importance_label == "Major hub"
  group_participation <- group_participation[order(-group_participation$importance_score), , drop = FALSE]
  rownames(group_participation) <- NULL

  pathway_composition <- .composition_table(ct_comm$lr_table, key = "pathway")
  lr_composition <- .composition_table(ct_comm$lr_table, key = "lr_pair", extra = c("pathway", "ligand", "receptor"))

  res <- list(
    cell_participation = cell_participation,
    group_participation = group_participation,
    pathway_composition = pathway_composition,
    lr_composition = lr_composition,
    params = list(min_active_events = min_active_events,
                  n_cells = length(labels),
                  n_cell_types = length(cell_types))
  )
  class(res) <- "LogicCommParticipation"
  if (isTRUE(verbose)) {
    n_hub <- sum(group_participation$is_major_hub, na.rm = TRUE)
    message(sprintf("[Participation] Done. %d / %d cell types flagged as major communication hubs.",
                    n_hub, nrow(group_participation)))
  }
  res
}

#' Composition of active communication by pathway or L-R pair and direction
#' @keywords internal
.composition_table <- function(lr_table, key = "pathway", extra = character(0)) {
  empty <- data.frame(cell_type = character(), direction = character(),
                      stringsAsFactors = FALSE)
  if (is.null(lr_table) || !nrow(lr_table)) return(empty)
  act <- lr_table[lr_table$active %in% TRUE & is.finite(lr_table$lcs) & lr_table$lcs > 0, , drop = FALSE]
  if (!nrow(act)) return(empty)

  one_direction <- function(df, type_col, direction) {
    if (!nrow(df)) return(NULL)
    grp <- split(seq_len(nrow(df)), list(df[[type_col]], df[[key]]), drop = TRUE)
    rows <- lapply(grp, function(ii) {
      sub <- df[ii, , drop = FALSE]
      out <- data.frame(
        cell_type = sub[[type_col]][1],
        direction = direction,
        key_value = as.character(sub[[key]][1]),
        n_active_lr = nrow(sub),
        sum_lcs = sum(sub$lcs, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      for (e in intersect(extra, names(sub))) out[[e]] <- as.character(sub[[e]][1])
      out
    })
    do.call(rbind, rows)
  }

  out_tab <- one_direction(act, "sender_type", "outgoing")
  in_tab <- one_direction(act, "receiver_type", "incoming")
  comb <- do.call(rbind, list(out_tab, in_tab))
  if (is.null(comb) || !nrow(comb)) return(empty)

  # Share within each cell_type x direction (sums to 1).
  comb$share <- stats::ave(comb$sum_lcs, paste(comb$cell_type, comb$direction), FUN = function(x) {
    s <- sum(x, na.rm = TRUE); if (s > 0) x / s else rep(0, length(x))
  })
  names(comb)[names(comb) == "key_value"] <- key
  comb <- comb[order(comb$cell_type, comb$direction, -comb$share), , drop = FALSE]
  rownames(comb) <- NULL
  comb
}

#' @export
print.LogicCommParticipation <- function(x, ...) {
  gp <- x$group_participation
  cat(sprintf("LogicCommParticipation | %d cell types | %d major hub(s)\n",
              nrow(gp), sum(gp$is_major_hub, na.rm = TRUE)))
  show <- intersect(c("cell_type", "n_cells", "frac_communicating",
                      "frac_sender_active", "frac_receiver_active",
                      "total_active_lr_events", "importance_label"), names(gp))
  print(utils::head(gp[, show, drop = FALSE], 10), row.names = FALSE)
  invisible(x)
}

#' Plot Per-Cell-Type Communication Participation
#'
#' Shows the fraction of cells in each cell type that express the
#' ligand/receptor machinery for communication (as sender, receiver, or either).
#'
#' @param participation A \code{LogicCommParticipation} object.
#' @param mode Which fractions to show: \code{"both"} (sender and receiver),
#'   \code{"communicating"} (combined only), \code{"sender"}, or \code{"receiver"}.
#' @param top_n Optional number of top cell types (by importance) to show.
#' @param title Optional plot title.
#' @return A ggplot2 object.
#' @export
plot_celltype_participation <- function(participation,
                                        mode = c("both", "communicating", "sender", "receiver"),
                                        top_n = NULL,
                                        title = NULL) {
  stopifnot(inherits(participation, "LogicCommParticipation"))
  mode <- match.arg(mode)
  gp <- participation$group_participation
  if (!is.null(top_n)) gp <- utils::head(gp, top_n)
  cols <- switch(mode,
    both = c(frac_sender_active = "Sender-active", frac_receiver_active = "Receiver-active"),
    communicating = c(frac_communicating = "Communicating"),
    sender = c(frac_sender_active = "Sender-active"),
    receiver = c(frac_receiver_active = "Receiver-active"))
  long <- do.call(rbind, lapply(names(cols), function(cl) {
    data.frame(cell_type = gp$cell_type, participation_type = cols[[cl]],
               fraction = gp[[cl]], stringsAsFactors = FALSE)
  }))
  long$cell_type <- factor(long$cell_type, levels = rev(gp$cell_type))
  if (is.null(title)) title <- "Fraction of cells participating in communication"
  ggplot2::ggplot(long, ggplot2::aes(x = fraction, y = cell_type, fill = participation_type)) +
    ggplot2::geom_col(position = if (length(cols) > 1) ggplot2::position_dodge(width = 0.7) else "stack",
                      width = 0.7) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1), name = "Fraction of cells") +
    ggplot2::scale_fill_manual(values = c("Sender-active" = "#d6604d", "Receiver-active" = "#4393c3",
                                          "Communicating" = "#5aae61"), name = NULL) +
    ggplot2::labs(title = title, y = NULL) +
    theme_logiccomm()
}

#' Plot Pathway (or L-R) Composition of Cell-Type Communication
#'
#' Stacked composition of each cell type's active communication, showing what
#' share of its outgoing or incoming communication strength runs through each
#' signaling pathway (or ligand-receptor pair).
#'
#' @param participation A \code{LogicCommParticipation} object.
#' @param direction \code{"outgoing"} (cell type as sender) or \code{"incoming"}
#'   (cell type as receiver).
#' @param level \code{"pathway"} or \code{"lr_pair"}.
#' @param top_n Number of top pathways/pairs to show; the rest are grouped as
#'   \code{"Other"}. Default: \code{8}.
#' @param title Optional plot title.
#' @return A ggplot2 object.
#' @export
plot_celltype_pathway_composition <- function(participation,
                                              direction = c("outgoing", "incoming"),
                                              level = c("pathway", "lr_pair"),
                                              top_n = 8,
                                              title = NULL) {
  stopifnot(inherits(participation, "LogicCommParticipation"))
  direction <- match.arg(direction)
  level <- match.arg(level)
  tab <- if (level == "pathway") participation$pathway_composition else participation$lr_composition
  if (is.null(tab) || !nrow(tab)) stop("No active communication composition to plot.")
  tab <- tab[tab$direction == direction, , drop = FALSE]
  if (!nrow(tab)) stop("No ", direction, " communication composition to plot.")
  tab$grp_key <- tab[[level]]

  # Keep top_n keys per cell type, collapse the rest into "Other".
  tab <- do.call(rbind, lapply(split(tab, tab$cell_type), function(d) {
    d <- d[order(-d$share), , drop = FALSE]
    if (nrow(d) > top_n) {
      other <- data.frame(cell_type = d$cell_type[1], direction = direction,
                          grp_key = "Other", share = sum(d$share[(top_n + 1):nrow(d)]),
                          stringsAsFactors = FALSE)
      d <- rbind(d[seq_len(top_n), c("cell_type", "direction", "grp_key", "share")], other)
    } else {
      d <- d[, c("cell_type", "direction", "grp_key", "share")]
    }
    d
  }))
  keys <- unique(tab$grp_key[tab$grp_key != "Other"])
  tab$grp_key <- factor(tab$grp_key, levels = c(keys, "Other"))
  if (is.null(title)) {
    dir_label <- paste0(toupper(substring(direction, 1, 1)), substring(direction, 2))
    title <- sprintf("%s communication composition by %s",
                     dir_label, if (level == "pathway") "pathway" else "L-R pair")
  }
  ggplot2::ggplot(tab, ggplot2::aes(x = share, y = cell_type, fill = grp_key)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(), name = "Share of communication strength") +
    ggplot2::labs(title = title, y = NULL, fill = if (level == "pathway") "Pathway" else "L-R pair") +
    theme_logiccomm()
}

#' Plot a Single Cell Type's Communication Profile
#'
#' A one-cell-type "profile card": the top ligand-receptor pairs that make up the
#' cell type's outgoing (as sender) and incoming (as receiver) communication,
#' coloured by pathway, with participation and hub-importance reported in the
#' subtitle. Answers, for a chosen subgroup, which pathways and which genes its
#' communication runs through.
#'
#' @param participation A \code{LogicCommParticipation} object.
#' @param cell_type The cell type to profile.
#' @param top_n Number of top L-R pairs per direction. Default: \code{10}.
#' @return A ggplot2 object.
#' @export
plot_celltype_communication_profile <- function(participation, cell_type, top_n = 10) {
  stopifnot(inherits(participation, "LogicCommParticipation"))
  gp <- participation$group_participation
  if (!cell_type %in% gp$cell_type) stop("cell_type '", cell_type, "' not found.")
  lc <- participation$lr_composition
  lc <- lc[lc$cell_type == cell_type, , drop = FALSE]
  if (!nrow(lc)) stop("No active communication for cell type '", cell_type, "'.")
  lc$direction <- factor(ifelse(lc$direction == "outgoing", "Sending (ligand)", "Receiving (receptor)"),
                         levels = c("Sending (ligand)", "Receiving (receptor)"))
  lc <- do.call(rbind, lapply(split(lc, lc$direction, drop = TRUE), function(d) {
    utils::head(d[order(-d$share), , drop = FALSE], top_n)
  }))
  lc$lr_pair <- stats::reorder(lc$lr_pair, lc$share)
  row <- gp[gp$cell_type == cell_type, ]
  subtitle <- sprintf(
    "%s | %.0f%% of cells communicating (%.0f%% sender, %.0f%% receiver) | %d active L-R events | %s",
    cell_type, 100 * row$frac_communicating, 100 * row$frac_sender_active,
    100 * row$frac_receiver_active, row$total_active_lr_events, row$importance_label)
  ggplot2::ggplot(lc, ggplot2::aes(x = share, y = lr_pair, fill = pathway)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_wrap(~ direction, scales = "free_y") +
    ggplot2::scale_x_continuous(labels = scales::percent_format(), name = "Share of communication strength") +
    ggplot2::labs(title = paste("Communication profile:", cell_type), subtitle = subtitle,
                  y = NULL, fill = "Pathway") +
    theme_logiccomm()
}
