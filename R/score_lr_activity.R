# R/score_lr_activity.R

#' Score Per-Cell Ligand-Receptor Activity
#'
#' Assigns each cell a binary activity score for each L-R pair in the database,
#' based on the REO logic matrix. This is useful for dimensionality reduction,
#' trajectory analysis, or cell ranking by communication potential.
#'
#' @details
#' Two activity types are scored per cell per L-R pair:
#' \itemize{
#'   \item \strong{Sender score}: ligand complex is active (all subunits above anchor)
#'   \item \strong{Receiver score}: receptor complex is active (all subunits above anchor)
#' }
#' An aggregate \strong{communication score} is computed per cell as the fraction
#' of L-R pairs for which the cell is either a sender or a receiver.
#'
#' @param reo_mat REO binary matrix (genes x cells) from \code{\link{calc_REO_matrix}}.
#' @param lr_db LR database. Defaults to \code{lr_pairs_human}.
#' @param mode Character. \code{"sender"}, \code{"receiver"}, or \code{"both"}.
#'   Default: \code{"both"}.
#' @param aggregate Logical. If \code{TRUE}, also returns a per-cell aggregate
#'   communication score. Default: \code{TRUE}.
#' @param verbose Logical. Default: \code{TRUE}.
#'
#' @return A list with:
#'   \describe{
#'     \item{sender_mat}{Logical sparse matrix (lr_pairs x cells) for sender activity}
#'     \item{receiver_mat}{Logical sparse matrix (lr_pairs x cells) for receiver activity}
#'     \item{sender_score}{Numeric vector (cells): fraction of L-R pairs as sender}
#'     \item{receiver_score}{Numeric vector (cells): fraction of L-R pairs as receiver}
#'     \item{comm_score}{Numeric vector (cells): combined communication potential}
#'   }
#'
#' @examples
#' \dontrun{
#' data(lr_pairs_human)
#' reo <- calc_REO_matrix(my_seurat, lr_genes = all_lr_genes(lr_pairs_human))
#' scores <- score_lr_activity(reo)
#'
#' # Add to Seurat metadata
#' my_seurat$comm_score <- scores$comm_score
#' }
#'
#' @importFrom Matrix Matrix
#' @export
score_lr_activity <- function(reo_mat,
                               lr_db     = lr_pairs_human,
                               mode      = "both",
                               aggregate = TRUE,
                               verbose   = TRUE) {

  mode <- match.arg(mode, c("sender", "receiver", "both"))

  n_pairs <- nrow(lr_db)
  n_cells <- ncol(reo_mat)
  cell_names <- colnames(reo_mat)

  if (verbose) message(sprintf(
    "[Score] Scoring %d L-R pairs for %d cells...", n_pairs, n_cells))

  sender_list   <- vector("list", n_pairs)
  receiver_list <- vector("list", n_pairs)

  for (i in seq_len(n_pairs)) {
    lig_genes <- lr_db$ligand_genes[[i]]
    rec_genes <- lr_db$receptor_genes[[i]]

    if (mode %in% c("sender","both")) {
      sender_list[[i]] <- as.integer(
        .resolve_complex_logic(lig_genes, reo_mat))
    }
    if (mode %in% c("receiver","both")) {
      receiver_list[[i]] <- as.integer(
        .resolve_complex_logic(rec_genes, reo_mat))
    }
  }

  result <- list()

  if (mode %in% c("sender","both")) {
    sender_mat <- do.call(rbind, sender_list)
    rownames(sender_mat) <- lr_db$lr_pair
    colnames(sender_mat) <- cell_names
    result$sender_mat   <- Matrix::Matrix(sender_mat, sparse = TRUE)
    result$sender_score <- setNames(
      colMeans(sender_mat, na.rm = TRUE), cell_names)
  }

  if (mode %in% c("receiver","both")) {
    receiver_mat <- do.call(rbind, receiver_list)
    rownames(receiver_mat) <- lr_db$lr_pair
    colnames(receiver_mat) <- cell_names
    result$receiver_mat   <- Matrix::Matrix(receiver_mat, sparse = TRUE)
    result$receiver_score <- setNames(
      colMeans(receiver_mat, na.rm = TRUE), cell_names)
  }

  if (aggregate && mode == "both") {
    # Combined: fraction of LR pairs where cell is active as either sender or receiver
    combined <- (sender_mat > 0) | (receiver_mat > 0)
    result$comm_score <- setNames(
      colMeans(combined, na.rm = TRUE), cell_names)
    if (verbose) {
      q <- stats::quantile(result$comm_score, c(0.5, 0.9, 0.99))
      message(sprintf("[Score] Comm score quantiles: median=%.3f, 90th=%.3f, 99th=%.3f",
                      q[1], q[2], q[3]))
    }
  }

  if (verbose) message("[Score] Done.")
  result
}

#' Rank Cells by Communication Potential
#'
#' Ranks cells by their aggregate communication score and returns the top senders,
#' top receivers, and top communicators.
#'
#' @param activity_scores Output of \code{\link{score_lr_activity}}.
#' @param n Top cells to return per category. Default: \code{50}.
#' @param cell_labels Optional named vector mapping cell names to cluster labels.
#'
#' @return A data frame with cell, sender_score, receiver_score, comm_score, rank,
#'   and optionally cluster.
#' @export
rank_comm_cells <- function(activity_scores, n = 50, cell_labels = NULL) {
  stopifnot(is.list(activity_scores))
  if (is.null(activity_scores$comm_score)) {
    stop("activity_scores must contain comm_score. Run score_lr_activity(..., mode = 'both', aggregate = TRUE).")
  }

  cells <- names(activity_scores$comm_score)
  df <- data.frame(
    cell           = cells,
    sender_score   = activity_scores$sender_score[cells],
    receiver_score = activity_scores$receiver_score[cells],
    comm_score     = activity_scores$comm_score[cells],
    stringsAsFactors = FALSE
  )
  df <- df[order(-df$comm_score), ]
  df$rank <- seq_len(nrow(df))

  if (!is.null(cell_labels)) {
    cl_vec <- as.character(cell_labels)
    if (!is.null(names(cell_labels))) {
      cl_vec <- cell_labels[df$cell]
    }
    df$cluster <- cl_vec
  }

  utils::head(df, n)
}
