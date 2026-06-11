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

  # Accept a LogicCommREOResult (the output of calc_REO_matrix(..., return_rank =
  # TRUE)) by using its binary logic matrix, mirroring the other scoring entry
  # points. Then require a 2-D genes x cells matrix with names: a single-row or
  # single-column slice that dropped to a vector (subset without drop = FALSE), a
  # 1-D array, or a bare list would otherwise fail later inside
  # .resolve_complex_logic() with a cryptic "invalid 'times' argument".
  if (is.list(reo_mat) && !is.null(reo_mat$logic)) reo_mat <- reo_mat$logic
  if (is.null(dim(reo_mat)) || length(dim(reo_mat)) != 2L) {
    detail <- if (is.null(dim(reo_mat))) {
      sprintf("a %s of length %d (a single-gene or single-cell slice drops to a vector; subset with drop = FALSE)",
              class(reo_mat)[1], length(reo_mat))
    } else {
      sprintf("a %s with %d dimension(s)", class(reo_mat)[1], length(dim(reo_mat)))
    }
    stop("reo_mat must be a 2D genes x cells REO matrix; got ", detail,
         ". If you built it with calc_REO_matrix(..., return_rank = TRUE), pass its $logic element.",
         call. = FALSE)
  }
  if (is.null(rownames(reo_mat)) || is.null(colnames(reo_mat))) {
    stop("reo_mat must have gene rownames and cell colnames.", call. = FALSE)
  }

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
    if (!is.null(names(cell_labels))) {
      cl_vec <- as.character(cell_labels)[match(df$cell, names(cell_labels))]
    } else {
      # df has been reordered by comm_score, so unnamed labels (assumed to be in
      # the original scored-cell order) must be aligned to df$cell rather than
      # pasted on in sorted order.
      if (length(cell_labels) != length(cells)) {
        stop("Unnamed cell_labels must have length equal to the number of scored cells.")
      }
      cl_vec <- stats::setNames(as.character(cell_labels), cells)[df$cell]
    }
    df$cluster <- unname(cl_vec)
  }

  utils::head(df, n)
}

#' Score Receiver Downstream Response Evidence
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param reo_mat Binary REO matrix or LogicCommREOResult.
#' @param response_db Data frame with \code{lr_pair} and \code{response_genes}.
#' @param response_mode How to combine response genes: \code{"any"} or \code{"all"}.
#' @return A data frame with receiver-response evidence per L-R row.
#' @export
score_receiver_response <- function(ct_comm,
                                    reo_mat,
                                    response_db,
                                    response_mode = c("any", "all")) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  response_mode <- match.arg(response_mode)
  if (is.list(reo_mat) && !is.null(reo_mat$logic)) reo_mat <- reo_mat$logic
  if (!all(c("lr_pair", "response_genes") %in% names(response_db))) {
    stop("response_db must contain lr_pair and response_genes columns.")
  }
  # Align labels and reo_mat to their shared cells. ct_comm$cell_labels may cover
  # fewer cells than reo_mat (summarize_celltype_communication() filters cells
  # with missing/empty labels), and reo_mat may carry extra cells; indexing
  # cell_labels by colnames(reo_mat) would otherwise inject NA labels (and NA
  # names) that propagate into the receiver-cell selection.
  labels <- ct_comm$cell_labels
  if (is.null(names(labels))) stop("ct_comm$cell_labels must be named by cell.")
  common <- intersect(names(labels), colnames(reo_mat))
  if (!length(common)) {
    stop("reo_mat shares no cells with ct_comm$cell_labels; recompute on the same cells.")
  }
  labels <- labels[common]
  reo_mat <- reo_mat[, common, drop = FALSE]
  response_map <- stats::setNames(lapply(as.character(response_db$response_genes), .parse_response_genes), as.character(response_db$lr_pair))
  df <- ct_comm$lr_table
  response_score <- numeric(nrow(df))
  response_gene_count <- integer(nrow(df))
  response_active_frac <- numeric(nrow(df))

  for (i in seq_len(nrow(df))) {
    genes <- response_map[[df$lr_pair[i]]]
    genes <- intersect(genes %||% character(0), rownames(reo_mat))
    response_gene_count[i] <- length(genes)
    receiver_cells <- names(labels)[labels == df$receiver_type[i]]
    if (length(genes) == 0 || length(receiver_cells) == 0) {
      response_score[i] <- NA_real_
      response_active_frac[i] <- NA_real_
      next
    }
    mat <- reo_mat[genes, receiver_cells, drop = FALSE]
    if (inherits(mat, "sparseMatrix")) {
      sums <- Matrix::colSums(mat)
    } else {
      sums <- colSums(as.matrix(mat), na.rm = TRUE)
    }
    active <- if (response_mode == "all") sums == length(genes) else sums > 0
    response_active_frac[i] <- mean(active, na.rm = TRUE)
    response_score[i] <- response_active_frac[i]
  }

  out <- df
  out$response_gene_count <- response_gene_count
  out$response_active_frac <- response_active_frac
  out$receiver_response_score <- response_score
  out$response_integrated_score <- out$lcs * ifelse(is.na(response_score), 0, response_score)
  out
}

#' Add Receiver Response Scores to a Cell-Type Communication Object
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param reo_mat Binary REO matrix or LogicCommREOResult.
#' @param response_db Data frame with \code{lr_pair} and \code{response_genes}.
#' @param response_mode How to combine response genes: \code{"any"} or \code{"all"}.
#' @return Updated \code{LogicCommCellTypeComm} object.
#' @export
add_receiver_response_score <- function(ct_comm,
                                        reo_mat,
                                        response_db,
                                        response_mode = c("any", "all")) {
  response_mode <- match.arg(response_mode)
  scored <- score_receiver_response(ct_comm, reo_mat, response_db, response_mode = response_mode)
  cols <- c("response_gene_count", "response_active_frac", "receiver_response_score", "response_integrated_score")
  ct_comm$lr_table[, cols] <- scored[, cols]
  ct_comm
}

.parse_response_genes <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("\\s+", "", x)
  out <- unlist(strsplit(x, "[;,|+]", perl = TRUE), use.names = FALSE)
  unique(out[nzchar(out) & !is.na(out)])
}
