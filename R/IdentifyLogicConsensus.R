# R/IdentifyLogicConsensus.R

#' Identify Logic Consensus Scores (LCS) for all L-R pairs
#'
#' For each ligand-receptor pair in the database, computes a Logic Consensus
#' Score (LCS) that measures how consistently the complete L-R signaling logic
#' is realized across the cell population (or its neighborhoods).
#'
#' @details
#' **Two scoring modes:**
#'
#' \describe{
#'   \item{\strong{Neighborhood mode} (recommended, requires KNN)}{
#'     For each directed KNN edge (sender → receiver):
#'     \deqn{LCS = \frac{\text{edges where sender has L=TRUE AND receiver has R=TRUE}}{\text{total edges}}}
#'     This captures the *spatial* logic of cell-cell communication — the ligand
#'     must be expressed in the sending cell, and the receptor in the (neighboring)
#'     receiving cell.
#'   }
#'   \item{\strong{Global mode} (no KNN required)}{
#'     \deqn{LCS = \frac{\text{cells with both L=TRUE AND R=TRUE}}{n_{cells}}}
#'     Simpler co-expression score; ignores spatial arrangement.
#'   }
#' }
#'
#' **Multimer complexes:** A complex like "TGFBR1+TGFBR2" requires ALL subunit
#' genes to be active (\code{TRUE}) for the receptor to be counted as expressed.
#'
#' @param reo_mat Binary REO matrix (genes x cells), output of
#'   \code{\link{calc_REO_matrix}}.
#' @param seurat_obj Optional Seurat object. If provided and \code{knn_mat} is
#'   \code{NULL}, the KNN graph is extracted from
#'   \code{seurat_obj@graphs} automatically.
#' @param knn_mat Optional KNN adjacency sparse matrix (cells x cells). Overrides
#'   automatic extraction. Non-zero entries indicate neighbor relationships.
#'   Can be a \code{dgCMatrix} or binary matrix.
#' @param lr_db Ligand-receptor database. Data frame with columns
#'   \code{lr_pair}, \code{ligand_genes}, \code{receptor_genes} (list columns).
#'   Defaults to the built-in \code{lr_pairs_human}.
#' @param graph_name Name of the KNN graph slot in \code{seurat_obj}. If
#'   \code{NULL}, auto-selected.
#' @param lcs_threshold Minimum LCS for a pair to be considered "active" in
#'   downstream comparisons. Default: \code{0.01}.
#' @param remove_self_edges Logical; if \code{TRUE}, diagonal/self-loop entries
#'   in KNN graphs are removed before scoring. Default: \code{TRUE}.
#' @param graph_symmetrize How to symmetrize the graph before scoring: \code{"none"}
#'   keeps the matrix as supplied; \code{"or"} uses binary union; \code{"max"}
#'   keeps the maximum weight in either direction. By default, binary KNN scoring
#'   uses \code{"or"} when self-loops are removed; weighted scoring and explicit
#'   self-loop scoring keep supplied directed weights.
#' @param edge_weight_mode \code{"binary"} treats every non-zero graph entry
#'   as one edge; \code{"weighted"} uses graph weights in the LCS denominator.
#' @param verbose Print progress messages. Default: \code{TRUE}.
#'
#' @return A named numeric vector of LCS values, one per L-R pair in
#'   \code{lr_db}. Names are \code{lr_pair} identifiers. Pairs where neither
#'   the ligand nor receptor is found in \code{reo_mat} get \code{NA}.
#'
#' @seealso \code{\link{calc_REO_matrix}}, \code{\link{CompareLogicGroups}}
#'
#' @examples
#' reo <- Matrix::Matrix(
#'   matrix(c(1, 0,
#'            0, 1), nrow = 2, byrow = TRUE,
#'          dimnames = list(c("L", "R"), c("C1", "C2"))),
#'   sparse = TRUE
#' )
#' lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
#'                     stringsAsFactors = FALSE)
#' lr_db$ligand_genes <- list("L")
#' lr_db$receptor_genes <- list("R")
#' lcs <- IdentifyLogicConsensus(reo, lr_db = lr_db, verbose = FALSE)
#' sort(lcs, decreasing = TRUE)
#'
#' @importFrom Matrix summary
#' @export
IdentifyLogicConsensus <- function(reo_mat,
                                    seurat_obj    = NULL,
                                    knn_mat       = NULL,
                                    lr_db         = lr_pairs_human,
                                    graph_name    = NULL,
                                    lcs_threshold = 0.01,
                                    remove_self_edges = TRUE,
                                    graph_symmetrize = c("or", "none", "max"),
                                    edge_weight_mode = c("binary", "weighted"),
                                    verbose       = TRUE) {

  stopifnot(inherits(reo_mat, c("dgCMatrix","REOMatrix","sparseMatrix","matrix")))
  edge_weight_mode <- match.arg(edge_weight_mode)
  if (missing(graph_symmetrize)) {
    graph_symmetrize <- if (edge_weight_mode == "binary" && isTRUE(remove_self_edges)) "or" else "none"
  } else {
    graph_symmetrize <- match.arg(graph_symmetrize)
  }

  # ── 1. Resolve KNN ─────────────────────────────────────────────────────────
  use_knn   <- FALSE
  knn_edges <- NULL

  if (!is.null(knn_mat)) {
    use_knn   <- TRUE
    knn_mat   <- .validate_and_align_knn(knn_mat, colnames(reo_mat))
    knn_mat   <- .symmetrize_sparse_graph(knn_mat, graph_symmetrize)
    knn_edges <- .sparse_to_edges(knn_mat, remove_self_edges = remove_self_edges,
                                  edge_weight_mode = edge_weight_mode)
  } else if (!is.null(seurat_obj)) {
    tryCatch({
      knn_mat   <- .extract_knn(seurat_obj, graph_name)
      knn_mat   <- .validate_and_align_knn(knn_mat, colnames(reo_mat))
      knn_mat   <- .symmetrize_sparse_graph(knn_mat, graph_symmetrize)
      knn_edges <- .sparse_to_edges(knn_mat, remove_self_edges = remove_self_edges,
                                    edge_weight_mode = edge_weight_mode)
      use_knn   <- TRUE
      if (verbose) message("[LCS] Neighborhood mode: KNN graph loaded (",
                           length(knn_edges$i), " edges).")
    }, error = function(e) {
      if (verbose) message("[LCS] KNN graph not found, falling back to global mode. ",
                           "(Run FindNeighbors() to enable neighborhood mode.)")
    })
  }

  if (use_knn && length(knn_edges$i) == 0) {
    if (verbose) message("[LCS] KNN graph has no usable edges after filtering; falling back to global mode.")
    use_knn <- FALSE
  }

  if (!use_knn && verbose) {
    message("[LCS] Global co-expression mode (no KNN).")
  }

  # ── 2. Pre-calculate complex logic ────────────────────────────────────────
  if (verbose) message("[LCS] Resolving complex logic for unique gene sets...")
  
  unique_complexes <- unique(c(lr_db$ligand_genes, lr_db$receptor_genes))
  complex_keys <- vapply(unique_complexes, function(gs) paste(sort(gs), collapse = "|"), character(1))
  
  # Map each unique complex to its logic vector
  complex_logic_map <- lapply(unique_complexes, function(gs) {
    .resolve_complex_logic(gs, reo_mat)
  })
  names(complex_logic_map) <- complex_keys
  
  # Helper to get pre-calculated logic
  get_logic <- function(gs) {
    key <- paste(sort(gs), collapse = "|")
    complex_logic_map[[key]]
  }

  # ── 3. Compute LCS per L-R pair ───────────────────────────────────────────
  n_pairs <- nrow(lr_db)
  lcs_vec <- setNames(rep(NA_real_, n_pairs), lr_db$lr_pair)

  if (verbose) message(sprintf("[LCS] Scoring %d L-R pairs...", n_pairs))

  for (i in seq_len(n_pairs)) {
    lig_genes <- lr_db$ligand_genes[[i]]
    rec_genes <- lr_db$receptor_genes[[i]]

    # Skip if no subunit found in REO matrix
    if (!any(lig_genes %in% rownames(reo_mat)) || !any(rec_genes %in% rownames(reo_mat))) {
      next
    }

    lig_logic <- get_logic(lig_genes)
    rec_logic <- get_logic(rec_genes)

    if (use_knn) {
      if (edge_weight_mode == "weighted") {
        lcs_vec[i] <- .compute_lcs_knn_weighted(lig_logic, rec_logic,
                                                knn_edges$i, knn_edges$j, knn_edges$w)
      } else {
        lcs_vec[i] <- compute_lcs_knn(lig_logic, rec_logic,
                                      knn_edges$i, knn_edges$j)
      }
    } else {
      lcs_vec[i] <- compute_lcs_global(lig_logic, rec_logic)
    }
  }

  n_active <- sum(!is.na(lcs_vec) & lcs_vec >= lcs_threshold, na.rm = TRUE)
  if (verbose) message(sprintf(
    "[LCS] Done. %d / %d pairs active (LCS >= %.3f).",
    n_active, n_pairs, lcs_threshold))

  class(lcs_vec) <- c("LCSVector", class(lcs_vec))
  lcs_vec
}

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Validate and align KNN adjacency matrix to REO cell order
#' @keywords internal
.validate_and_align_knn <- function(knn_mat, cell_names) {
  if (is.null(dim(knn_mat)) || nrow(knn_mat) != ncol(knn_mat)) {
    stop("knn_mat must be a square cells x cells adjacency matrix.")
  }

  if (!is.null(rownames(knn_mat)) && !is.null(colnames(knn_mat))) {
    missing_rows <- setdiff(cell_names, rownames(knn_mat))
    missing_cols <- setdiff(cell_names, colnames(knn_mat))
    if (length(missing_rows) > 0 || length(missing_cols) > 0) {
      stop("knn_mat row/column names must contain all reo_mat column names.")
    }
    knn_mat <- knn_mat[cell_names, cell_names, drop = FALSE]
  } else if (nrow(knn_mat) != length(cell_names)) {
    stop("knn_mat dimensions must match ncol(reo_mat) when names are absent.")
  } else {
    rownames(knn_mat) <- colnames(knn_mat) <- cell_names
  }

  if (anyNA(knn_mat)) stop("knn_mat must not contain NA values.")
  if (!inherits(knn_mat, "sparseMatrix")) {
    knn_mat <- Matrix::Matrix(knn_mat, sparse = TRUE)
  }
  knn_mat
}

#' Symmetrize a sparse graph while preserving optional edge weights
#' @keywords internal
.symmetrize_sparse_graph <- function(knn_mat, mode = c("none", "or", "max")) {
  mode <- match.arg(mode)
  if (!inherits(knn_mat, "sparseMatrix")) {
    knn_mat <- Matrix::Matrix(knn_mat, sparse = TRUE)
  }
  if (mode == "none") return(knn_mat)
  s <- Matrix::summary(knn_mat)
  if (nrow(s) == 0) return(knn_mat)
  df <- rbind(
    data.frame(i = as.integer(s$i), j = as.integer(s$j), x = as.numeric(s$x)),
    data.frame(i = as.integer(s$j), j = as.integer(s$i), x = as.numeric(s$x))
  )
  if (mode == "or") {
    df$x <- 1
    agg <- stats::aggregate(x ~ i + j, df, max)
  } else {
    agg <- stats::aggregate(x ~ i + j, df, max)
  }
  Matrix::sparseMatrix(
    i = agg$i, j = agg$j, x = agg$x,
    dims = dim(knn_mat), dimnames = dimnames(knn_mat)
  )
}

#' Convert a KNN matrix to directed edge index vectors (1-based)
#' @keywords internal
.sparse_to_edges <- function(knn_mat, remove_self_edges = TRUE,
                             edge_weight_mode = c("binary", "weighted")) {
  edge_weight_mode <- match.arg(edge_weight_mode)
  if (!inherits(knn_mat, "sparseMatrix")) {
    knn_mat <- Matrix::Matrix(knn_mat, sparse = TRUE)
  }
  if (isTRUE(remove_self_edges) && nrow(knn_mat) > 0) {
    diag(knn_mat) <- 0
    knn_mat <- Matrix::drop0(knn_mat)
  }
  s <- Matrix::summary(knn_mat)
  if (nrow(s) == 0) {
    return(list(i = integer(0), j = integer(0), w = numeric(0)))
  }
  w <- if (edge_weight_mode == "weighted") as.numeric(s$x) else rep(1, nrow(s))
  w[!is.finite(w)] <- 0
  keep <- w > 0
  list(i = as.integer(s$i[keep]), j = as.integer(s$j[keep]), w = as.numeric(w[keep]))
}

#' Weighted KNN LCS fallback used when weighted SNN/KNN edges are requested
#' @keywords internal
.compute_lcs_knn_weighted <- function(lig_logic, rec_logic, edge_i, edge_j, edge_w) {
  if (length(edge_i) == 0) return(NA_real_)
  edge_w <- as.numeric(edge_w)
  denom <- sum(edge_w, na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  active <- lig_logic[edge_i] & rec_logic[edge_j]
  sum(edge_w[active], na.rm = TRUE) / denom
}

#' Print method for LCSVector
#' @param x LCSVector.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}, called for side effects.
#' @examples
#' lcs <- structure(c(L_R = 0.25, X_Y = NA_real_), class = "LCSVector")
#' print(lcs)
#' @export
print.LCSVector <- function(x, ...) {
  n_total  <- length(x)
  n_scored <- sum(!is.na(x))
  n_active <- sum(x > 0.01, na.rm = TRUE)
  cat(sprintf(
    "LCSVector: %d pairs | scored: %d | active (>0.01): %d\n",
    n_total, n_scored, n_active))
  cat("Top 5:\n")
  top <- head(sort(x[!is.na(x)], decreasing = TRUE), 5)
  print(round(top, 4))
  invisible(x)
}
