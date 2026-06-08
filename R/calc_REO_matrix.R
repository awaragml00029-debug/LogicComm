# R/calc_REO_matrix.R

#' Calculate the Relative Expression Ordering (REO) Matrix
#'
#' Converts a gene expression matrix into a binary logical matrix by comparing
#' each retained gene's expression in each cell against that cell's dynamic
#' background anchor. A gene is considered "actively expressed" (\code{TRUE})
#' when its expression exceeds the cell's anchor, which is a per-cell quantile
#' computed over the anchor gene universe.
#'
#' @details
#' **Why REO?** Raw count data cannot be compared across samples or even across
#' cells in a single sample without normalization. REO sidesteps this by using a
#' *within-cell* relative ordering: a gene whose expression ranks above the chosen
#' within-cell quantile is treated as active, regardless of sequencing depth.
#'
#' **Important implementation detail:** by default the anchor is computed over
#' all genes in the expression matrix, while \code{lr_genes} only controls which
#' genes are returned in the final REO matrix. This avoids making the binary state
#' depend on the size or composition of a ligand-receptor database.
#'
#' **Multimer handling:** Downstream L-R scoring treats composite
#' receptors/ligands (e.g. TGFBR1+TGFBR2) as active only when *all* subunit genes
#' are above anchor (logical AND).
#'
#' @param expr_mat Gene expression matrix (genes x cells). Accepts a numeric
#'   matrix, a sparse \code{Matrix}, a \code{Seurat} object, or a \code{BPCells}
#'   \code{IterableMatrix}.
#' @param lr_genes Character vector of genes to retain in the output. If
#'   \code{NULL}, all genes in the matrix are retained.
#' @param rank_threshold Per-cell quantile used as the expression anchor.
#'   Default: \code{0.5} (median of expressed genes). Increase to be more
#'   stringent (e.g. \code{0.75}).
#' @param layer For Seurat objects: which assay layer/slot to use. Default:
#'   \code{"counts"}.
#' @param chunk_size Number of cells to process per chunk. Default: \code{5000}.
#' @param return_rank Logical; if \code{TRUE}, return a list containing both the
#'   sparse binary REO matrix and a dense within-cell rank-percentile matrix for
#'   retained genes.
#' @param rank_output Rank output type. \code{"percentile"} returns values in
#'   \code{[0, 1]} where larger values indicate higher within-cell rank;
#'   \code{"none"} suppresses rank output even when \code{return_rank = TRUE}.
#' @param anchor_genes Optional character vector defining the gene universe used
#'   to compute per-cell anchors and rank percentiles. The default \code{NULL}
#'   uses all genes in \code{expr_mat}. Use with care: anchors computed on a small
#'   subset are less comparable across databases or studies.
#' @param gene_background Cell-type-free ambient guard. \code{"none"} (default)
#'   keeps the original within-cell REO rule. \code{"quantile"} additionally
#'   requires each gene's count to exceed its own background, taken as the
#'   \code{gene_background_quantile} quantile of that gene over the cells that
#'   detect it (nonzero counts). This removes ambient/dropout counts that clear
#'   the within-cell anchor but are low for the gene itself (for example ambient
#'   \code{CD8A} in non-CD8 cells), and stays effective even when the gene is
#'   detected in fewer than half of cells. It needs no cell-type annotation. For
#'   rigorous ambient removal, prefer an upstream tool such as SoupX, DecontX, or
#'   CellBender.
#' @param gene_background_quantile Quantile of each gene's detected (nonzero)
#'   counts used as its background when \code{gene_background = "quantile"}.
#'   Default: \code{0.5} (the median detected count).
#' @param verbose Print progress messages. Default: \code{TRUE}.
#'
#' @return By default, a sparse \code{dgCMatrix}-like matrix (genes x cells) where
#'   stored 1 values indicate active expression above the cell's anchor. If
#'   \code{return_rank = TRUE}, a list with \code{logic} and \code{rank} matrices
#'   is returned and given class \code{LogicCommREOResult}.
#'
#' @seealso \code{\link{IdentifyLogicConsensus}}, \code{\link{all_lr_genes}}
#'
#' @examples
#' expr <- matrix(
#'   c(10, 0,
#'     0, 10,
#'     1, 1),
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(c("L", "R", "BG"), c("C1", "C2"))
#' )
#' reo <- calc_REO_matrix(expr, lr_genes = c("L", "R"), verbose = FALSE)
#' reo
#'
#' @export
calc_REO_matrix <- function(expr_mat,
                            lr_genes       = NULL,
                            rank_threshold = 0.5,
                            layer          = "counts",
                            chunk_size     = 5000,
                            return_rank    = FALSE,
                            rank_output    = c("percentile", "none"),
                            anchor_genes   = NULL,
                            gene_background = c("none", "quantile"),
                            gene_background_quantile = 0.5,
                            verbose        = TRUE) {
  rank_output <- match.arg(rank_output)
  gene_background <- match.arg(gene_background)
  stopifnot(is.numeric(rank_threshold), length(rank_threshold) == 1,
            rank_threshold > 0, rank_threshold < 1)
  stopifnot(is.numeric(gene_background_quantile), length(gene_background_quantile) == 1,
            gene_background_quantile >= 0, gene_background_quantile < 1)
  stopifnot(is.numeric(chunk_size), length(chunk_size) == 1, chunk_size >= 1)
  chunk_size <- as.integer(chunk_size)

  is_bpcells <- inherits(expr_mat, "IterableMatrix")
  if (is_bpcells) {
    if (!requireNamespace("BPCells", quietly = TRUE)) {
      stop("BPCells package required for IterableMatrix input.")
    }
    return(.calc_REO_bpcells(expr_mat, lr_genes, rank_threshold, chunk_size,
                             return_rank, rank_output, anchor_genes,
                             gene_background, gene_background_quantile, verbose))
  }

  mat_full <- .extract_matrix(expr_mat, layer = layer)
  if (is.null(rownames(mat_full)) || is.null(colnames(mat_full))) {
    stop("expr_mat must have gene rownames and cell colnames.")
  }

  output_genes <- .select_output_genes(rownames(mat_full), lr_genes, "matrix", verbose)
  anchor_gene_set <- .select_anchor_genes(rownames(mat_full), anchor_genes, "matrix")

  n_genes <- length(output_genes)
  n_cells <- ncol(mat_full)
  if (verbose) message(sprintf(
    "[REO] Computing REO for %d retained genes x %d cells (anchors from %d genes; %.0f%% quantile%s).",
    n_genes, n_cells, length(anchor_gene_set), rank_threshold * 100,
    if (gene_background == "quantile") sprintf("; gene-background gate q=%.2f", gene_background_quantile) else ""))

  # Cell-type-free ambient guard: a gene is active in a cell only if it is also
  # above its own across-cell background (quantile over all cells). This removes
  # ambient/dropout counts that clear the within-cell anchor but are low for the
  # gene itself, without requiring any cell-type annotation.
  gene_bg <- if (gene_background == "quantile") {
    .gene_background_threshold(mat_full, output_genes, gene_background_quantile)
  } else NULL

  n_chunks <- ceiling(n_cells / chunk_size)
  result_cols <- vector("list", n_chunks)
  make_rank <- isTRUE(return_rank) && rank_output == "percentile"
  rank_cols <- if (make_rank) vector("list", n_chunks) else NULL

  for (ch in seq_len(n_chunks)) {
    s <- (ch - 1) * chunk_size + 1
    e <- min(ch * chunk_size, n_cells)

    anchor_chunk <- as.matrix(mat_full[anchor_gene_set, s:e, drop = FALSE])
    anchors <- cell_quantiles(anchor_chunk, prob = rank_threshold)

    expr_chunk <- as.matrix(mat_full[output_genes, s:e, drop = FALSE])
    logic_chunk <- fast_logic_scan(expr_chunk, anchors)
    if (!is.null(gene_bg)) logic_chunk <- logic_chunk & (expr_chunk > gene_bg)
    result_cols[[ch]] <- logic_chunk
    if (make_rank) rank_cols[[ch]] <- .rank_percentile_subset(anchor_chunk, output_genes)
  }

  logic_mat <- do.call(cbind, result_cols)
  rownames(logic_mat) <- output_genes
  colnames(logic_mat) <- colnames(mat_full)

  sparse_out <- Matrix::Matrix(logic_mat * 1L, sparse = TRUE)
  attr(sparse_out, "logic_type") <- "REO"
  attr(sparse_out, "anchor_n_genes") <- length(anchor_gene_set)
  attr(sparse_out, "anchor_source") <- if (is.null(anchor_genes)) "all_genes" else "anchor_genes"
  attr(sparse_out, "gene_background") <- gene_background

  if (verbose) {
    active_frac <- mean(logic_mat)
    message(sprintf("[REO] Done. Active fraction: %.1f%% of retained gene-cell entries.",
                    active_frac * 100))
  }

  if (isTRUE(return_rank)) {
    rank_mat <- NULL
    if (make_rank) {
      rank_mat <- do.call(cbind, rank_cols)
      rownames(rank_mat) <- output_genes
      colnames(rank_mat) <- colnames(mat_full)
    }
    return(structure(list(logic = sparse_out, rank = rank_mat),
                     class = "LogicCommREOResult"))
  }
  sparse_out
}

#' Select output genes for REO calculation
#' @keywords internal
.select_output_genes <- function(all_genes, lr_genes, label = "matrix", verbose = TRUE) {
  if (is.null(lr_genes)) return(all_genes)
  lr_genes <- unique(as.character(lr_genes))
  lr_genes <- lr_genes[nzchar(lr_genes) & !is.na(lr_genes)]
  keep <- intersect(all_genes, lr_genes)
  if (length(keep) == 0) stop("None of lr_genes found in ", label, " rownames.")
  if (verbose) message(sprintf("[REO] Retaining %d / %d requested genes.",
                               length(keep), length(lr_genes)))
  keep
}

#' Select anchor genes for REO calculation
#' @keywords internal
.select_anchor_genes <- function(all_genes, anchor_genes = NULL, label = "matrix") {
  if (is.null(anchor_genes)) return(all_genes)
  anchor_genes <- unique(as.character(anchor_genes))
  anchor_genes <- anchor_genes[nzchar(anchor_genes) & !is.na(anchor_genes)]
  keep <- intersect(all_genes, anchor_genes)
  if (length(keep) == 0) stop("None of anchor_genes found in ", label, " rownames.")
  keep
}

#' Per-gene across-cell background threshold (cell-type-free ambient guard)
#'
#' Returns, for each retained gene, the \code{prob} quantile of its expression
#' across all cells. For large matrices the cells are deterministically thinned
#' to \code{max_cells} so the estimate stays memory-safe. No cell-type labels are
#' used.
#' @keywords internal
.gene_background_threshold <- function(mat_full, output_genes, prob, max_cells = 30000L) {
  n_cells <- ncol(mat_full)
  idx <- if (n_cells > max_cells) unique(round(seq(1, n_cells, length.out = max_cells))) else seq_len(n_cells)
  sub <- as.matrix(mat_full[output_genes, idx, drop = FALSE])
  # Threshold from the cells that DETECT the gene (nonzero), not all cells:
  # a median over all cells collapses to 0 whenever a gene is detected in fewer
  # than half of cells (the common case), which would disable the guard.
  bg <- apply(sub, 1, function(x) {
    nz <- x[x > 0]
    if (!length(nz)) return(0)
    as.numeric(stats::quantile(nz, probs = prob, names = FALSE, type = 7))
  })
  bg <- as.numeric(bg)
  bg[!is.finite(bg)] <- 0
  bg
}

#' Rank retained genes against a within-cell universe
#' @keywords internal
.rank_percentile_subset <- function(universe_chunk, output_genes) {
  out <- matrix(NA_real_, nrow = length(output_genes), ncol = ncol(universe_chunk),
                dimnames = list(output_genes, colnames(universe_chunk)))
  if (nrow(universe_chunk) == 0 || ncol(universe_chunk) == 0) return(out)

  ranked <- apply(universe_chunk, 2, function(x) {
    if (length(x) == 1) return(1)
    rank(x, ties.method = "average") / length(x)
  })
  if (is.null(dim(ranked))) ranked <- matrix(ranked, ncol = 1)
  rownames(ranked) <- rownames(universe_chunk)
  colnames(ranked) <- colnames(universe_chunk)

  present <- intersect(output_genes, rownames(ranked))
  if (length(present) > 0) out[present, ] <- ranked[present, , drop = FALSE]
  out
}

#' BPCells streaming variant
#' @keywords internal
.calc_REO_bpcells <- function(bp_mat, lr_genes, rank_threshold, chunk_size,
                              return_rank = FALSE, rank_output = "percentile",
                              anchor_genes = NULL, gene_background = "none",
                              gene_background_quantile = 0.5, verbose = TRUE) {
  all_genes <- rownames(bp_mat)
  all_cells <- colnames(bp_mat)
  if (is.null(all_genes) || is.null(all_cells)) {
    stop("BPCells matrix must have gene rownames and cell colnames.")
  }

  output_genes <- .select_output_genes(all_genes, lr_genes, "BPCells matrix", verbose)
  anchor_gene_set <- .select_anchor_genes(all_genes, anchor_genes, "BPCells matrix")

  n_cells <- length(all_cells)
  n_chunks <- ceiling(n_cells / chunk_size)
  if (verbose) message(sprintf(
    "[REO/BPCells] Streaming %d chunks for %d cells; retaining %d genes; anchors from %d genes%s.",
    n_chunks, n_cells, length(output_genes), length(anchor_gene_set),
    if (gene_background == "quantile") sprintf("; gene-background gate q=%.2f", gene_background_quantile) else ""))

  gene_bg <- if (gene_background == "quantile") {
    bg <- .gene_background_threshold(bp_mat, output_genes, gene_background_quantile)
    stats::setNames(bg, output_genes)
  } else NULL

  result_cols <- vector("list", n_chunks)
  make_rank <- isTRUE(return_rank) && rank_output == "percentile"
  rank_cols <- if (make_rank) vector("list", n_chunks) else NULL

  for (ch in seq_len(n_chunks)) {
    s <- (ch - 1) * chunk_size + 1
    e <- min(ch * chunk_size, n_cells)

    anchor_chunk <- as.matrix(bp_mat[anchor_gene_set, s:e])
    rownames(anchor_chunk) <- anchor_gene_set
    colnames(anchor_chunk) <- all_cells[s:e]
    anchors <- cell_quantiles(anchor_chunk, prob = rank_threshold)

    expr_chunk <- as.matrix(bp_mat[output_genes, s:e])
    rownames(expr_chunk) <- output_genes
    colnames(expr_chunk) <- all_cells[s:e]
    logic_chunk <- fast_logic_scan(expr_chunk, anchors)
    if (!is.null(gene_bg)) logic_chunk <- logic_chunk & (expr_chunk > gene_bg)
    result_cols[[ch]] <- logic_chunk
    if (make_rank) rank_cols[[ch]] <- .rank_percentile_subset(anchor_chunk, output_genes)
  }

  logic_mat <- do.call(cbind, result_cols)
  rownames(logic_mat) <- output_genes
  colnames(logic_mat) <- all_cells

  out <- Matrix::Matrix(logic_mat * 1L, sparse = TRUE)
  attr(out, "logic_type") <- "REO"
  attr(out, "anchor_n_genes") <- length(anchor_gene_set)
  attr(out, "anchor_source") <- if (is.null(anchor_genes)) "all_genes" else "anchor_genes"
  attr(out, "gene_background") <- gene_background

  if (isTRUE(return_rank)) {
    rank_mat <- NULL
    if (make_rank) {
      rank_mat <- do.call(cbind, rank_cols)
      rownames(rank_mat) <- output_genes
      colnames(rank_mat) <- all_cells
    }
    return(structure(list(logic = out, rank = rank_mat),
                     class = "LogicCommREOResult"))
  }
  out
}

#' Print method for LogicCommREOResult
#' @param x LogicCommREOResult object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}, called for side effects.
#' @examples
#' reo <- structure(
#'   list(logic = Matrix::Matrix(1, nrow = 1, ncol = 1, sparse = TRUE),
#'        rank = matrix(0.8, nrow = 1, ncol = 1)),
#'   class = "LogicCommREOResult"
#' )
#' print(reo)
#' @export
print.LogicCommREOResult <- function(x, ...) {
  rank_msg <- if (is.null(x$rank)) "rank not returned" else
    sprintf("rank %d genes x %d cells", nrow(x$rank), ncol(x$rank))
  cat(sprintf("LogicCommREOResult: logic %d genes x %d cells | %s\n",
              nrow(x$logic), ncol(x$logic), rank_msg))
  invisible(x)
}
