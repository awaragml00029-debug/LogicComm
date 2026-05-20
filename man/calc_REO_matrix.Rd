\name{calc_REO_matrix}
\alias{calc_REO_matrix}
\title{Calculate the Relative Expression Ordering matrix}
\description{Convert a genes x cells expression matrix into a sparse binary REO logic matrix.}
\usage{calc_REO_matrix(expr_mat, lr_genes = NULL, rank_threshold = 0.5,
  layer = "counts", chunk_size = 5000, return_rank = FALSE,
  rank_output = c("percentile", "none"), anchor_genes = NULL,
  verbose = TRUE)}
\arguments{
\item{expr_mat}{Expression input: matrix, sparse \pkg{Matrix}, Seurat object, or BPCells IterableMatrix.}
\item{lr_genes}{Optional genes retained in the output matrix.}
\item{rank_threshold}{Per-cell quantile used as the activity anchor.}
\item{layer}{Seurat layer or slot to extract.}
\item{chunk_size}{Number of cells processed per chunk.}
\item{return_rank}{Whether to return rank percentiles in addition to logic states.}
\item{rank_output}{Rank output type.}
\item{anchor_genes}{Optional gene universe used to compute anchors; \code{NULL} uses all genes.}
\item{verbose}{Whether to print progress messages.}
}
\details{By default, anchors are computed from all genes in \code{expr_mat}; \code{lr_genes} only controls retained output rows.}
\value{A sparse logic matrix, or a \code{LogicCommREOResult} list with \code{logic} and \code{rank}.}
\examples{\dontrun{
data(lr_pairs_human)
lr_genes <- all_lr_genes(lr_pairs_human)
reo <- calc_REO_matrix(count_matrix, lr_genes = lr_genes)
}}
