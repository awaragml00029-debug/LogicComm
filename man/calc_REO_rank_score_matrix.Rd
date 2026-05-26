\name{calc_REO_rank_score_matrix}
\alias{calc_REO_rank_score_matrix}
\title{Calculate a REO rank score matrix}
\description{Return within-cell rank percentiles for ligand-receptor genes.}
\usage{calc_REO_rank_score_matrix(expr_mat, genes = NULL, lr_db = lr_pairs_human,
  layer = "counts", chunk_size = 5000, anchor_genes = NULL, verbose = TRUE)}
\arguments{
\item{expr_mat}{Expression input accepted by \code{\link{calc_REO_matrix}}.}
\item{genes}{Optional genes to retain. If \code{NULL}, genes are taken from \code{lr_db}.}
\item{lr_db}{Ligand-receptor database used when \code{genes = NULL}.}
\item{layer}{Seurat layer or slot.}
\item{chunk_size}{Number of cells processed per chunk.}
\item{anchor_genes}{Optional anchor gene universe.}
\item{verbose}{Whether to print progress messages.}
}
\details{The output keeps continuous within-cell rank information that is lost when REO is binarized.}
\value{Dense numeric genes x cells matrix with values in \code{[0, 1]}.}
\examples{\dontrun{
rank_mat <- calc_REO_rank_score_matrix(count_matrix, lr_db = lr_pairs_human)
}}
