\name{IdentifyRankLogicConsensus}
\alias{IdentifyRankLogicConsensus}
\title{Identify rank-aware logic consensus evidence}
\description{Compute binary LCS plus continuous rank-dominance, rank-margin, stability, and specificity evidence.}
\usage{IdentifyRankLogicConsensus(reo_mat = NULL, rank_mat = NULL,
  expr_mat = NULL, seurat_obj = NULL, knn_mat = NULL, lr_db = lr_pairs_human,
  graph_name = NULL, rank_threshold = 0.5,
  threshold_grid = c(0.4, 0.5, 0.6, 0.7), cell_labels = NULL,
  layer = "counts", chunk_size = 5000, anchor_genes = NULL,
  remove_self_edges = TRUE, graph_symmetrize = c("or", "none", "max"),
  edge_weight_mode = c("binary", "weighted"),
  complex_aggregate = c("min", "mean", "max"), lcs_threshold = 0.01,
  verbose = TRUE)}
\arguments{
\item{reo_mat}{Optional binary REO matrix or \code{LogicCommREOResult}.}
\item{rank_mat}{Optional rank-percentile matrix.}
\item{expr_mat}{Optional expression input used to calculate ranks.}
\item{seurat_obj}{Optional Seurat object for graph extraction.}
\item{knn_mat}{Optional cells x cells adjacency matrix.}
\item{lr_db}{Ligand-receptor database.}
\item{graph_name}{Optional Seurat graph name.}
\item{rank_threshold}{Rank threshold for binary activity when \code{reo_mat} is absent.}
\item{threshold_grid}{Rank thresholds used for stability scoring.}
\item{cell_labels}{Optional named cell labels for specificity scoring.}
\item{layer}{Seurat layer used when \code{expr_mat} is provided.}
\item{chunk_size}{Chunk size used when \code{expr_mat} is provided.}
\item{anchor_genes}{Optional anchor genes used when \code{expr_mat} is provided.}
\item{remove_self_edges}{Whether to remove diagonal graph edges.}
\item{graph_symmetrize}{Graph symmetrization mode.}
\item{edge_weight_mode}{Whether graph edges are binary or weighted.}
\item{complex_aggregate}{How to combine multi-subunit ranks.}
\item{lcs_threshold}{Threshold used for evidence-tier labels.}
\item{verbose}{Whether to print progress messages.}
}
\details{This function complements \code{\link{IdentifyLogicConsensus}} by preserving rank strength and threshold stability instead of only reporting active frequency.}
\value{A data frame with one row per ligand-receptor pair.}
\examples{\dontrun{
reo <- calc_REO_matrix(count_matrix, lr_genes = all_lr_genes(lr_pairs_human),
  return_rank = TRUE)
rank_ev <- IdentifyRankLogicConsensus(reo)
}}
