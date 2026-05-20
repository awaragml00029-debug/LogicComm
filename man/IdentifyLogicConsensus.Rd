\name{IdentifyLogicConsensus}
\alias{IdentifyLogicConsensus}
\title{Identify Logic Consensus Scores}
\description{Compute ligand-receptor Logic Consensus Scores (LCS) for one sample.}
\usage{IdentifyLogicConsensus(reo_mat, seurat_obj = NULL, knn_mat = NULL,
  lr_db = lr_pairs_human, graph_name = NULL, lcs_threshold = 0.01,
  remove_self_edges = TRUE, graph_symmetrize = c("none", "or", "max"),
  edge_weight_mode = c("binary", "weighted"), verbose = TRUE)}
\arguments{
\item{reo_mat}{Binary REO matrix returned by \code{calc_REO_matrix()}.}
\item{seurat_obj}{Optional Seurat object from which to extract a graph.}
\item{knn_mat}{Optional cells x cells KNN/SNN adjacency matrix.}
\item{lr_db}{Ligand-receptor database in LogicComm format.}
\item{graph_name}{Optional graph name for Seurat objects.}
\item{lcs_threshold}{Threshold used for active-pair summaries.}
\item{remove_self_edges}{Whether to remove diagonal self-loop graph entries.}
\item{graph_symmetrize}{How to symmetrize graph entries before scoring: \code{"none"}, \code{"or"}, or \code{"max"}.}
\item{edge_weight_mode}{\code{"binary"} treats each non-zero graph entry as one edge; \code{"weighted"} uses graph weights.}
\item{verbose}{Whether to print progress messages.}
}
\details{
In neighborhood mode, a graph edge is interpreted as a potential sender to
receiver neighborhood. In weighted mode, LCS is the sum of active edge weights
divided by the total edge weight.
}
\value{A named numeric vector of class \code{LCSVector}.}
\examples{\dontrun{lcs <- IdentifyLogicConsensus(reo, lr_db = lr_pairs_human,
  graph_symmetrize = "max", edge_weight_mode = "weighted")}}
