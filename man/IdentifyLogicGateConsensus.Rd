\name{IdentifyLogicGateConsensus}
\alias{IdentifyLogicGateConsensus}
\title{Identify rank-aware logic-gate consensus scores}
\description{Extend base LCS with optional CellChatDB agonist, antagonist, co-activating receptor, and co-inhibitory receptor logic gates.}
\usage{IdentifyLogicGateConsensus(reo_mat, rank_mat = NULL, seurat_obj = NULL,
  knn_mat = NULL, lr_db = lr_pairs_human, graph_name = NULL,
  positive_gate = c("all", "any", "ignore"),
  negative_gate = c("rank_block", "active_block", "ignore"),
  agonist_side = c("receiver", "sender", "either", "ignore"),
  antagonist_side = c("both", "sender", "receiver", "either", "ignore"),
  lcs_threshold = 0.01, remove_self_edges = TRUE,
  graph_symmetrize = c("none", "or", "max"),
  edge_weight_mode = c("binary", "weighted"), verbose = TRUE)}
\arguments{
\item{reo_mat}{Binary REO matrix or \code{LogicCommREOResult}.}
\item{rank_mat}{Optional rank-percentile matrix. Required for \code{negative_gate = "rank_block"} unless included in \code{reo_mat}.}
\item{seurat_obj}{Optional Seurat object from which to extract a graph.}
\item{knn_mat}{Optional cells x cells KNN/SNN adjacency matrix.}
\item{lr_db}{Ligand-receptor database with optional modulator list columns.}
\item{graph_name}{Optional graph name for Seurat objects.}
\item{positive_gate}{How positive modulators are handled.}
\item{negative_gate}{How negative modulators are handled.}
\item{agonist_side}{Sender/receiver side for agonist gates.}
\item{antagonist_side}{Sender/receiver side for antagonist gates.}
\item{lcs_threshold}{Threshold used for summaries.}
\item{remove_self_edges}{Whether to remove diagonal self-loop graph entries.}
\item{graph_symmetrize}{How to symmetrize graph entries before scoring.}
\item{edge_weight_mode}{\code{"binary"} or \code{"weighted"}.}
\item{verbose}{Whether to print progress messages.}
}
\value{A data frame of class \code{LogicGateResult}.}
\examples{\dontrun{gate <- IdentifyLogicGateConsensus(reo_result, lr_db = cellchat_lr,
  edge_weight_mode = "weighted")}}
