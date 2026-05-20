\name{summarize_celltype_communication}
\alias{summarize_celltype_communication}
\alias{celltype_comm_to_lcs}
\title{Cell-type-resolved LogicComm communication summaries}
\description{
Summarize ligand-receptor Logic Consensus Scores (LCS) by directed
sender-cell-type to receiver-cell-type relationships, and convert the resulting
features to named vectors for sample-level comparison.
}
\usage{
summarize_celltype_communication(reo_mat, cell_labels = NULL,
  seurat_obj = NULL, label_col = NULL, knn_mat = NULL,
  lr_db = lr_pairs_human, graph_name = NULL, lcs_threshold = 0.01,
  min_edges = 20, min_active_edges = 1,
  mode = c("auto", "neighborhood", "global"),
  remove_self_edges = TRUE, graph_symmetrize = c("none", "or", "max"),
  edge_weight_mode = c("binary", "weighted"),
  min_role_hub_quantile = 0.2, min_role_event_count = 5,
  include_self = TRUE, verbose = TRUE)

celltype_comm_to_lcs(ct_comm,
  level = c("celltype_lr", "celltype_pair", "pathway_pair"),
  metric = NULL, active_only = FALSE)
}
\arguments{
\item{reo_mat}{Binary REO matrix, genes x cells, or a \code{LogicCommREOResult} containing \code{$logic}.}
\item{cell_labels}{Named vector/factor mapping cells to cell types, or a Seurat object.}
\item{seurat_obj}{Optional Seurat object used for KNN graph extraction and/or labels.}
\item{label_col}{Optional Seurat metadata column used as cell type labels.}
\item{knn_mat}{Optional cells x cells KNN/SNN adjacency matrix.}
\item{lr_db}{LogicComm ligand-receptor database.}
\item{graph_name}{Optional Seurat graph name, for example \code{"RNA_nn"} or \code{"RNA_snn"}.}
\item{lcs_threshold}{Minimum cell-type LCS for an active call.}
\item{min_edges}{Minimum number of sender-to-receiver cell edges required for active calls.}
\item{min_active_edges}{Minimum active edges required for an active call.}
\item{mode}{\code{"auto"}, \code{"neighborhood"}, or \code{"global"}.}
\item{remove_self_edges}{Remove diagonal/self-loop KNN entries before scoring.}
\item{graph_symmetrize}{How to symmetrize the graph before scoring: \code{"none"}, \code{"or"}, or \code{"max"}.}
\item{edge_weight_mode}{\code{"binary"} treats each non-zero graph entry as one edge; \code{"weighted"} uses graph weights in LCS.}
\item{min_role_hub_quantile}{Positive-hub quantile below which cell types are labeled \code{Low-communication}.}
\item{min_role_event_count}{Minimum incoming plus outgoing active L-R events required for a non-low role label.}
\item{include_self}{Include same-cell-type communication pairs.}
\item{verbose}{Whether to print progress messages.}
\item{ct_comm}{A \code{LogicCommCellTypeComm} object returned by \code{summarize_celltype_communication()}.}
\item{level}{Feature level to export as a named vector.}
\item{metric}{Numeric column to extract. Defaults to \code{"lcs"} for celltype-LR rows and \code{"sum_lcs"} otherwise.}
\item{active_only}{Whether to keep only active rows. For group comparison, \code{FALSE} is usually recommended.}
}
\details{
Neighborhood mode treats each non-zero KNN/SNN entry as a directed sender cell to
receiver cell edge. In weighted mode, the numerator and denominator are weighted
by graph edge weights. The output reports both the graph-edge denominator
(\code{n_edges} and \code{edge_weight_sum}) and the active evidence
(\code{n_active_edges}, \code{active_edge_weight_sum}, \code{lcs}, and
\code{lcs_unweighted}).

The role table distinguishes three count concepts: active L-R event counts,
unique L-R counts, and active edge support. Role labels include low-communication
filtering and confidence values, so weak cell types are not forced into a sender,
receiver, mediator, or influencer label.
}
\value{
\code{summarize_celltype_communication()} returns a list of class
\code{LogicCommCellTypeComm} with \code{lr_table}, \code{pair_summary},
\code{pathway_summary}, \code{role_summary}, \code{celltype_sizes},
\code{adjacency_strength}, \code{adjacency_count}, and
\code{adjacency_active_edge_support}.
\code{celltype_comm_to_lcs()} returns a named numeric vector.
}
\examples{
\dontrun{
ct_comm <- summarize_celltype_communication(
  reo, cell_labels = celltype_labels, seurat_obj = seurat_obj,
  graph_name = "RNA_snn", edge_weight_mode = "weighted",
  graph_symmetrize = "max", lr_db = lr_pairs_human
)
plot_celltype_heatmap(ct_comm, metric = "sum_lcs")
vec <- celltype_comm_to_lcs(ct_comm, level = "celltype_lr")
}
}
