\name{publication_extensions}
\alias{extract_spatial_coordinates}
\alias{build_spatial_graph}
\alias{summarize_spatial_communication}
\alias{plot_spatial_logic}
\alias{summarize_communication_dynamics}
\alias{plot_communication_dynamics}
\alias{fit_celltype_comm_glm}
\alias{plot_celltype_glm_volcano}
\alias{summarize_communication_findings}
\alias{write_communication_report}
\alias{plot_communication_qc}
\alias{plot_role_confidence}
\alias{print.LogicCommDynamics}
\alias{print.LogicCommFindings}
\alias{print.LogicCommGLMResult}
\title{Publication-Level LogicComm Extensions}
\description{
Spatial neighborhood construction, pseudotime/time-bin communication dynamics,
sample-level quasibinomial differential modeling, and publication-oriented
summary/report helpers for LogicComm cell-type communication results.
}
\usage{
extract_spatial_coordinates(object, coord_cols = NULL, reduction = NULL, image = NULL)

build_spatial_graph(coords, coord_cols = NULL, mode = c("knn", "radius"),
  k = 6, radius = NULL, directed = FALSE,
  distance_weight = c("binary", "inverse", "gaussian"), sigma = NULL,
  verbose = TRUE)

summarize_spatial_communication(reo_mat, coords, cell_labels,
  lr_db = lr_pairs_human, coord_cols = NULL, spatial_mode = c("knn", "radius"),
  k = 6, radius = NULL, directed = FALSE,
  distance_weight = c("binary", "inverse", "gaussian"), sigma = NULL,
  graph_symmetrize = c("none", "or", "max"), edge_weight_mode = NULL,
  verbose = TRUE, ...)

plot_spatial_logic(reo_mat, coords, lr_pair, lr_db = lr_pairs_human,
  coord_cols = NULL, cell_labels = NULL, sender_type = NULL, receiver_type = NULL,
  pt_size = 1.2, title = NULL)

summarize_communication_dynamics(reo_mat, pseudotime, cell_labels,
  lr_db = lr_pairs_human, knn_mat = NULL, n_bins = 5,
  bin_method = c("quantile", "equal_width", "categorical"),
  min_cells_per_bin = 20, verbose = TRUE, ...)

plot_communication_dynamics(dynamics, level = c("pair", "role", "lr"),
  metric = NULL, sender = NULL, receiver = NULL, lr_pair = NULL,
  cell_type = NULL, top_n = 10, title = NULL)

fit_celltype_comm_glm(sample_ct_list, sample_metadata, design = ~group,
  coef = NULL, level = c("celltype_lr", "celltype_pair"), min_samples = 3,
  min_total_trials = 20, fdr_method = "BH", verbose = TRUE)

plot_celltype_glm_volcano(glm_result, top_n_label = 15, title = NULL)

summarize_communication_findings(ct_comm, top_n = 10)

write_communication_report(ct_comm, file = "LogicComm_report.md",
  title = "LogicComm communication report", top_n = 10)

plot_communication_qc(ct_comm, title = NULL)

plot_role_confidence(ct_comm, title = NULL)
}
\arguments{
\item{object, coords}{Coordinate matrix/data frame or a Seurat object.}
\item{coord_cols}{Optional two coordinate columns.}
\item{reduction, image}{Optional Seurat reduction or image source.}
\item{mode, spatial_mode}{Spatial graph mode: k-nearest neighbors or fixed radius.}
\item{k}{Number of spatial nearest neighbors.}
\item{radius}{Radius for fixed-radius spatial graphs.}
\item{directed}{Whether to keep directed spatial edges.}
\item{distance_weight}{How to convert spatial distance into edge weights.}
\item{sigma}{Gaussian distance kernel bandwidth.}
\item{reo_mat}{Binary REO matrix or LogicCommREOResult.}
\item{cell_labels}{Named cell-type labels.}
\item{lr_db}{LogicComm ligand-receptor database.}
\item{graph_symmetrize, edge_weight_mode}{Graph handling options passed to cell-type communication scoring.}
\item{lr_pair}{Ligand-receptor pair to visualize.}
\item{sender_type, receiver_type}{Optional cell-type highlights for spatial plots.}
\item{pseudotime}{Named numeric pseudotime or categorical time/stage vector.}
\item{knn_mat}{Optional KNN/SNN/spatial graph.}
\item{n_bins, bin_method, min_cells_per_bin}{Time-bin construction controls.}
\item{dynamics}{Output from \code{summarize_communication_dynamics()}.}
\item{level}{Feature level for dynamics or GLM.}
\item{metric}{Metric to plot.}
\item{sender, receiver, cell_type}{Optional filters.}
\item{top_n, top_n_label}{Number of top rows or labels to show.}
\item{sample_ct_list}{Named list of \code{LogicCommCellTypeComm} objects, one per biological sample.}
\item{sample_metadata}{Data frame with rownames matching sample names.}
\item{design}{Model formula for sample-level GLM, for example \code{~ group + batch}.}
\item{coef}{Coefficient name or index to report from the model matrix.}
\item{min_samples, min_total_trials}{Feature-level model filters.}
\item{fdr_method}{P-value adjustment method.}
\item{ct_comm}{Output from \code{summarize_celltype_communication()}.}
\item{file, title, pt_size, verbose, ...}{Output, plotting, and additional arguments.}
}
\details{
Spatial graphs are ordinary sparse adjacency matrices and can be passed to
\code{IdentifyLogicConsensus()} or \code{summarize_celltype_communication()}.
When distance weights are used, set \code{edge_weight_mode = "weighted"} so LCS is
computed as a weighted active-edge fraction.

The dynamics helper scores each pseudotime or time bin independently and returns
bin-level cell-type pair, role, and L-R trajectories. The GLM helper treats
biological samples as the independent statistical units and models active edge
counts as successes with available edge opportunities as trials.
}
\value{
Depending on the function: a sparse graph, a \code{LogicCommCellTypeComm}, a
\code{LogicCommDynamics} object, a GLM result data frame, a ggplot object, a
\code{LogicCommFindings} list, or the path to a markdown report.
}
\examples{
\dontrun{
sp_graph <- build_spatial_graph(coords, mode = "knn", k = 6,
                                distance_weight = "gaussian")
ct_sp <- summarize_celltype_communication(reo, cell_labels = labels,
                                          knn_mat = sp_graph,
                                          edge_weight_mode = "weighted")
plot_spatial_logic(reo, coords, lr_pair = "VEGFA_KDR")

findings <- summarize_communication_findings(ct_sp)
write_communication_report(ct_sp, "LogicComm_report.md")
}
}
