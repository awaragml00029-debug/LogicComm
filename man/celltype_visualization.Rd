\name{celltype_visualization}
\alias{plot_celltype_heatmap}
\alias{plot_celltype_network}
\alias{plot_celltype_roles}
\alias{plot_lr_bubble_by_celltype}
\alias{plot_pathway_heatmap}
\alias{explain_celltype_interaction}
\alias{plot_celltype_role_heatmap}
\alias{plot_celltype_role_radar}
\alias{plot_lr_activity_balance}
\alias{plot_lr_evidence}
\title{Visualization helpers for cell-type-resolved LogicComm results}
\description{
Visualization and explanation functions for \code{LogicCommCellTypeComm} objects.
}
\usage{
plot_celltype_heatmap(ct_comm,
  metric = c("n_active_lr", "active_lr_event_count", "sum_lcs",
    "mean_lcs_active", "sum_lcs_all", "n_edges", "edge_weight_sum",
    "sum_active_edges", "sum_active_edge_weight"),
  title = NULL, rotate_x = 45)

plot_celltype_network(ct_comm, metric = "sum_lcs", min_value = 0,
  top_n_edges = NULL, show_self = FALSE, title = NULL)

plot_celltype_roles(ct_comm, label = TRUE, title = NULL)

plot_celltype_role_heatmap(ct_comm, role_cols = c("sender_role_score",
  "receiver_role_score", "mediator_role_score", "influencer_role_score"),
  scaled = TRUE, title = NULL)

plot_celltype_role_radar(ct_comm, cell_types = NULL, scaled = TRUE,
  title = NULL)

plot_lr_bubble_by_celltype(ct_comm, sender = NULL, receiver = NULL,
  top_n = 30, active_only = TRUE,
  color_by = c("pathway", "lcs", "n_active_edges",
    "ligand_active_frac_sender", "receptor_active_frac_receiver"),
  title = NULL)

plot_lr_activity_balance(ct_comm, sender = NULL, receiver = NULL,
  active_only = TRUE, top_n = 100, title = NULL)

plot_lr_evidence(ct_comm, sender, receiver, lr_pair, title = NULL)

plot_pathway_heatmap(ct_comm,
  metric = c("sum_lcs", "n_active_lr", "mean_lcs_active"),
  sender = NULL, receiver = NULL, top_n_pathways = 30, title = NULL)

explain_celltype_interaction(ct_comm, sender, receiver, lr_pair,
  reo_mat = NULL, seurat_obj = NULL, lr_db = NULL,
  cell_labels = NULL, pt_size = 0.7)
}
\arguments{
\item{ct_comm}{A \code{LogicCommCellTypeComm} object from \code{summarize_celltype_communication()}.}
\item{metric}{Metric to display. The accepted values depend on the plot.}
\item{title}{Optional plot title.}
\item{rotate_x}{Rotation angle for x-axis labels.}
\item{min_value}{Minimum edge value for network plotting.}
\item{top_n_edges}{Maximum number of strongest network edges to show.}
\item{show_self}{Whether to draw same-cell-type network edges.}
\item{label}{Whether to label cell types in the role plot.}
\item{role_cols}{Role score columns to display.}
\item{scaled}{Whether to z-score or rescale role scores for visualization.}
\item{cell_types}{Optional cell types for the radar plot.}
\item{sender}{Optional sender cell type filter, or the selected sender in explanation functions.}
\item{receiver}{Optional receiver cell type filter, or the selected receiver in explanation functions.}
\item{top_n}{Number of strongest L-R rows to show.}
\item{active_only}{Whether to show only active cell-type L-R rows.}
\item{color_by}{Variable used to color L-R bubbles.}
\item{top_n_pathways}{Number of top pathways to display.}
\item{lr_pair}{Selected L-R pair for explanation.}
\item{reo_mat}{Optional REO matrix used for UMAP logic plotting.}
\item{seurat_obj}{Optional Seurat object with UMAP coordinates.}
\item{lr_db}{Optional LR database; defaults to the database stored in \code{ct_comm}.}
\item{cell_labels}{Optional named cell labels; defaults to labels stored in \code{ct_comm}.}
\item{pt_size}{UMAP point size.}
}
\value{
Plot functions return ggplot2 objects. \code{explain_celltype_interaction()} returns a list with an evidence table and, when possible, a ggplot object.
}
\examples{
\dontrun{
plot_celltype_heatmap(ct_comm, metric = "n_active_lr")
plot_celltype_roles(ct_comm)
plot_celltype_role_heatmap(ct_comm)
plot_lr_activity_balance(ct_comm, sender = "DC", receiver = "CD8 T")
plot_lr_evidence(ct_comm, sender = "DC", receiver = "CD8 T", lr_pair = "CD80_CD28")
}
}
