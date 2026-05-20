\name{score_communication_specificity}
\alias{score_communication_specificity}
\alias{interpret_celltype_roles}
\alias{rank_communication_axes}
\alias{diagnose_celltype_communication}
\alias{diagnose_permutation_resolution}
\alias{plot_pathway_dominance}
\alias{plot_celltype_role_dotplot}
\alias{plot_specificity_stability}
\alias{plot_celltype_network_publication}
\alias{plot_pathway_network}
\title{Specificity, diagnostics, and publication-oriented interpretation}
\usage{
score_communication_specificity(
  ct_comm,
  active_only = TRUE,
  broad_pair_specificity_cutoff = 0.25,
  ubiquitous_pair_fraction = 0.25,
  identity_pathway_patterns = c("^MHC", "ANTIGEN", "BCR", "TCR"),
  identity_lr_patterns = c("^HLA", "_CD4$", "_CD8", "B2M"),
  verbose = TRUE
)

interpret_celltype_roles(ct_comm, cell_type = NULL)

rank_communication_axes(ct_comm, null_pair = NULL, sens = NULL, top_n = NULL)

diagnose_celltype_communication(
  ct_comm,
  null_pair = NULL,
  sens = NULL,
  min_cells_per_type = 30,
  dominant_pathway_fraction_warning = 0.45,
  min_n_perm_publication = 1000
)

diagnose_permutation_resolution(null_result, min_n_perm_publication = 1000)

plot_pathway_dominance(ct_comm, metric = "sum_lcs", top_n = 15, title = NULL)

plot_celltype_role_dotplot(ct_comm, scaled = TRUE, title = NULL)

plot_specificity_stability(sens, ct_comm = NULL, specificity = NULL,
  top_n_label = 15, title = NULL)

plot_celltype_network_publication(ct_comm, metric = "sum_lcs",
  top_n_edges = 20, min_value = 0,
  edge_label = c("top_pathway", "none", "n_active_lr"),
  show_self = FALSE, title = NULL)

plot_pathway_network(ct_comm, pathway, metric = "sum_lcs", top_n_edges = 20,
  title = NULL)
}
\arguments{
\item{ct_comm}{A \code{LogicCommCellTypeComm} object from \code{summarize_celltype_communication()}.}
\item{active_only}{Whether specificity should be calculated from active L-R rows only.}
\item{broad_pair_specificity_cutoff}{Pair-specificity threshold below which an L-R pair is considered broad.}
\item{ubiquitous_pair_fraction}{Fraction of possible sender-receiver pairs above which an L-R pair is considered ubiquitous.}
\item{identity_pathway_patterns, identity_lr_patterns}{Regular expressions for identity-associated pathways or L-R names, for example MHC/HLA axes in PBMC.}
\item{verbose}{Whether to print progress messages.}
\item{cell_type}{Optional cell-type subset for role interpretation.}
\item{null_pair, null_result}{Permutation-null result from \code{permute_celltype_communication()}.}
\item{sens}{Sensitivity-analysis result from \code{sensitivity_REO_threshold()}.}
\item{top_n, top_n_label, top_n_edges, top_n_pathways}{Number of rows, labels, or edges to display.}
\item{min_cells_per_type}{Small-cell-type warning threshold.}
\item{dominant_pathway_fraction_warning}{Pathway-dominance warning threshold.}
\item{min_n_perm_publication}{Recommended minimum number of permutations for publication-level empirical p-values.}
\item{metric}{Metric to plot or rank.}
\item{title}{Optional plot title.}
\item{scaled}{Whether to z-score role scores before plotting.}
\item{specificity}{Optional specificity summary table.}
\item{min_value}{Minimum edge metric for network plotting.}
\item{edge_label}{Edge labels in the publication network.}
\item{show_self}{Whether to show self/self cell-type edges.}
\item{pathway}{Pathway name to display in a pathway-filtered network.}
}
\description{
These helpers convert LogicComm cell-type communication results into
publication-oriented interpretation layers. They distinguish pair-specific
communication axes from broad or identity-associated programs, split cell-type
role interpretation into role separation and communication evidence, diagnose
permutation p-value resolution and degenerate nulls, and provide figure-ready
visualizations for roles, pathway dominance, specificity-stability landscapes,
and top-edge networks.
}
\details{
\code{score_communication_specificity()} computes entropy-based sender,
receiver, and sender-receiver-pair specificity for each active L-R pair. A stable
but low-specificity MHC/HLA axis in PBMC should usually be described as a broad
immune-recognition or identity-associated program unless it also has high
cell-type-pair specificity or receiver-response evidence.

\code{diagnose_permutation_resolution()} reports the minimum empirical p-value
implied by the number of permutations. When \code{null_sd == 0} and the observed
score is larger than the null mean, LogicComm treats the z-score as \code{Inf}
and flags a degenerate positive null. Biologically this means the shuffled labels
never produced comparable communication under the current graph and REO states;
this is a structural-null signal that should be interpreted together with edge
support, cell abundance, specificity, and independent replication.
}
\value{
\code{score_communication_specificity()} returns an augmented
\code{LogicCommCellTypeComm} object. \code{interpret_celltype_roles()},
\code{rank_communication_axes()}, \code{diagnose_permutation_resolution()}, and
\code{diagnose_celltype_communication()} return data frames or diagnostic lists.
Plot functions return \code{ggplot2} objects.
}
\examples{
\dontrun{
ct_comm <- summarize_celltype_communication(reo, cell_labels, knn_mat = snn)
ct_comm <- score_communication_specificity(ct_comm)
interpret_celltype_roles(ct_comm)
plot_celltype_role_dotplot(ct_comm)
plot_pathway_dominance(ct_comm)

null_pair <- permute_celltype_communication(ct_comm, reo_mat = reo, knn_mat = snn,
                                            n_perm = 1000)
diagnose_permutation_resolution(null_pair)

diag <- diagnose_celltype_communication(ct_comm, null_pair = null_pair)
diag
}
}
