\name{differential_celltype_plots}
\alias{plot_differential_celltype_heatmap}
\alias{plot_differential_celltype_volcano}
\title{Differential cell-type communication plots}
\description{
Plot differential communication results generated from cell-type-level vectors
and \code{CompareLogicGroups()}.
}
\usage{
plot_differential_celltype_heatmap(result,
  metric = c("asymmetry", "log2fc_lcs", "case_mean_lcs", "ctrl_mean_lcs"),
  agg = c("sum", "mean", "max_abs"), title = NULL)

plot_differential_celltype_volcano(result,
  x = c("log2fc_lcs", "asymmetry"),
  p_col = c("fdr_fisher", "fdr_wilcox", "p_fisher", "p_wilcox"),
  color_by = c("pathway", "feature_type"), top_n_label = 15,
  title = NULL)
}
\arguments{
\item{result}{A \code{LogicCommResult}, usually produced from \code{celltype_comm_to_lcs()} vectors.}
\item{metric}{Comparison metric for the heatmap.}
\item{agg}{How to aggregate multiple features per sender-receiver pair.}
\item{title}{Optional plot title.}
\item{x}{Volcano x-axis column.}
\item{p_col}{P-value or FDR column for the y-axis.}
\item{color_by}{Column used for colors.}
\item{top_n_label}{Number of strongest points to label.}
}
\value{ggplot2 objects.}
\examples{
\dontrun{
vecs <- lapply(ct_comm_list, celltype_comm_to_lcs, level = "celltype_lr")
res <- CompareLogicGroups(vecs, group_info)
plot_differential_celltype_heatmap(res)
plot_differential_celltype_volcano(res)
}
}
