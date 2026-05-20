\name{plot_lcs_heatmap}
\alias{plot_lcs_heatmap}
\title{Heatmap of LCS across samples}
\description{Plot LCS values for top ligand-receptor pairs across samples.}
\usage{plot_lcs_heatmap(result, group_info, top_n = 30, fdr_cutoff = 0.1,
  scale_rows = TRUE, color_palette = NULL, ...)}
\arguments{
\item{result}{A \code{LogicCommResult}.}
\item{group_info}{Named group labels.}
\item{top_n}{Number of top pairs to plot.}
\item{fdr_cutoff}{FDR cutoff used to select pairs.}
\item{scale_rows}{Whether to Z-score rows.}
\item{color_palette}{Optional color palette.}
\item{...}{Additional arguments passed to \code{pheatmap::pheatmap()}.}
}
\value{Invisibly returns the pheatmap object.}
