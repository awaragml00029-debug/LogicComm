\name{plot_lcs_bubble}
\alias{plot_lcs_bubble}
\title{Bubble plot of Logic Consensus Scores}
\description{Visualize Case versus Control ligand-receptor activity frequencies.}
\usage{plot_lcs_bubble(result, top_n = 20, fdr_cutoff = 0.05,
  min_asymmetry = 0.2, color_by = "log2fc", pathway_filter = NULL,
  point_alpha = 0.75, label_size = 3, title = NULL)}
\arguments{
\item{result}{A \code{LogicCommResult}.}
\item{top_n}{Maximum number of labels.}
\item{fdr_cutoff}{FDR cutoff for labels.}
\item{min_asymmetry}{Minimum asymmetry for labels.}
\item{color_by}{Color variable, \code{"log2fc"} or \code{"asymmetry"}.}
\item{pathway_filter}{Optional pathways to highlight.}
\item{point_alpha}{Point transparency.}
\item{label_size}{Label size.}
\item{title}{Optional plot title.}
}
\value{A ggplot object.}
