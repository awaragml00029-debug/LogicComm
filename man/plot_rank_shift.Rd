\name{plot_rank_shift}
\alias{plot_rank_shift}
\title{Rank-shift volcano plot}
\description{Plot gene rank shifts and FDR values.}
\usage{plot_rank_shift(rs_result, top_n_label = 20, fdr_cutoff = 0.05,
  min_shift_label = 0.1, color_by_role = TRUE, highlight_genes = NULL,
  point_size = 2, label_size = 3, title = NULL)}
\arguments{
\item{rs_result}{A \code{RankShiftResult}.}
\item{top_n_label}{Maximum number of genes to label.}
\item{fdr_cutoff}{FDR cutoff.}
\item{min_shift_label}{Minimum absolute shift for labeling.}
\item{color_by_role}{Whether to color by significance direction.}
\item{highlight_genes}{Optional genes to label.}
\item{point_size}{Point size.}
\item{label_size}{Label size.}
\item{title}{Optional plot title.}
}
\value{A ggplot object.}
