\name{plot_shift_lcs_combined}
\alias{plot_shift_lcs_combined}
\title{Combined rank-shift and LCS plot}
\description{Show top differentially active LR pairs together with ligand and receptor rank shifts.}
\usage{plot_shift_lcs_combined(rs_result, lcs_result, top_n = 20,
  fdr_cutoff = 0.1)}
\arguments{
\item{rs_result}{A \code{RankShiftResult}.}
\item{lcs_result}{A \code{LogicCommResult}.}
\item{top_n}{Number of top pairs to display.}
\item{fdr_cutoff}{FDR cutoff for LCS pairs.}
}
\value{A ggplot object.}
