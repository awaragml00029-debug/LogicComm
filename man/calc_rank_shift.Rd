\name{calc_rank_shift}
\alias{calc_rank_shift}
\title{Compute gene expression rank shifts between groups}
\description{Detect genes whose within-cell expression rank changes between Case and Control groups.}
\usage{calc_rank_shift(sample_list, group_info, genes = NULL,
  lr_db = lr_pairs_human, case_label = "Case", ctrl_label = "Ctrl",
  min_detection_frac = 0.05, layer = "counts", verbose = TRUE)}
\arguments{
\item{sample_list}{Named list of expression matrices or Seurat objects.}
\item{group_info}{Group labels for samples.}
\item{genes}{Optional genes to analyze.}
\item{lr_db}{LR database used when \code{genes = NULL}.}
\item{case_label}{Case group label.}
\item{ctrl_label}{Control group label.}
\item{min_detection_frac}{Minimum cell detection fraction per sample.}
\item{layer}{Seurat layer or slot.}
\item{verbose}{Whether to print progress messages.}
}
\value{A data frame of class \code{RankShiftResult}.}
\examples{\dontrun{rs <- calc_rank_shift(sample_list, group_info)}}
