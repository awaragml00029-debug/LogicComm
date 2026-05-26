\name{CompareRankLogicGroups}
\alias{CompareRankLogicGroups}
\title{Compare rank-aware logic evidence between groups}
\description{Compare continuous rank-aware LR evidence across biological samples.}
\usage{CompareRankLogicGroups(rank_result_list, group_info,
  score_col = "specificity_weighted_rank_lcs", case_label = "Case",
  ctrl_label = "Ctrl", rank_score_threshold = 0.05,
  min_samples_per_group = 2, lr_db = NULL, verbose = TRUE)}
\arguments{
\item{rank_result_list}{Named list of results from \code{\link{IdentifyRankLogicConsensus}}.}
\item{group_info}{Named vector mapping samples to group labels.}
\item{score_col}{Rank evidence column to compare.}
\item{case_label}{Case group label.}
\item{ctrl_label}{Control group label.}
\item{rank_score_threshold}{Threshold for rank-positive sample calls.}
\item{min_samples_per_group}{Minimum samples required per group.}
\item{lr_db}{Optional LR database for metadata.}
\item{verbose}{Whether to print progress messages.}
}
\value{A data frame sorted by case-control rank evidence difference.}
\examples{\dontrun{
cmp <- CompareRankLogicGroups(rank_result_list, group_info)
}}
