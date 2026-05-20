\name{CompareLogicGroups}
\alias{CompareLogicGroups}
\title{Compare LCS values between two sample groups}
\description{Compare per-sample LCS vectors between Case and Control groups.}
\usage{CompareLogicGroups(lcs_list, group_info, case_label = "Case",
  ctrl_label = "Ctrl", lcs_threshold = 0.01, min_case_freq = 0,
  lr_db = NULL, min_samples_per_group = 1, verbose = TRUE)}
\arguments{
\item{lcs_list}{Named list of per-sample LCS vectors.}
\item{group_info}{Named group labels for samples.}
\item{case_label}{Case group label.}
\item{ctrl_label}{Control group label.}
\item{lcs_threshold}{Threshold for calling a sample positive.}
\item{min_case_freq}{Minimum Case frequency retained in output.}
\item{lr_db}{Optional LR database used to append metadata.}
\item{min_samples_per_group}{Minimum non-missing samples required per group for each pair.}
\item{verbose}{Whether to print progress messages.}
}
\details{Missing LCS values are treated as unavailable, not as negative calls.}
\value{A data frame of class \code{LogicCommResult}.}
\examples{\dontrun{res <- CompareLogicGroups(lcs_list, group_info)}}
