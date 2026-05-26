\name{score_lr_rank_activity}
\alias{score_lr_rank_activity}
\title{Score per-cell ligand-receptor rank activity}
\description{Convert a gene-level REO rank matrix into continuous LR sender and receiver scores.}
\usage{score_lr_rank_activity(rank_mat, lr_db = lr_pairs_human,
  mode = c("both", "sender", "receiver"),
  complex_aggregate = c("min", "mean", "max"), aggregate = TRUE,
  verbose = TRUE)}
\arguments{
\item{rank_mat}{Rank-percentile matrix or \code{LogicCommREOResult}.}
\item{lr_db}{Ligand-receptor database.}
\item{mode}{Whether to score sender, receiver, or both roles.}
\item{complex_aggregate}{How to combine multi-subunit complex ranks.}
\item{aggregate}{Whether to return per-cell aggregate scores.}
\item{verbose}{Whether to print progress messages.}
}
\value{List containing LR-pair x cell sender/receiver rank matrices and optional per-cell aggregate scores.}
\examples{\dontrun{
rank_scores <- score_lr_rank_activity(rank_mat)
}}
