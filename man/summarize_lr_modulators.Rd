\name{summarize_lr_modulators}
\alias{summarize_lr_modulators}
\title{Summarize LR modulator activity}
\description{Summarize CellChatDB-style agonist, antagonist, and co-receptor activity in a LogicComm dataset.}
\usage{summarize_lr_modulators(reo_mat, lr_db, rank_mat = NULL)}
\arguments{
\item{reo_mat}{Binary REO matrix.}
\item{lr_db}{LR database with optional modulator columns.}
\item{rank_mat}{Optional rank-percentile matrix aligned to \code{reo_mat}.}
}
\value{A data frame with modulator presence, activity, and rank-reversal summaries.}
\examples{\dontrun{summarize_lr_modulators(reo, lr_db)}}
