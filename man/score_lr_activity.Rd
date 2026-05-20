\name{score_lr_activity}
\alias{score_lr_activity}
\alias{rank_comm_cells}
\title{Score per-cell ligand-receptor activity}
\description{Compute sender, receiver, and communication activity scores for cells.}
\usage{score_lr_activity(reo_mat, lr_db = lr_pairs_human, mode = "both",
  aggregate = TRUE, verbose = TRUE)
rank_comm_cells(activity_scores, n = 50, cell_labels = NULL)}
\arguments{
\item{reo_mat}{Binary REO matrix.}
\item{lr_db}{Ligand-receptor database.}
\item{mode}{One of \code{"sender"}, \code{"receiver"}, or \code{"both"}.}
\item{aggregate}{Whether to compute aggregate per-cell scores.}
\item{verbose}{Whether to print progress messages.}
\item{activity_scores}{Output from \code{score_lr_activity()}.}
\item{n}{Number of top cells returned.}
\item{cell_labels}{Optional labels to append to the ranked cells table.}
}
\value{\code{score_lr_activity()} returns a list of sparse matrices and score vectors; \code{rank_comm_cells()} returns a data frame.}
\examples{\dontrun{scores <- score_lr_activity(reo); rank_comm_cells(scores)}}
