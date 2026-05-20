\name{receiver_response}
\alias{score_receiver_response}
\alias{add_receiver_response_score}
\title{Add receiver downstream response evidence to LogicComm communication}
\description{
Scores downstream response-gene logic in receiver cell types for each
celltype-resolved L-R event.
}
\usage{
score_receiver_response(ct_comm, reo_mat, response_db,
  response_col = NULL, response_mode = c("any", "all"))

add_receiver_response_score(ct_comm, reo_mat, response_db,
  response_col = NULL, response_mode = c("any", "all"))
}
\arguments{
\item{ct_comm}{A \code{LogicCommCellTypeComm} object.}
\item{reo_mat}{Binary REO matrix or \code{LogicCommREOResult}.}
\item{response_db}{Data frame containing \code{lr_pair} and a response gene column.}
\item{response_col}{Response gene column name. If \code{NULL}, tries \code{response_genes}, \code{target_genes}, then \code{genes}.}
\item{response_mode}{\code{"any"} requires at least one response gene active; \code{"all"} requires all available response genes active.}
}
\details{
The returned response columns include \code{response_active_frac_receiver} and
\code{response_integrated_score = lcs * response_active_frac_receiver}. This
adds evidence that a ligand-receptor potential is accompanied by a receiver-side
functional state.
}
\value{
\code{score_receiver_response()} returns a data frame matching \code{ct_comm$lr_table}.
\code{add_receiver_response_score()} returns an updated \code{LogicCommCellTypeComm} object.
}
