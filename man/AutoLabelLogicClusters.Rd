\name{AutoLabelLogicClusters}
\alias{AutoLabelLogicClusters}
\title{Label clusters by ligand-receptor logic activity}
\description{Identify clusters enriched for sender or receiver logic for selected ligand-receptor pairs.}
\usage{AutoLabelLogicClusters(reo_mat, cell_labels, lr_db = lr_pairs_human,
  lr_pairs = NULL, min_sender_frac = 0.1, min_receiver_frac = 0.1,
  verbose = TRUE)}
\arguments{
\item{reo_mat}{Binary REO matrix.}
\item{cell_labels}{Named cluster/cell labels or a Seurat object.}
\item{lr_db}{Ligand-receptor database.}
\item{lr_pairs}{Optional subset of LR pair names.}
\item{min_sender_frac}{Minimum sender fraction for reporting.}
\item{min_receiver_frac}{Minimum receiver fraction for reporting.}
\item{verbose}{Whether to print progress messages.}
}
\value{A data frame summarizing sender and receiver cluster enrichment.}
\examples{\dontrun{labels <- AutoLabelLogicClusters(reo, cell_labels)}}
