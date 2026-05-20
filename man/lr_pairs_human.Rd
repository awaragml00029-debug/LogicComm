\name{lr_pairs_human}
\alias{lr_pairs_human}
\docType{data}
\title{Built-in human ligand-receptor pair database}
\description{A compact human ligand-receptor seed database in LogicComm format.}
\format{A data frame with columns \code{lr_pair}, \code{ligand}, \code{receptor}, list columns \code{ligand_genes} and \code{receptor_genes}, and metadata columns such as \code{pathway} and \code{annotation}.}
\usage{data(lr_pairs_human)}
\details{The seed set was compiled from public ligand-receptor resources and normalized into explicit ligand and receptor gene list columns for REO/LCS scoring.}
\examples{data(lr_pairs_human)
head(lr_pairs_human[, c("lr_pair", "ligand", "receptor")])
}
\keyword{datasets}
