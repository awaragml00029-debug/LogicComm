\name{filter_lcs}
\alias{filter_lcs}
\title{Filter LogicComm comparison results}
\description{Filter a \code{LogicCommResult} by FDR, asymmetry, and direction.}
\usage{filter_lcs(result, min_asymmetry = 0.3, max_fdr = 0.05,
  direction = "up")}
\arguments{
\item{result}{A \code{LogicCommResult}.}
\item{min_asymmetry}{Minimum absolute asymmetry.}
\item{max_fdr}{Maximum Fisher FDR.}
\item{direction}{One of \code{"up"}, \code{"down"}, or \code{"both"}.}
}
\value{A filtered \code{LogicCommResult}.}
