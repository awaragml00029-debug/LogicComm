\name{print-methods}
\alias{print.LCSVector}
\alias{print.LogicCommResult}
\alias{print.LogicCommREOResult}
\alias{print.LogicGateResult}
\alias{print.LogicCommMulti}
\alias{print.LogicCommCellTypeComm}
\alias{print.RankShiftResult}
\title{Print methods for LogicComm objects}
\description{Compact print methods for LogicComm result classes.}
\usage{\method{print}{LCSVector}(x, ...)
\method{print}{LogicCommResult}(x, n = 10, ...)
\method{print}{LogicCommREOResult}(x, ...)
\method{print}{LogicGateResult}(x, ...)
\method{print}{LogicCommMulti}(x, ...)
\method{print}{LogicCommCellTypeComm}(x, ...)
\method{print}{RankShiftResult}(x, n = 10, ...)}
\arguments{
\item{x}{LogicComm object to print.}
\item{n}{Number of rows to display when applicable.}
\item{...}{Ignored.}
}
\value{Invisibly returns \code{x}.}
