% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneExpressionSimulation.R
\name{simulateMicroArray}
\alias{simulateMicroArray}
\title{Bulk MicroArray Error Model}
\usage{
simulateMicroArray(meanExp, params)
}
\arguments{
\item{meanExp}{matrix of mean expression values}

\item{params}{GenExpressionParams object}
}
\value{
matrix with error model incorporated
}
\description{
adds a normally distributed error across the entire expression
matrix
}
\keyword{internal}
