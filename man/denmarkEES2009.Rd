% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{denmarkEES2009}
\alias{denmarkEES2009}
\title{Danish Module of the 2009 European Election Study (EES)}
\format{
A matrix with 1,000 rows and 8 columns:
\describe{
\item{V1}{The propensity to vote for Party 1.}
\item{V2}{The propensity to vote for Party 2.}
\item{V3}{The propensity to vote for Party 3.}
\item{V4}{The propensity to vote for Party 4.}
\item{V5}{The propensity to vote for Party 5.}
\item{V6}{The propensity to vote for Party 6.}
\item{V7}{The propensity to vote for Party 7.}
\item{V8}{The propensity to vote for Party 8.}
}
}
\usage{
data(denmarkEES2009)
}
\description{
The \code{denmarkEES2009} dataset contains data from the Danish module of the 2009 European Election Study (EES).
This dataset includes responses from 1,000 Danish participants who rated their propensity to vote for each of
eight political parties on a 0–10 point scale. A score of 0 denotes "not at all possible," while a score of 10 denotes "very probable."
}
\details{
The dataset can be used to demonstrate the \code{smacofRect()} function.

There are only 61 missing ratings in the 1000×8 matrix.

The dataset is particularly useful for demonstrating the use of the \code{smacofRect()} function in multidimensional scaling analysis.
}
\examples{
data(denmarkEES2009) 
}
\seealso{
\code{\link[smacof]{smacofRect}} for more details on the \code{smacofRect} function.
}
\keyword{datasets}
