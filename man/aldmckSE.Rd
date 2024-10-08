% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aldmckSE.R
\name{aldmckSE}
\alias{aldmckSE}
\title{Calculate Standard Errors for Aldrich-McKelvey Scaling}
\usage{
aldmckSE(obj, data, ...)
}
\arguments{
\item{obj}{An object containing the Aldrich-McKelvey scaling results. This object should contain respondent parameters (e.g., alpha and beta) and stimuli positions.}

\item{data}{A matrix or dataframe containing the original data used in the Aldrich-McKelvey scaling.}

\item{...}{Additional arguments (currently not used).}
}
\value{
A numeric vector of standard errors for each stimulus.
}
\description{
This function calculates the standard errors for the Aldrich-McKelvey scaling results based on the provided object and data.
}
\examples{
\dontrun{
data(result.france)
data(franceEES2009)
se_values <- aldmckSE(result.france, franceEES2009)
print(se_values)
}
}
