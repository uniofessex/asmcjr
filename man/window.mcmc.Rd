% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{window.mcmc}
\alias{window.mcmc}
\title{Subset an MCMC Object by Specifying a Time Window}
\usage{
\method{window}{mcmc}(x, start, end, thin, ...)
}
\arguments{
\item{x}{An object of class \code{mcmc}, representing a single MCMC chain.}

\item{start}{The starting iteration (time) for the subset. If not provided, the starting time of the original chain is used.}

\item{end}{The ending iteration (time) for the subset. If not provided, the ending time of the original chain is used.}

\item{thin}{The thinning interval to be applied. If not provided, the original thinning interval is used. If \code{thin} is not a multiple of the original thinning interval, a warning is issued and the original thinning interval is retained.}

\item{...}{Additional arguments (currently not used).}
}
\value{
A new \code{mcmc} object representing the subsetted chain with the specified start, end, and thinning parameters.
}
\description{
The \code{window.mcmc} function subsets a Markov Chain Monte Carlo (MCMC) object by selecting iterations within a specified time window.
}
\details{
This function allows you to focus on a specific window of iterations within an MCMC chain. It adjusts the start and end points according to the provided values, ensuring they align with the chain's iterations. The thinning interval can also be adjusted, but it must be a multiple of the original thinning interval; otherwise, a warning is issued and the original interval is kept.
}
\examples{
# Assume `mcmc_chain` is an object of class `mcmc`
# Subset the MCMC chain to include only iterations between 500 and 1000
windowed_chain <- window.mcmc(mcmc_chain, start = 500, end = 1000, thin = 2)

}
