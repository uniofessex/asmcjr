% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{france4}
\alias{france4}
\title{Roll Call Data from the French Fourth Republic}
\format{
A data frame with the following variables:
\describe{
\item{CASEID}{Unique identifier for each deputy's party affiliation. Changes if the deputy switches parties.}
\item{MID}{Unique ID for each deputy that remains constant even if the deputy switches party.}
\item{NAME}{Name of the deputy (legislator).}
\item{PAR}{Party affiliation of the deputy.}
\item{PARSEQ}{Sequence number of the party affiliation for deputies who switched parties.}
\item{V1001, V1002, ..., V3172}{Results of the roll call votes, with each variable representing the result of a specific roll call vote.}
}
}
\source{
Data from Rosenthal and Voeten (2004). The original dataset and documentation can be found online at
\url{http://www9.georgetown.edu/faculty/ev42/france.htm}.
}
\usage{
data(france4)
}
\description{
The \code{france4} dataset contains roll call data from the French Fourth Republic, as analyzed by Rosenthal and Voeten (2004).
This dataset was used to estimate a party-switcher model where a separate ideal point is estimated each time a legislator
changes party affiliation. The dataset includes legislator-specific variables and roll call votes. Rosenthal and Voeten (2004)
found that the latent ideological space remained stable over the course of the French Fourth Republic, so the roll call data
is not segmented by legislative session.
}
\details{
The dataset includes 1,416 separate ideal points, as each party-switching deputy has a separate entry for each party affiliation.
The first five columns of the dataset are legislator-specific variables: NAME (deputy name), MID (constant deputy ID), CASEID
(unique ID for each party affiliation), PAR (party affiliation), and PARSEQ (sequence number of party affiliations for party-switching deputies).
The remaining columns represent the roll call votes.

The dataset and its extensive documentation were originally made available online by Rosenthal and Voeten (2004).
For more information and access to the original data, visit \url{http://www9.georgetown.edu/faculty/ev42/france.htm}.
}
\examples{
\dontrun{
data(france4)
}

}
\keyword{datasets}
