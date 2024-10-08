% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_oc_rollcall.R
\name{plot_oc_rollcall}
\alias{plot_oc_rollcall}
\title{Plot Optimal Classification (OC) Coordinates for a Specific Roll Call}
\usage{
plot_oc_rollcall(
  obj,
  data,
  shapeVar = NULL,
  rcnum,
  rotMat = diag(2),
  dropNV = FALSE,
  onlyErrors = FALSE,
  ptSize = 4
)
}
\arguments{
\item{obj}{An object of class \code{OCobject} containing the results of an Optimal Classification analysis, including legislators' ideal points.}

\item{data}{A \code{rollcall} object containing the voting data corresponding to the OC analysis.}

\item{shapeVar}{An optional variable used to differentiate groups of legislators by shape and color in the plot. If provided, the points will be shaped and colored according to this variable.}

\item{rcnum}{An integer specifying the roll call number to be visualized.}

\item{rotMat}{A 2x2 rotation matrix used to rotate the coordinates. The default is the identity matrix (no rotation).}

\item{dropNV}{A logical value indicating whether to drop non-voters (NA values) from the plot. Default is \code{FALSE}.}

\item{onlyErrors}{A logical value indicating whether to plot only those observations where classification errors occurred. Default is \code{FALSE}.}

\item{ptSize}{A numeric value specifying the size of the points in the plot. Default is \code{4}.}
}
\value{
A \code{ggplot} object displaying the legislators' coordinates for the specified roll call, optionally with different shapes and colors for different groups.
}
\description{
This function creates a scatter plot of legislators' ideal points based on the results of an Optimal Classification (OC) analysis for a specific roll call. The coordinates can be optionally rotated using a rotation matrix and plotted with different shapes and colors for different groups.
}
\examples{
\dontrun{
# Assuming `oc_result` is an OCobject, `rc_data` is a rollcall object,
# and `legislator_data` contains a grouping variable
g <- plot_oc_rollcall(oc_result, data = rc_data, 
shapeVar = legislator_data$group, rcnum = 5, dropNV = TRUE)
print(g)
}
}
