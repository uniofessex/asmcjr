% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doubleCenterRect.R
\name{doubleCenterRect}
\alias{doubleCenterRect}
\title{Double Centering of a Rectangular Matrix}
\usage{
doubleCenterRect(x)
}
\arguments{
\item{x}{A numeric rectangular matrix that you want to double center. The matrix can have different numbers of rows and columns.}
}
\value{
A numeric matrix of the same dimensions as \code{x}, with both rows and columns centered so that their means are zero.
}
\description{
This function performs double centering on a rectangular matrix. Double centering is a technique used to adjust a matrix so that both the row and column means are zero, even if the matrix is not square.
}
\details{
The function calculates the double-centered matrix by subtracting the row means, column means, and adjusting for the grand mean. The result is a matrix where both row and column means are zero.
}
\examples{
\dontrun{
# Load the example dataset
data("interest1981")

# Prepare the input matrix from columns 9 to 38
input <- as.matrix(interest1981[, 9:38])

# Set a cutoff for the minimum number of non-NA values per row
cutoff <- 5
input <- input[rowSums(!is.na(input)) >= cutoff, ]

# Transform the input data
input <- (100 - input) / 50

# Square the transformed input values
input2 <- input * input

# Replace NA values with the square of the mean of non-NA values
input2[is.na(input)] <- (mean(input, na.rm = TRUE))^2

# Perform double centering on the rectangular matrix
inputDC <- doubleCenterRect(input2)

# Print the double-centered matrix
print(inputDC)
}

}
