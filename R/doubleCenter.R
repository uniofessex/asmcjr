#' Double Centering of a Matrix
#'
#' This function performs double centering on a matrix. Double centering is a common preprocessing step in multivariate data analysis, where the matrix is adjusted so that both the row and column means are zero.
#'
#' @param x A numeric matrix that you want to double center. The matrix should be of size \eqn{p \times n}, where \eqn{p} is the number of rows and \eqn{n} is the number of columns.
#' 
#' @return A numeric matrix of the same dimensions as `x`, where the rows and columns have been centered so that their means are zero.
#'
#' @details The function computes the double-centered matrix by subtracting the row and column means and adjusting for the grand mean. The resulting matrix will have row and column means equal to zero.
#'
#' @examples
#' \dontrun{
#' # Load the example dataset
#' data(nations)
#' 
#' # Print the nations dataset
#' print(nations)
#' 
#' # Calculate squared differences
#' d <- (9 - nations)^2
#' 
#' # Perform double centering on the squared differences matrix
#' D <- doubleCenter(d)
#' 
#' # Print the double-centered matrix
#' print(D)
#' }
#' @export
doubleCenter <- function(x){
  p <- dim(x)[1]
  n <- dim(x)[2]
  -(x - matrix(apply(x,1,mean),nrow=p,ncol=n) -
      t(matrix(apply(x,2,mean),nrow=n,ncol=p)) + mean(x)) / 2
}
