#' Generate a Gray Color Palette
#'
#' This function generates a palette of gray colors.
#'
#' @param n Integer, the number of gray shades to generate.
#' @param lower Numeric, the lower bound for the gray scale (default is 0.3).
#' @param upper Numeric, the upper bound for the gray scale (default is 0.7).
#' @return A vector of gray colors.
#' @importFrom grDevices rgb
#' @examples
#' gray.palette(5)
#' @export

gray.palette <- function(n, lower = .3, upper = .7) {
  s <- seq(lower, upper, length = n)
  rgb(matrix(rep(s, each = 3), ncol = 3, byrow = TRUE))
}
