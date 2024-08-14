#' Plot Optimal Classification Cutlines in 2D
#'
#' The `plot_OCcutlines2` function plots cutlines for a 2D Optimal Classification (OC) object. It visualizes the decision boundaries (cutlines) in a two-dimensional space based on the OC analysis results.
#'
#' @param x An object of class `OCobject`. This object should contain the results of an Optimal Classification analysis, including roll call data and norm vectors for plotting cutlines.
#' @param main.title A character string specifying the title of the plot. Defaults to `"Cutting Lines"`.
#' @param d1.title A character string specifying the label for the first dimension (x-axis). Defaults to `"First Dimension"`.
#' @param d2.title A character string specifying the label for the second dimension (y-axis). Defaults to `"Second Dimension"`.
#' @param lines An integer specifying the number of cutlines to plot. If `lines` is a vector, it selects specific cutlines by their row indices. If `lines` is a single integer, it randomly selects that many cutlines to plot. Defaults to `50`.
#' @param dims A numeric vector of length 2 specifying the dimensions to be plotted. Defaults to `c(1, 2)`.
#' @param lwd A numeric value specifying the line width of the cutlines. Defaults to `2`.
#' @param ... Additional graphical parameters passed to the `symbols` function.
#'
#' @return This function does not return a value but produces a object for plot.
#'
#' @details
#' The function checks that the input object `x` is of class `OCobject` and that it contains two dimensions for plotting. If the `lines` parameter is a single integer, it will randomly select that many cutlines to display. If `lines` is a vector, it will use those specific rows from the roll call data to plot the cutlines.#' @importFrom yourpackage add_OCcutline
#' @examples
#' \dontrun{
#' # Assuming `oc_result` is an OCobject with 2D results
#' plot_OCcutlines2(oc_result, lines = 10)
#' }
#'
#' @export
plot_OCcutlines2 <- function (x, main.title = "Cutting Lines", d1.title = "First Dimension", 
                              d2.title = "Second Dimension", lines = 50, dims = c(1, 2), 
                              lwd = 2, ...) 
{
  if (!class(x) == "OCobject") 
    stop("Input is not of class 'OCobject'.")
  if (x$dimensions == 1) 
    stop("All angles in 1D Optimal Classification are 90 degrees.")
  if (length(dims) != 2) 
    stop("'dims' must be an integer vector of length 2.")
  if(length(lines) == 1){
    if (lines < 1) {
      stop("'Lines' must be less than 1.")
    }
  }
  cutlineData <- cbind(x$rollcalls[, paste("normVector", dims[1], 
                                           "D", sep = "")], x$rollcalls[, paste("normVector", dims[2], 
                                                                                "D", sep = "")], x$rollcalls[, "midpoints"])
  cutlineData <- na.omit(cutlineData)
  suppressWarnings(symbols(x = 0, y = 0, circles = 1, inches = FALSE, 
                           asp = 1, main = main.title, xlab = d1.title, ylab = d2.title, 
                           xlim = c(-1, 1), ylim = c(-1, 1), cex.main = 1.2, cex.lab = 1.2, 
                           font.main = 2, lwd = 2, fg = "grey", frame.plot = FALSE, 
                           ...))
  if (length(lines) == 1){
    if(lines < dim(cutlineData)[1]) {
      cutlineData <- cutlineData[sample(1:dim(cutlineData)[1], lines), ]
    }
  } 
  if(length(lines) > 1){
    cutlineData <- cutlineData[lines, ]
  }
  suppressWarnings(apply(cutlineData, 1, add_OCcutline, lwd = lwd))
}  

