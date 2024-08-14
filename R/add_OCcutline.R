#' Add Cutline to Plot in Optimal Classification (OC) Analysis
#'
#' The `add_OCcutline` function calculates and adds a cutline to a plot based on the results of an Optimal Classification (OC) analysis. The cutline is determined by the slope and intercept derived from the provided `cutData`.
#'
#' @param cutData A numeric vector containing the parameters needed to compute the cutline. This typically includes coefficients and midpoint values from the OC analysis.
#' @param lwd A numeric value specifying the line width of the cutline. The default value is `2`.
#'
#' @details
#' The function calculates the slope and intercept of the cutline based on the values provided in `cutData`. If the slope is not defined (i.e., it is `NA`), the function will attempt to plot a vertical line at the specified intercept. If the slope is defined, the function calculates the points where the line intersects the unit circle and plots the resulting line segment.
#'
#' This function is intended for internal use within the package and is primarily designed to assist with visualizing OC analysis results.
#' @importFrom graphics lines
#' @keywords internal
#' @export
add_OCcutline <- function(cutData,lwd=2) {
  
  slope <- -cutData[1]/cutData[2]
  if (is.na(slope)) {
    x <- c(cutData[1]*cutData[3],cutData[1]*cutData[3])
    y <- c(sqrt(1-(cutData[1]*cutData[3])^2),-sqrt(1-(cutData[1]*cutData[3])^2))
    slope <- NA
    intercept <- NA
  }
  else {
    intercept <- -slope*cutData[1]*cutData[3]+cutData[2]*cutData[3]
    x <- c( (-slope*intercept + sqrt( (slope*intercept)^2 -
                                        (1+slope*slope)*(intercept*intercept-1)))/(1+slope*slope),
            (-slope*intercept - sqrt( (slope*intercept)^2 - 
                                        (1+slope*slope)*(intercept*intercept-1)))/(1+slope*slope) )
    if (is.na(x[1])) {
      warning("Couldn't solve for points on the unit circle!\n")
      x<-NA
      y<-NA
      slope<-NA
      intercept<-NA  
    }             
    else {
      y <- intercept + slope*x
      y[y < -1] <- -sqrt(1-x[y<1]^2)
      y[y >  1] <-  sqrt(1-x[y>1]^2)
    }
  }
  lines(x,y,lwd=lwd)
}
