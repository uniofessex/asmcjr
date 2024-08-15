
#' Plot OC Coordinates with Rotation
#'
#' This function plots the OC (Optimal Classification) coordinates after applying an optional rotation matrix.
#'
#' @param obj A list object containing the OC results. Typically, this comes from a function like `oc`.
#' @param shapeVar Optional. A variable to define the shape and color of points in the plot, typically a factor.
#' @param dropNV Logical. If `TRUE`, rows with missing values in the `shapeVar` will be dropped from the plot.
#' @param ptSize Numeric. The size of the points in the plot.
#' @param rotMat A 2x2 rotation matrix. Defaults to the identity matrix (`diag(2)`), which means no rotation.
#' 
#' @return A ggplot object representing the OC coordinates plot.
#' 
#' @examples
#' # Example usage:
#' # Assume `rc` is an object containing legislative data
#' result <- oc(rc, dims=2, minvotes=20, lop=0.025, polarity=c(2,2), verbose=FALSE)
#'
#' # Define a 45-degree rotation matrix
#' deg2rad <- function(x) x * pi / 180
#' rad45 <- deg2rad(45)
#' A <- matrix(c(cos(rad45), -sin(rad45), sin(rad45), cos(rad45)), nrow=2, ncol=2, byrow=TRUE)
#'
#' # Plot OC coordinates without rotation
#' plot1 <- plot_oc_coords(result)
#' print(plot1)
#'
#' # Plot OC coordinates with a 45-degree rotation
#' plot2 <- plot_oc_coords(result, rotMat=A)
#'
#' @import ggplot2
#' @export
plot_oc_coords <- function(obj, shapeVar=NULL, dropNV=FALSE, ptSize=4, rotMat=diag(2)){
  A <- rotMat
  rot.oc <- cbind(obj$legislators$coord1D, obj$legislators$coord2D)
  for (i in 1:nrow(rot.oc)){
    rot.oc[i,] <- rot.oc[i,] %*% A
  }
  oc.dat <- data.frame(
    X1 = rot.oc[,1], 
    X2 = rot.oc[,2])
  if(!is.null(shapeVar)){
    oc.dat$group <- shapeVar
  }
  if(dropNV){
    oc.dat <- na.omit(oc.dat)
  }
  if(is.null(shapeVar)){
    g <- ggplot(oc.dat, aes_string(x="X1", y="X2")) + geom_point(size=ptSize) 
  }
  if(!is.null(shapeVar)){
    g <- ggplot(oc.dat, aes_string(x="X1", y="X2", colour="group", shape="group")) + geom_point(size=ptSize) 
  }
  g
}
