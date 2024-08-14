#' Plot W-NOMINATE Coordinates with Optional Confidence Intervals
#'
#' This function creates a scatter plot of legislators' ideal points based on W-NOMINATE analysis results. The plot can optionally include confidence ellipses around the points, representing uncertainty in the coordinates.
#'
#' @param obj An object of class `wnominate` containing the results of a W-NOMINATE analysis, including legislators' ideal points and associated uncertainties.
#' @param shapeVar An optional variable used to differentiate groups of legislators by shape and color in the plot. If provided, the points will be shaped and colored according to this variable.
#' @param dropNV A logical value indicating whether to drop non-voters (NA values) from the plot. Default is `FALSE`.
#' @param ptSize A numeric value specifying the size of the points in the plot. Default is `4`.
#' @param ci A logical value indicating whether to include confidence ellipses around the points. Default is `FALSE`.
#' @param level A numeric value specifying the confidence level for the ellipses. Default is `0.95`.
#' 
#' @return A `ggplot` object displaying the legislators' coordinates, optionally with confidence ellipses.
#' 
#' @examples
#' \dontrun{
#' # Assuming `result` is a wnominate object and `legislator_data` contains a grouping variable
#' g <- plot_wnom_coords(result, dropNV = TRUE, ci = TRUE)
#' print(g)
#' }
#' @import ggplot2
#' @import ellipse
#' @export
plot_wnom_coords <- function(obj, shapeVar=NULL, dropNV=FALSE, ptSize=4, ci=FALSE, level=.95){
  weight <-  obj$weights[2]/obj$weights[1] 
  wnom.dat <- data.frame(
    X1 = obj$legislators$coord1D, 
    X2 = obj$legislators$coord2D*weight)
  if(!is.null(shapeVar)){
    wnom.dat$group <- shapeVar
  }
  if(dropNV){
    wnom.dat <- na.omit(wnom.dat)
  }
  if(!ci){
    if(is.null(shapeVar)){
      g <- ggplot(wnom.dat) + geom_point(aes_string(x="X1", y="X2"), size=ptSize) 
    }
    if(!is.null(shapeVar)){
      g <- ggplot(wnom.dat) + geom_point(aes_string(x="X1", y="X2", colour="group", shape="group"), size=ptSize) 
    }
  }
  if(ci){
    wnom.dat <- data.frame(
      X1 = obj$legislators$coord1D, 
      X2 = obj$legislators$coord2D*weight, 
      se1 = obj$legislators$se1D, 
      se2 = obj$legislators$se2D, 
      r = obj$legislators$corr.1)
    if(!is.null(shapeVar)){
      wnom.dat$group <- shapeVar
    }
    if(dropNV){
      wnom.dat <- na.omit(wnom.dat)
    }
    elldat <- NULL
    for(i in 1:nrow(wnom.dat)){
      R <- diag(2)
      R[2,1] <- R[1,2] <- wnom.dat$r[i]
      scl <- c(wnom.dat$se1[i], wnom.dat$se2[i])
      cent <- c(wnom.dat$X1[i], wnom.dat$X2[i])
      ell <- ellipse(R, scale=scl, centre=cent, level=level)
      tmp <- data.frame(x = ell[,1], y = ell[,2], stim = i)
      elldat <- rbind(elldat, tmp)
    }
    if(is.null(shapeVar)){
      g <- ggplot(wnom.dat) + geom_path(data=elldat, aes_string(x="x", y="y", group="stim"), col="gray50", alpha=.25) + 
        geom_point(aes_string(x="X1", y="X2"), size=ptSize) 
    }
    if(!is.null(shapeVar)){
      g <- ggplot(wnom.dat) + geom_path(data=elldat, aes_string(x="x", y="y", group="stim"), col="gray50", alpha=.25)+ 
        geom_point(aes_string(x="X1", y="X2", colour="group", shape="group"), size=ptSize) 
    }
  }
  g
}
