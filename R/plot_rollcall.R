#' Plot Roll Call Results
#'
#' This function creates a scatter plot of legislators' ideal points for a given roll call vote. It allows for the visualization of voting patterns, highlighting errors, and optionally adding shape variation for different groups of legislators.
#'
#' @param obj An object of class `nomObject` containing the results of a NOMINATE analysis, including legislators' ideal points and roll call data.
#' @param data A data frame containing the roll call data, including vote outcomes and codes for "yea" and "nay" votes.
#' @param gdat A data frame containing additional data related to the legislators. (This parameter is currently not used in the function but can be included for future extensions.)
#' @param rcnum An integer specifying the roll call number to be plotted.
#' @param shapeVar An optional variable used to differentiate groups of legislators by shape in the plot. If provided, the points will be shaped according to this variable.
#' @param dropNV A logical value indicating whether to drop non-voters (NA values) from the plot. Default is `FALSE`.
#' @param onlyErrors A logical value indicating whether to plot only the legislators who made voting errors. Default is `FALSE`.
#' @param ptSize A numeric value specifying the size of the points in the plot. Default is `4`.
#' 
#' @return A `ggplot` object showing the roll call results, including the voting patterns of legislators.
#' 
#' @examples
#' \dontrun{
#' g <- plot_rollcall(result, hr108, gdat = NULL, rcnum = 528)
#' print(g)
#' }
#' @import ggplot2
#' @export
plot_rollcall <- function(obj, data, gdat, rcnum, 
                          shapeVar = NULL, dropNV=FALSE, onlyErrors=FALSE, ptSize=4){
  WEIGHT <- (obj$weights[2])/(obj$weights[1])
  X1 <- obj$legislators$coord1D
  X2 <- (obj$legislators$coord2D)*WEIGHT
  vote <- as.integer(data$votes[,rcnum])
  DL1 <- obj$rollcalls[rcnum,7]
  DL2 <- obj$rollcalls[rcnum,8]
  ZM1 <- obj$rollcalls[rcnum,9]
  ZM2 <- obj$rollcalls[rcnum,10]
  YEA1 <- ZM1-DL1
  YEA2W <- (ZM2-DL2)*WEIGHT
  NAY1 <- ZM1+DL1
  NAY2W <- (ZM2+DL2)*WEIGHT
  A1 <- NAY1 - YEA1
  A2 <- NAY2W - YEA2W
  ALENGTH <- sqrt(A1*A1 + A2*A2)
  N1W <- A1/ALENGTH
  N2W <- A2/ALENGTH
  if (N1W < 0){
    N1W <- -N1W
    N2W <- -N2W
  }
  ws <- N1W*ZM1 + N2W*ZM2*WEIGHT
  xws <- ws*N1W
  yws <- ws*N2W
  errors <- rc.errors(obj, data, rcnum)$errors
  df <- data.frame(X1=X1, X2=X2)
  df$vote <- NA
  df$vote[which(vote %in% data$codes$yea)] <- 2
  df$vote[which(vote %in% data$codes$nay)] <- 1
  df$vote <- factor(df$vote, levels=1:2, labels=c("Nay", "Yea"))
  if(!is.null(shapeVar)){
    df$group <- shapeVar
  }
  if(onlyErrors){
    df <- df[which(errors[,1] | errors[,2]), ]
  }
  if(dropNV){
    df <- na.omit(df)
  }
  if(is.null(shapeVar)){
    g <- ggplot(df, aes_string(x="X1", y="X2", colour="vote")) + geom_point(size=ptSize) + 
      scale_colour_manual(values=c("gray67", "gray33"), name="Vote")
  }
  if(!is.null(shapeVar)){
    g <- ggplot(df, aes_string(x="X1", y="X2", colour="vote", shape="group")) + geom_point(size=ptSize) + 
      scale_colour_manual(values=c("gray67", "gray33"), name="Vote")
  }
  g + geom_segment(aes(x=xws+N2W, y=yws-N1W, xend=xws-N2W, yend=yws+N1W))
}