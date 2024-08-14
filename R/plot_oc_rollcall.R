#' Plot Optimal Classification (OC) Coordinates for a Specific Roll Call
#'
#' This function creates a scatter plot of legislators' ideal points based on the results of an Optimal Classification (OC) analysis for a specific roll call. The coordinates can be optionally rotated using a rotation matrix and plotted with different shapes and colors for different groups.
#'
#' @param obj An object of class `OCobject` containing the results of an Optimal Classification analysis, including legislators' ideal points.
#' @param data A `rollcall` object containing the voting data corresponding to the OC analysis.
#' @param shapeVar An optional variable used to differentiate groups of legislators by shape and color in the plot. If provided, the points will be shaped and colored according to this variable.
#' @param rcnum An integer specifying the roll call number to be visualized.
#' @param rotMat A 2x2 rotation matrix used to rotate the coordinates. The default is the identity matrix (no rotation).
#' @param dropNV A logical value indicating whether to drop non-voters (NA values) from the plot. Default is `FALSE`.
#' @param onlyErrors A logical value indicating whether to plot only those observations where classification errors occurred. Default is `FALSE`.
#' @param ptSize A numeric value specifying the size of the points in the plot. Default is `4`.
#' 
#' @return A `ggplot` object displaying the legislators' coordinates for the specified roll call, optionally with different shapes and colors for different groups.
#' 
#' @examples
#' \dontrun{
#' # Assuming `oc_result` is an OCobject, `rc_data` is a rollcall object,
#' # and `legislator_data` contains a grouping variable
#' g <- plot_oc_rollcall(oc_result, data = rc_data, 
#' shapeVar = legislator_data$group, rcnum = 5, dropNV = TRUE)
#' print(g)
#' }
#' @import ggplot2
#' @export
plot_oc_rollcall <- function(obj, data, shapeVar = NULL, rcnum, rotMat=diag(2),
                             dropNV=FALSE, onlyErrors=FALSE, ptSize=4){
  nrollcall <- rcnum
  vote <- as.integer(data$votes[,nrollcall])
  A <- rotMat
  rot.oc <- cbind(obj$legislators$coord1D, obj$legislators$coord2D)
  for (i in 1:nrow(rot.oc)){
    rot.oc[i,] <- rot.oc[i,] %*% A
  }
  oc1 <- rot.oc[,1]
  oc2 <- rot.oc[,2]
  rot.NV <- cbind(obj$rollcalls[,6], obj$rollcalls[,7])
  for (i in 1:nrow(rot.NV)){
    rot.NV[i,] <- rot.NV[i,] %*% A
  }
  N1 <- rot.NV[,1]
  N2 <- rot.NV[,2]
  ws <- obj$rollcalls[,8]
  xws <- ws[nrollcall]*N1[nrollcall]
  yws <- ws[nrollcall]*N2[nrollcall]
  N1W <- N1[nrollcall]
  N2W <- N2[nrollcall]
  errors <- rc.errors(obj, data, rcnum)$errors
  df <- data.frame(X1=oc1, X2=oc2)
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
