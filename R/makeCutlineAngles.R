#' Calculate Cutline Angles
#'
#' The `makeCutlineAngles` function calculates the cutline angles for a given object of class `wnominate`, which can either be a `nomObject` or `OCobject`. The function processes the roll call data and computes the angles associated with the cutlines.
#'
#' @param obj An object of class `wnominate`, which can either be a `nomObject` or `OCobject`. This object contains roll call data and other relevant information needed to calculate the cutline angles.
#'
#' @return A data frame with three columns:
#' \describe{
#'   \item{angle}{The calculated cutline angles in degrees.}
#'   \item{N1W}{The adjusted first component of the normal vector to the cutline.}
#'   \item{N2W}{The adjusted second component of the normal vector to the cutline.}
#' }
#'
#' @details
#' The function differentiates between objects of class `nomObject` and `OCobject` to perform appropriate calculations. It uses the roll call data within the object to compute the normal vectors and their corresponding angles.
#'
#' @examples
#' \dontrun{
#' # Assuming `wnominate_obj` is a valid wnominate object of class `nomObject` or `OCobject`
#' result <- makeCutlineAngles(wnominate_obj)
#' print(result)
#' }
#'
#' @keywords datasets
#' @export
makeCutlineAngles <- function(obj){
  if("nomObject" %in% class(obj)){
    WEIGHT <- (obj$weights[2])/(obj$weights[1])
    DL1 <- obj$rollcalls[,7]
    DL2 <- obj$rollcalls[,8]
    ZM1 <- obj$rollcalls[,9]
    ZM2 <- obj$rollcalls[,10]
    YEA1 <- ZM1 - DL1
    YEA2W <- (ZM2 - DL2) * WEIGHT
    NAY1 <- ZM1 + DL1
    NAY2W <- (ZM2 + DL2) * WEIGHT
    A1 <- NAY1 - YEA1
    A2 <- NAY2W - YEA2W
    ALENGTH <- sqrt(A1*A1 + A2*A2)
    N1W <- A1 / ALENGTH
    N2W <- A2 / ALENGTH
    for (i in 1:nrow(obj$rollcalls)){
      if (N1W[i] < 0 & !is.na(N2W[i])) N2W[i] <- -N2W[i]
      if (N1W[i] < 0 & !is.na(N1W[i])) N1W[i] <- -N1W[i]
    }
    C1 <- N2W
    C2 <- -N1W
    for (i in 1:nrow(obj$rollcalls)){
      if (C1[i] < 0 & !is.na(C2[i])) C2[i] <- -C2[i]
      if (C1[i] < 0 & !is.na(C1[i])) C1[i] <- -C1[i]
    }
    theta <- atan2(C2,C1)
    theta4 <- theta * (180/pi)
  }
  if("OCobject" %in% class(obj)){
    oc1 <- obj$legislators[,7]
    oc2 <- obj$legislators[,8]
    PRE <- obj$rollcalls[,5]
    N1 <- obj$rollcalls[,6]
    N2 <- obj$rollcalls[,7]
    ws <- obj$rollcalls[,8]
    xws <- ws * N1
    yws <- ws * N2
    C1 <- N2
    C2 <- -N1
    for (i in 1:nrow(obj$rollcalls)){
      if (C1[i] < 0 & !is.na(C2[i])) C2[i] <- -C2[i]
      if (C1[i] < 0 & !is.na(C1[i])) C1[i] <- -C1[i]
    }
    theta <- atan2(C2,C1)
    theta4 <- theta * (180/pi)
    N1W <- N1
    N2W <- N2
  }
  res <- data.frame(angle = theta4, N1W = N1W, N2W = N2W)
  return(res)
}
