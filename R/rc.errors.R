#' Calculate Roll Call Errors and Proportional Reduction in Error (PRE)
#'
#' This function calculates the roll call errors and the Proportional Reduction in Error (PRE) statistic for a given roll call vote. The function can handle both "nomObject" and "OCobject" classes, typically representing NOMINATE or Optimal Classification results.
#'
#' @param obj An object of class `nomObject` or `OCobject` containing the results of a NOMINATE or Optimal Classification analysis.
#' @param data A data frame containing the roll call data, including vote outcomes and codes for "yea" and "nay" votes.
#' @param rcnum An integer specifying the roll call number for which errors and PRE should be calculated.
#' @param rotMat A 2x2 rotation matrix used to rotate the coordinates in the case of an `OCobject`. The default is the identity matrix.
#' 
#' @return A list with the following components:
#' \itemize{
#'   \item \code{tot.errors}: A named vector with the total number of "yea" and "nay" errors.
#'   \item \code{PRE}: The Proportional Reduction in Error statistic, calculated as the proportion of errors reduced relative to the worst-case scenario.
#'   \item \code{errors}: A matrix indicating which legislators made errors in their votes, with columns for "yea" and "nay" errors.
#' }
#'
#' @examples
#' \dontrun{
#' data(hr108)
#' library(wnominate)
#' result <- wnominate(hr108, ubeta=15, uweights=0.5, dims=2,
#' minvotes=20, lop=0.025, trials=1, polarity=c(1,5), verbose=FALSE)
#' rc.errors(result, hr108, 528)[c("tot.errors", "PRE")]
#' }
#' @export
rc.errors <- function(obj, data, rcnum, rotMat = diag(2)){
  if("nomObject" %in% class(obj)){
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
    polarity <- X1*N1W + X2*N2W - ws
    errors1 <- (vote %in% data$codes$yea) & polarity >= 0
    errors2 <- (vote %in% data$codes$nay) & polarity <= 0
    errors3 <- (vote %in% data$codes$yea) & polarity <= 0
    errors4 <- (vote %in% data$codes$nay) & polarity >= 0
    kerrors12 <- sum(errors1==1,na.rm=T)+sum(errors2==1,na.rm=T)
    kerrors34 <- sum(errors3==1,na.rm=T)+sum(errors4==1,na.rm=T)
    if (kerrors12 >= kerrors34){
      yeaerror <- errors3
      nayerror <- errors4
    }
    if (kerrors12 < kerrors34){
      yeaerror <- errors1
      nayerror <- errors2
    }
    kerrorsmin <- min(kerrors12,kerrors34)
    errors = cbind(yeaerror, nayerror)
    colnames(errors) <- c("yea_error", "nay_error")
    kpyea <- sum(vote %in% data$codes$yea)
    kpnay <- sum(vote %in% data$codes$nay)
    PRE <- (min(kpyea,kpnay) - kerrorsmin) / min(kpyea,kpnay)
  }
  if("OCobject" %in% class(obj)){
    vote <- as.integer(data$votes[,rcnum])
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
    xws <- ws[rcnum]*N1[rcnum]
    yws <- ws[rcnum]*N2[rcnum]
    N1W <- N1[rcnum]
    N2W <- N2[rcnum]
    ws <- obj$rollcalls[,8]
    xws <- ws[rcnum]*N1[rcnum]
    yws <- ws[rcnum]*N2[rcnum]
    N1W <- N1[rcnum]
    N2W <- N2[rcnum]
    kpyea <- sum(vote==1)
    kpnay <- sum(vote==6)
    polarity <- oc1*N1W + oc2*N2W - ws[rcnum]
    errors1 <- (vote %in% data$codes$yea) & polarity >= 0
    errors2 <- (vote %in% data$codes$nay) & polarity <= 0
    errors3 <- (vote %in% data$codes$yea) & polarity <= 0
    errors4 <- (vote %in% data$codes$nay) & polarity >= 0
    kerrors12 <- sum(errors1==1,na.rm=T) + sum(errors2==1,na.rm=T)
    kerrors34 <- sum(errors3==1,na.rm=T) + sum(errors4==1,na.rm=T)
    if (kerrors12 >= kerrors34){
      yeaerror <- errors3
      nayerror <- errors4
    }
    if (kerrors12 < kerrors34){
      yeaerror <- errors1
      nayerror <- errors2
    }
    kerrorsmin <- min(kerrors12,kerrors34)
    errors = cbind(yeaerror, nayerror)
    PRE <- (min(kpyea,kpnay)-kerrorsmin) / min(kpyea,kpnay)
  }
  ret <- list(tot.errors = colSums(errors), PRE = PRE, errors=errors)
  return(ret)
}

