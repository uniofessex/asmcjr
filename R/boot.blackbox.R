#' Bootstrap Analysis for Blackbox Model
#'
#' This function performs a bootstrap analysis on the Blackbox model.
#'
#' @param data A dataframe containing the data to be bootstrapped.
#' @param missing Vector of missing values.
#' @param dims Integer, the number of dimensions to be considered in the analysis.
#' @param minscale Minimum scale value.
#' @param verbose Logical, whether to print detailed information (default is FALSE).
#' @param posStimulus Integer, the position of the stimulus to check for sign correction (default is 1).
#' @param R Integer, the number of bootstrap replicates (default is 100).
#' @return An array of bootstrapped results with class \code{"bootbb"}.
#' @importFrom basicspace blackbox
#' @examples
#' \dontrun{
#' data(issues.sweden)
#' boot.blackbox(issues.sweden, missing=8, dims=3, minscale=5, verbose=FALSE, posStimulus=13)}
#' @export
boot.blackbox <- function(data, missing, dims=1, minscale, verbose=FALSE, posStimulus = 1, R=100){
  dot.args <- as.list(match.call(expand.dots = FALSE)$`...`)
  orig <- blackbox(data, missing=missing, dims=dims, minscale=minscale, verbose=verbose)
  if(orig$individuals[[dims]][posStimulus, 1] < 0){
    orig$individuals[[dims]][,1] <- -orig$individuals[[dims]][,1]
  }
  sample.dat <- lapply(1:R, function(i)data[,sample(1:ncol(data), ncol(data), replace=TRUE)])
  for(i in 1:length(sample.dat))colnames(sample.dat[[i]]) <- 1:ncol(sample.dat[[i]])
  out <- array(dim=c(nrow(data), dims, R))
  for(i in 1:R){
    tmp <- blackbox(sample.dat[[i]], missing=missing, dims=dims, minscale=minscale, verbose=verbose)
    if(tmp$individuals[[dims]][posStimulus, 1] < 0){
      tmp$individuals[[dims]][,1] <- -tmp$individuals[[dims]][,1]
    }
    out[,,i] <- as.matrix(tmp$individuals[[dims]])
  }
  class(out) <- "bootbb"
  invisible(out)
}