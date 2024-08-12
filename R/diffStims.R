#' Compute the Probabilities of Differences Between Stimuli
#'
#' This function calculates the probabilities that one stimulus is greater than another based on MCMC or resampled values. 
#' It can handle input data in the form of a matrix, data frame, or `mcmc.list`. The function computes pairwise comparisons 
#' between the specified stimuli.
#'
#' @param x A matrix, data frame, or `mcmc.list` containing MCMC or resampled values. The rows represent different samples, 
#' and the columns represent different stimuli or variables.
#' @param stims A vector of stimuli to compare. This can be a numeric vector of column indices or a character vector of column names.
#' @param digits An integer specifying the number of decimal places to display in the output. Default is 3.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with two columns: `Comparison`, which contains the comparison labels (e.g., "Pr(stim2 > stim1)"), 
#' and `Probability`, which contains the probabilities formatted to the specified number of digits.
#' @importFrom utils combn
#' @examples
#' \dontrun{
#' data(bam.france)
#' diffStims(bam.france$zhat, stims=c(3,4))
#' }
#' @export
diffStims <- function(x, stims, digits=3, ...){
  if("mcmc.list" %in% class(x)){
    x <- do.call("rbind", x)
  }
  if(!(is.matrix(x) | is.data.frame(x))) stop("x must be a matrix or data frame of MCMC or resampled values\n")
  x <- as.matrix(x)
  if(!is.numeric(stims)){
    stims <- match(stims, colnames(x))
  }
  combs <- combn(stims, 2)[,,drop=FALSE]
  D <- matrix(0, ncol=ncol(combs), nrow=ncol(x))
  D[cbind(combs[1,], 1:ncol(combs))] <- -1
  D[cbind(combs[2,], 1:ncol(combs))] <- 1
  diffs <- x %*% D
  probs <- colMeans(diffs > 0)
  comps <- paste("Pr(", colnames(x)[combs[2, ]], " > ", colnames(x)[combs[1,]], ")", sep="")
  result <- data.frame('Comparison'=comps, 'Probability'=sprintf(paste0("%.", digits, "f"), probs))
  result
}
