#' Bootstrap Analysis for the Blackbox Transpose Model
#'
#' This function performs a bootstrap analysis on the Blackbox Transpose model, generating bootstrap replicates of the ideal points for stimuli.
#'
#' @param data A matrix containing the input data for the analysis.
#' @param missing Numeric vector indicating which values are treated as missing.
#' @param dims Integer, the number of dimensions to analyze (default is 1).
#' @param minscale Numeric, the minimum scale value used in the analysis.
#' @param verbose Logical, whether to print detailed information during execution (default is FALSE).
#' @param posStimulus Integer, the index of the stimulus used to check and correct the sign (default is 1).
#' @param R Integer, the number of bootstrap replicates (default is 100).
#' @return An array containing the bootstrapped results.
#' @importFrom basicspace blackbox_transpose
#' @examples
#' \dontrun{
#' data(issues.sweden)
#' result <- boot.blackbox_transpose(issues.sweden, missing = 8, 
#' dims = 2, minscale = 5, verbose = FALSE, posStimulus = 13, R = 100)
#' }
#' @export

boot.blackbox_transpose <- function(data, missing, dims = 1, minscale, verbose = FALSE, posStimulus = 1, R = 100) {
  out <- array(dim = c(ncol(data), dims, R))
  for (i in 1:R) {
    tmp <- data[sample(1:nrow(data), nrow(data), replace = TRUE), ]
    result <- blackbox_transpose(tmp, missing = missing, dims = dims, minscale = minscale, verbose = verbose)
    if (result$stimuli[[dims]][posStimulus, 2] > 0)
      result$stimuli[[dims]][, 2] <- -1 * result$stimuli[[dims]][, 2]
    out[, , i] <- as.matrix(result$stimuli[[dims]][, 2:((2 + dims) - 1)])
  }
  return(out)
}