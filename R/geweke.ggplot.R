#' Geweke Diagnostic Plot using ggplot2
#'
#' The `geweke.ggplot` function creates a Geweke diagnostic plot for MCMC chains using `ggplot2`. 
#' It calculates the Geweke diagnostic for segments of the MCMC chain and visualizes the Z-scores over the first iteration in each segment.
#'
#' @param x An object of class `mcmc` or `mcmc.list`. This object represents the MCMC chains for which the Geweke diagnostic will be calculated.
#' @param frac1 A numeric value specifying the fraction of the first segment of the chain to use in the Geweke diagnostic. Defaults to `0.1`.
#' @param frac2 A numeric value specifying the fraction of the last segment of the chain to use in the Geweke diagnostic. Defaults to `0.5`.
#' @param nbins An integer specifying the number of segments (bins) into which the chain should be divided for the diagnostic calculation. Defaults to `20`.
#' @param pvalue A numeric value specifying the p-value for the confidence limits in the Z-score plot. Defaults to `0.05`.
#'
#' @return A `ggplot` object that displays the Geweke diagnostic Z-scores for each segment of the MCMC chains. The plot includes horizontal lines indicating the confidence limits.
#'
#' @details
#' The function splits the MCMC chain into segments, computes the Geweke diagnostic for each segment, and then plots the Z-scores. The Z-scores are compared against the standard normal quantiles, and segments outside the confidence limits indicate potential non-convergence.
#'
#' @importFrom coda as.mcmc.list geweke.diag nvar nchain varnames chanames
#' @importFrom stats start end qnorm window
#' @importFrom graphics symbols
#' @examples
#' \dontrun{
#' library(coda)
#' data(line)
#' geweke_plot <- geweke.ggplot(line)
#' print(geweke_plot)
#' }
#' @export
geweke.ggplot <- function (x, frac1 = 0.1, frac2 = 0.5, nbins = 20, pvalue = 0.05)
{
  x <- as.mcmc.list(x)
  ystart <- seq(from = start(x), to = (start(x) + end(x))/2,
                length = nbins)
  gcd <- array(dim = c(length(ystart), nvar(x), nchain(x)),
               dimnames = list(ystart, varnames(x), chanames(x)))
  for (n in 1:length(ystart)) {
    geweke.out <- geweke.diag(window(x, start = ystart[n]),
                              frac1 = frac1, frac2 = frac2)
    for (k in 1:nchain(x)) gcd[n, , k] <- geweke.out[[k]]$z
  }
  climit <- qnorm(1 - pvalue/2)
  tmp.df <- data.frame(
    gcd = c(gcd), 
    chain = as.factor(rep(1:nchain(x), each=nrow(gcd)*ncol(gcd))), 
    var = factor(rep(rep(1:nvar(x), each=nrow(gcd)), dim(gcd)[3])), 
    x=rep(ystart, dim(gcd)[2]*dim(gcd)[3]))
  
  faclabs <- colnames(x[[1]])
  if(!is.null(faclabs)){
    levels(tmp.df$var) <- faclabs
  }
  
  g <- ggplot(tmp.df) + 
    geom_point(aes_string(x="x", y="gcd", colour="chain")) + 
    geom_hline(yintercept=c(-climit, climit), lty=2, size=.5) + 
    labs(x="First Iteration in Segment", y="Z-score")
  if(length(unique(tmp.df$var)) > 1){
    g <- g + facet_wrap(~var)
  }
  g
}
