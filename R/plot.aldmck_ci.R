#' Plot Confidence Intervals for Aldmck CI
#'
#' This function creates a plot of the ideal points and their confidence intervals for an object of class `aldmck_ci`.
#'
#' @param x An object of class `aldmck_ci` containing ideal points and confidence intervals.
#' @param ... Additional arguments passed to `ggplot2` functions.
#' @return A ggplot2 object representing the confidence intervals for the ideal points of stimuli.
#' @import ggplot2
#' @examples
#' \dontrun{
#' data(result.france)
#' plot.aldmck_ci(result.france$sumstats)
#' }
#' @export
plot.aldmck_ci <- function(x, ...) {
  g <- ggplot(x, aes(x = idealpt, y = stimulus)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
    theme_bw()
  return(g)
}
