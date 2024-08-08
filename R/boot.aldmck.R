#' Bootstrap Aldmck Analysis
#'
#' This function performs a bootstrap analysis on the Aldmck model, providing confidence intervals for the ideal points of stimuli.
#'
#' @param data A dataframe containing the data to be bootstrapped.
#' @param ... Additional arguments passed to the `aldmck` function.
#' @param boot.args A list of arguments to be passed to the `boot` function.
#' @param plot Logical, whether to plot the results (default is FALSE).
#' @return A list with two elements:
#' \describe{
#'   \item{sumstats}{A dataframe containing the summary statistics: stimulus, ideal point, standard deviation, lower and upper bounds of the confidence interval.}
#'   \item{bootres}{The bootstrap results from the `boot` function.}
#' }
#' @importFrom basicspace aldmck
#' @examples
#' \dontrun{
#' data(franceEES2009)
#' boot.france <- boot.aldmck(franceEES2009, polarity = 2, respondent = 1, missing = c(77, 88, 89), verbose = FALSE, boot.args = list(R = 100))
#' }
#' @export
boot.aldmck <- function(data, ..., boot.args = list(), plot = FALSE) {
  dot.args <- as.list(match.call(expand.dots = FALSE)$`...`)
  boot.fun <- function(data, inds, dot.args, ...) {
    tmp <- data[inds, ]
    dot.args$data <- tmp
    out <- do.call("aldmck", dot.args)
    out$stimuli
  }
  boot.args$data <- data
  boot.args$statistic = boot.fun
  boot.args$dot.args = dot.args
  b <- do.call("boot", boot.args)
  out <- data.frame(
    "stimulus" = factor(names(b$t0), levels = names(b$t0)[order(b$t0)]),
    "idealpt" = b$t0,
    "sd" = apply(b$t, 2, sd)
  )
  rownames(out) <- NULL
  out$lower <- out$idealpt - 1.96 * out$sd
  out$upper <- out$idealpt + 1.96 * out$sd
  out <- out[order(out$idealpt), ]
  class(out) <- c("aldmck_ci", "data.frame")
  return(list(sumstats = out, bootres = b))
}


#' Bootstrap Analysis for Aldmck Model
#'
#' This function performs a bootstrap analysis on the Aldmck model, estimating confidence intervals for the ideal points of stimuli.
#'
#' @param data A dataframe containing the input data for the bootstrap analysis.
#' @param ... Additional arguments to be passed to the `aldmck` function, which is used internally for analysis.
#' @param boot.args A list of arguments to be passed to the `boot` function, such as the number of bootstrap replicates (`R`).
#' @param plot Logical, whether to plot the bootstrap results (default is FALSE).
#' @return A list with two elements:
#' \describe{
#'   \item{sumstats}{A dataframe with summary statistics: stimulus names, estimated ideal points, standard deviations, and confidence intervals (lower and upper bounds).}
#'   \item{bootres}{The full bootstrap results returned by the `boot` function.}
#' }
#' @import basicspace 
#' @examples
#' \dontrun{
#' data(franceEES2009)
#' boot.france <- boot.aldmck(franceEES2009, polarity = 2, respondent = 1, missing = c(77, 88, 89), verbose = FALSE, boot.args = list(R = 100))
#' }
#' @export
boot.aldmck <- function(data, ..., boot.args = list(), plot = FALSE) {
  dot.args <- as.list(match.call(expand.dots = FALSE)$`...`)
  boot.fun <- function(data, inds, dot.args, ...) {
    tmp <- data[inds, ]
    dot.args$data <- tmp
    out <- do.call("aldmck", dot.args)
    out$stimuli
  }
  boot.args$data <- data
  boot.args$statistic = boot.fun
  boot.args$dot.args = dot.args
  b <- do.call("boot", boot.args)
  out <- data.frame(
    "stimulus" = factor(names(b$t0), levels = names(b$t0)[order(b$t0)]),
    "idealpt" = b$t0,
    "sd" = apply(b$t, 2, sd)
  )
  rownames(out) <- NULL
  out$lower <- out$idealpt - 1.96 * out$sd
  out$upper <- out$idealpt + 1.96 * out$sd
  out <- out[order(out$idealpt), ]
  class(out) <- c("aldmck_ci", "data.frame")
  return(list(sumstats = out, bootres = b))
}
