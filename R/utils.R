#' Generate a Gray Color Palette
#'
#' This function generates a palette of gray colors.
#'
#' @param n Integer, the number of gray shades to generate.
#' @param lower Numeric, the lower bound for the gray scale (default is 0.3).
#' @param upper Numeric, the upper bound for the gray scale (default is 0.7).
#' @return A vector of gray colors.
#' @importFrom grDevices rgb
#' @examples
#' gray.palette(5)
#' @export

gray.palette <- function(n, lower = .3, upper = .7) {
  s <- seq(lower, upper, length = n)
  rgb(matrix(rep(s, each = 3), ncol = 3, byrow = TRUE))
}


#' Get the End Time of an MCMC List
#'
#' This function returns the end time or the last element of the first chain in an MCMC list.
#'
#' @param x An object of class `mcmc.list`, which is a list of MCMC chains.
#' @param ... Additional arguments (currently not used).
#'
#' @return The result of calling the `end` function on the first element of the MCMC list.
#' 
#' @examples
#' # Assuming `mcmc_list` is an object of class `mcmc.list`
#' # end_time <- end.mcmc.list(mcmc_list)
#' 
#' @export
end.mcmc.list <- function (x, ...){
  end(x[[1]])
}


#' Get the Start Time of an MCMC List
#'
#' This function returns the start time or the first element of the first chain in an MCMC list.
#'
#' @param x An object of class `mcmc.list`, which is a list of MCMC chains.
#' @param ... Additional arguments (currently not used).
#'
#' @return The result of calling the `start` function on the first element of the MCMC list.
#' 
#' @examples
#' # Assuming `mcmc_list` is an object of class `mcmc.list`
#' # start_time <- start.mcmc.list(mcmc_list)
#' 
#' @export
start.mcmc.list <- function (x, ...){
  start(x[[1]])
}


#' Apply Windowing to an MCMC List
#'
#' This function applies the `window.mcmc` function to each chain in an MCMC list, allowing you to subset the chains based on the desired window.
#'
#' @param x An object of class `mcmc.list`, which is a list of MCMC chains.
#' @param ... Additional arguments passed to the `window.mcmc` function, such as `start`, `end`, and `thin`.
#'
#' @return An object of class `mcmc.list` with the same structure as the input but with each chain windowed according to the specified parameters.
#' 
#' @examples
#' # Assuming `mcmc_list` is an object of class `mcmc.list`
#' # Subset the MCMC list to include only iterations between 500 and 1000
#' windowed_mcmc_list <- window.mcmc.list(mcmc_list, start = 500, end = 1000)
#' 
#' @export
window.mcmc.list <- function (x, ...) {
  structure(lapply(x, window.mcmc, ...), class = "mcmc.list")
}

#' Subset an MCMC Object by Specifying a Time Window
#'
#' The `window.mcmc` function subsets a Markov Chain Monte Carlo (MCMC) object by selecting iterations within a specified time window.
#'
#' @param x An object of class `mcmc`, representing a single MCMC chain.
#' @param start The starting iteration (time) for the subset. If not provided, the starting time of the original chain is used.
#' @param end The ending iteration (time) for the subset. If not provided, the ending time of the original chain is used.
#' @param thin The thinning interval to be applied. If not provided, the original thinning interval is used. If `thin` is not a multiple of the original thinning interval, a warning is issued and the original thinning interval is retained.
#' @param ... Additional arguments (currently not used).
#'
#' @return A new `mcmc` object representing the subsetted chain with the specified start, end, and thinning parameters.
#'
#' @details
#' This function allows you to focus on a specific window of iterations within an MCMC chain. It adjusts the start and end points according to the provided values, ensuring they align with the chain's iterations. The thinning interval can also be adjusted, but it must be a multiple of the original thinning interval; otherwise, a warning is issued and the original interval is kept.
#'
#' @examples
#' # Assume `mcmc_chain` is an object of class `mcmc`
#' # Subset the MCMC chain to include only iterations between 500 and 1000
#' windowed_chain <- window.mcmc(mcmc_chain, start = 500, end = 1000, thin = 2)
#'
#' @export
window.mcmc <- function (x, start, end, thin, ...) {
  ts.eps <- getOption("ts.eps")
  xmcpar <- mcpar(x)
  xstart <- xmcpar[1]
  xend <- xmcpar[2]
  xthin <- xmcpar[3]
  if (missing(thin)) 
    thin <- xthin
  else if (thin%%xthin != 0) {
    thin <- xthin
    warning("Thin value not changed")
  }
  xtime <- as.vector(time(x))
  if (missing(start)) 
    start <- xstart
  else if (length(start) != 1) 
    stop("bad value for start")
  else if (start < xstart) {
    start <- xstart
    warning("start value not changed")
  }
  if (missing(end)) 
    end <- xend
  else if (length(end) != 1) 
    stop("bad value for end")
  else if (end > xend) {
    end <- xend
    warning("end value not changed")
  }
  if (start > end) 
    stop("start cannot be after end")
  if (all(abs(xtime - start) > abs(start) * ts.eps)) {
    start <- xtime[(xtime > start) & ((start + xthin) > xtime)]
  }
  if (all(abs(end - xtime) > abs(end) * ts.eps)) {
    end <- xtime[(xtime < end) & ((end - xthin) < xtime)]
  }
  use <- 1:niter(x)
  use <- use[use >= trunc((start - xstart)/xthin + 1.5) & use <= 
               trunc((end - xstart)/xthin + 1.5) & (use - trunc((start - 
                                                                   xstart)/xthin + 1.5))%%(thin%/%xthin) == 0]
  y <- if (is.matrix(x)) 
    x[use, , drop = FALSE]
  else x[use]
  return(mcmc(y, start = start, end = end, thin = thin))
}