#' #' Prepare Data for BAM (Bayesian Aldrich-McKelvey) Scaling
#' #'
#' #' This function preprocesses data for use in Bayesian Aldrich-McKelvey (BAM) scaling. It handles missing values, centers the data, and separates the self-placement from the stimuli placements.
#' #'
#' #' @param x A numeric matrix or dataframe containing the data to be preprocessed.
#' #' @param nmin Integer, the minimum number of non-missing values required per row (default is 1).
#' #' @param missing A vector of values that should be treated as missing (e.g., `c(77, 88, 99)`).
#' #' @param self Integer, the column number representing self-placement (default is 1).
#' #' @param midpt Numeric, the midpoint used for centering the data. If NULL, the midpoint is calculated as the mean of the minimum and maximum values.
#' #' @return A list with two components: `stims` (a matrix of stimuli placements) and `self` (a vector of self-placements).
#' #' @examples
#' #' \dontrun{
#' #' data(franceEES2009)
#' #' bamdata <- bamPrep(franceEES2009, missing=c(77,88,89), self=1)
#' #' }
#' #' @export
#' bamPrep <- function(x, nmin = 1, missing = NULL, self = 1, midpt = NULL) {
#'   x <- as.matrix(x)
#'   if (!is.numeric(x[, 1])) { stop("x must be a numeric data frame or matrix") }
#'   x[which(x %in% missing, arr.ind = TRUE)] <- NA
#'   if (is.null(midpt)) {
#'     x <- apply(x, 2, function(z) z - (min(z, na.rm = TRUE) + diff(range(z, na.rm = TRUE)) / 2))
#'   } else {
#'     x <- apply(x, 2, function(z) z - midpt)
#'   }
#'   nonmiss <- apply(x, 1, function(z) sum(!is.na(z)))
#'   x <- x[which(nonmiss >= nmin), ]
#'   out <- list(stims = x[, -self], self = x[, self])
#'   class(out) <- c("bamPrep", "list")
#'   out
#' }
