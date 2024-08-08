#' Print Aldmck Confidence Intervals
#'
#' This function formats and prints an `aldmck_ci` object, displaying the ideal points and their associated confidence intervals with specified precision.
#'
#' @param x An object of class `aldmck_ci` containing ideal points and confidence intervals.
#' @param ... Additional arguments passed to other methods (currently not used).
#' @param digits Integer, the number of decimal places to display (default is 3).
#' @return The function prints the formatted `aldmck_ci` object as a data frame.
#' @examples
#' \dontrun{
#' data(result.france)
#' print.aldmck_ci(result.france$sumstats, digits = 2)
#' }
#' @export
print.aldmck_ci <- function(x, ..., digits = 3) {
  x$idealpt <- sprintf(paste0("%.", digits, "f"), x$idealpt)
  x$sd <- sprintf(paste0("%.", digits, "f"), x$sd)
  x$lower <- sprintf(paste0("%.", digits, "f"), x$lower)
  x$upper <- sprintf(paste0("%.", digits, "f"), x$upper)
  print.data.frame(x)
}
