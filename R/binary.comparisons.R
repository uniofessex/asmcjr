#' Binary Comparisons Function
#'
#' This function performs binary comparisons between pairs of columns in a given matrix or data frame.
#' The comparisons are made by subtracting one column from another and determining the sign of the result.
#' The function is specifically designed for use with the `candidatetherms2008` dataset, but can be applied to other datasets with similar structure.
#'
#' @param obj A matrix or data frame where the rows represent observations and the columns represent different candidates or factors to compare. This object is specifically designed to be used with the `candidatetherms2008` dataset.
#' @return A matrix of the same number of rows as `obj`, with columns representing the comparisons between pairs of the original columns. The resulting values are encoded as:
#' The columns of the returned matrix are named according to the pair of columns being compared, joined by an underscore.
#' @examples
#' \dontrun{
#' data(candidatetherms2008)
#' X <- binary.comparisons(candidatetherms2008)
#' print(X[1:5,1:4])
#' }
#' @export
binary.comparisons <- function(obj){
  nr <- nrow(obj)
  nc <- ncol(obj)
  combs <- combn(nc, 2)
  res <- sapply(1:ncol(combs), function(i)
    sign(rowSums(cbind(obj[,combs[2,i]], -obj[,combs[1,i]]))))
  res[which(res == 1)] <- 6
  res[which(res == -1)] <- 1
  res[which(res == 0)] <- 9
  res[which(is.na(res))] <- 9
  colnames(res) <- combn(colnames(obj), 2, paste, collapse = "_")
  return(res)
}
