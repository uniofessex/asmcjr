#' Calculate Standard Errors for Aldrich-McKelvey Scaling
#'
#' This function calculates the standard errors for the Aldrich-McKelvey scaling results based on the provided object and data.
#'
#' @param obj An object containing the Aldrich-McKelvey scaling results. This object should contain respondent parameters (e.g., alpha and beta) and stimuli positions.
#' @param data A matrix or dataframe containing the original data used in the Aldrich-McKelvey scaling.
#' @param ... Additional arguments (currently not used).
#' @return A numeric vector of standard errors for each stimulus.
##' @importFrom stats na.omit
#' @examples
#' \dontrun{
#' data(result.france)
#' data(franceEES2009)
#' se_values <- aldmckSE(result.france, franceEES2009)
#' print(se_values)
#' }
#' @export

# aldmckSE <- function(obj, data, ...) {
#   tmp <- na.omit(cbind(obj$respondents[, 1:2], data))
#   alpha <- tmp[, 1]
#   beta <- tmp[, 2]
#   z <- tmp[, 3:ncol(tmp)]
#   zhat <- obj$stimuli
#   sigmaj <- rep(0, length(zhat))
#   
#   # Generate sigma_j
#   for (j in 1:length(zhat)) {
#     for (i in 1:length(alpha)) {
#       sigmaj[j] <- sigmaj[j] + ((alpha[i] + beta[i] * zhat[j]) - z[i, j])^2
#     }
#   }
#   
#   # Calculate standard errors
#   for (i in 1:length(zhat)) {
#     sigmaj[i] <- sqrt(sigmaj[i] / length(alpha))
#   }
#   
#   sigmaj
# }

aldmckSE <- function(obj, data, ...) {
  # 检查 obj 是否包含必要的字段
  if (!is.list(obj) || !all(c("respondents", "stimuli") %in% names(obj))) {
    stop("obj must be a list containing 'respondents' and 'stimuli'.")
  }
  
  # 检查 obj$respondents 是否为矩阵，且包含至少两列
  if (!is.matrix(obj$respondents) || ncol(obj$respondents) < 2) {
    stop("obj$respondents must be a matrix with at least two columns.")
  }
  
  # 检查 data 是否为数值型数据框或矩阵
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("data must be a numeric data frame or matrix.")
  }
  
  # 确保 data 是数值型
  if (!all(sapply(data, is.numeric))) {
    stop("All columns in data must be numeric.")
  }
  
  tmp <- na.omit(cbind(obj$respondents[, 1:2], data))
  alpha <- tmp[, 1]
  beta <- tmp[, 2]
  z <- tmp[, 3:ncol(tmp)]
  zhat <- obj$stimuli
  sigmaj <- rep(0, length(zhat))
  
  # Generate sigma_j
  for (j in 1:length(zhat)) {
    for (i in 1:length(alpha)) {
      sigmaj[j] <- sigmaj[j] + ((alpha[i] + beta[i] * zhat[j]) - z[i, j])^2
    }
  }
  
  # Calculate standard errors
  for (i in 1:length(zhat)) {
    sigmaj[i] <- sqrt(sigmaj[i] / length(alpha))
  }
  
  sigmaj
}
