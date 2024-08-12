#' Multidimensional Scaling with Missing Data (mlsmu6)
#'
#' The `mlsmu6` function performs multidimensional scaling (MDS) on data with missing values using an iterative approach to minimize the error sum of squares. It returns the estimated positions of stimuli and individuals in a lower-dimensional space.
#'
#' @param input A numeric matrix representing the data to be scaled. Rows correspond to individuals, and columns correspond to stimuli. Missing values (NA) are allowed.
#' @param ndim An integer specifying the number of dimensions to estimate. Default is 2.
#' @param cutoff An integer specifying the minimum number of non-missing values required for a row to be included in the analysis. Default is 5.
#' @param tol A numeric value specifying the convergence tolerance for the iterative procedure. Default is 0.0005.
#' @param maxit An integer specifying the maximum number of iterations allowed. Default is 50.
#' @param id An optional vector of identifiers for the rows (individuals) of the input matrix. If provided, it will be added to the output.
#'
#' @return A list of class `mlsmu6` containing:
#' \item{inds}{A data frame with the estimated positions of the individuals in the lower-dimensional space.}
#' \item{stims}{A data frame with the estimated positions of the stimuli in the lower-dimensional space.}
#' \item{iter}{A matrix tracking the error sum of squares over iterations.}
#'
#' @details The function starts by centering the input data, then applies singular value decomposition (SVD) to initialize the positions of the stimuli and individuals. It iteratively adjusts these positions to minimize the error sum of squares. The process stops when the change in the error sum of squares falls below the specified tolerance (`tol`) or the maximum number of iterations (`maxit`) is reached.
#'
#' @examples
#' \dontrun{
#' # Load the interest1981 dataset and prepare the input matrix
#' data(interest1981)
#' input <- as.matrix(interest1981[, 9:38])
#' cutoff <- 5
#' input <- input[rowSums(!is.na(input)) >= cutoff, ]
#' input <- (100 - input) / 50
#' input2 <- input * input
#' input2[is.na(input)] <- (mean(input, na.rm = TRUE))^2
#' inputDC <- doubleCenterRect(input2)
#' 
#' # Perform SVD and initialize positions
#' xsvd <- svd(inputDC)
#' ndim <- 2
#' stims <- xsvd$v[, 1:ndim]
#' inds <- xsvd$u[, 1:ndim]
#' for (i in 1:ndim) {
#'   stims[, i] <- stims[, i] * sqrt(xsvd$d[i])
#'   inds[, i] <- inds[, i] * sqrt(xsvd$d[i])
#' }
#'
#' # Run mlsmu6 to get the MDS solution
#' out <- mlsmu6(input = interest1981[, 9:38], ndim = 2, cutoff = 5,
#'               id = factor(interest1981$party, labels = c("D", "R")))
#' 
#' # Inspect the results
#' print(out$inds)
#' print(out$stims)
#' plot(out$iter, type = "b", main = "Error Sum of Squares Over Iterations")
#' }
#'
#' @export
mlsmu6 <- function(input, ndim = 2, cutoff = 5, tol = 0.0005, maxit = 50, id = NULL){
  rn <- rownames(input)
  cn <- colnames(input)
  input <- as.matrix(input)
  iter <- NULL
  keep <- which(rowSums(!is.na(input)) >= cutoff)
  T <- input[keep, ]
  id <- id[keep]
  #
  np <- nrow(T)
  nq <- ncol(T)
  xrow <- sapply(1:np, function(i) length(rep(1, nq)[!is.na(T[i, ])]))
  xcol <- sapply(1:nq, function(j) length((1:np)[!is.na(T[, j])]))
  #
  T <- (100 - T) / 50
  TT <- T
  TT[is.na(TT)] <- mean(T, na.rm = TRUE)
  TTSQ <- T * T
  TTSQ[is.na(T)] <- (mean(T, na.rm = TRUE))^2
  TEIGHT <- as.numeric(TT)
  dim(TEIGHT) <- c(np, nq)
  TTSQDC <- doubleCenterRect(TTSQ)
  #
  xsvd <- svd(TTSQDC)
  zz <- xsvd$v[, 1:ndim]
  xx <- rep(0, np * ndim)
  dim(xx) <- c(np, ndim)
  for (i in 1:ndim) {
    xx[, i] <- xsvd$u[, i] * sqrt(xsvd$d[i])
  }
  #
  sumaj <- function(i){
    jjj <- 1:nq
    j <- jjj[!is.na(T[i, ])]
    s <- (xx[i, 1] - zz[j, 1])^2 + (xx[i, 2] - zz[j, 2])^2  ### Note this needs to be expanded if estimating more than two dimensions
    s = sqrt(s)
    sx = TEIGHT[i, j]
    sum((s - sx)^2)
  }
  #
  sumvector <- sapply(1:np, sumaj)
  suma <- sum(sumvector)
  xxx <- xx
  zzz <- zz
  kp = 0
  ktp = 0
  sumalast <- 100000.00
  SAVEsumalast <- sumalast
  #
  loop <- done <- 1
  while (done >= tol & loop <= maxit) {
    xxx[, 1:ndim] <- 0
    zzz[, 1:ndim] <- 0
    kp = kp + 1
    ktp = ktp + 1
    #
    for (i in 1:np) {
      for (j in 1:nq) {
        if (!is.na(T[i, j])) {
          s = 0
          for (k in 1:ndim)
            s <- s + (xx[i, k] - zz[j, k])^2
          xc <- ifelse(s == 0, 1.0, TEIGHT[i, j] / sqrt(s))
          for (k in 1:ndim)
            zzz[j, k] = zzz[j, k] + (xx[i, k] - xc * (xx[i, k] - zz[j, k])) / xcol[j]
        }
      }
    }
    for (k in 1:ndim) {
      for (i in 1:np) {
        sw <- 0.0
        for (j in 1:nq) {
          if (!is.na(T[i, j])) {
            s = 0
            for (kk in 1:ndim)
              s <- s + (xx[i, kk] - zzz[j, kk])^2
            xc <- ifelse(s == 0, 1.0, TEIGHT[i, j] / sqrt(s))
            xxx[i, k] = xxx[i, k] + (zzz[j, k] - xc * (zzz[j, k] - xx[i, k]))
          }
        }
        xxx[i, k] = xxx[i, k] / xrow[i]
      }
    }
    xx <- xxx
    zz <- zzz
    sumvector <- sapply(1:np, sumaj)
    suma <- sum(sumvector)
    #
    iter <- c(iter, suma)
    done = ((sumalast - suma) / suma)
    sumalast = suma
    loop <- loop + 1
  }
  #
  # FLIP SPACE IF POLARITY WRONG
  if (xx[1, 1] < 0) {
    xx[, 1] <- -1 * xx[, 1]
    zz[, 1] <- -1 * zz[, 1]
  }
  rownames(xx) <- rn
  rownames(zz) <- cn
  xx <- as.data.frame(xx)
  zz <- as.data.frame(zz)
  names(xx) <- names(zz) <- paste0("Dim", 1:ndim)
  if (!is.null(id)) { xx$id <- id }
  ret <- list(inds = xx, stims = zz, iter = cbind(1:length(iter), iter))
  colnames(ret$iter) <- c("iter", "ErrorSS")
  class(ret) <- "mlsmu6"
  return(ret)
}
