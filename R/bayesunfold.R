#' Bayesian Unfolding with Slice Sampling and L-BFGS Optimization
#'
#' The `bayesunfold` function performs Bayesian unfolding using slice sampling and L-BFGS optimization. 
#' It estimates the positions of stimuli and individuals in a lower-dimensional space, 
#' allowing for the analysis of preference data with missing values.
#'
#' @param input A numeric matrix representing the input data, where rows correspond to respondents and 
#' columns correspond to stimuli. Missing values (NA) are allowed.
#' @param dims An integer specifying the number of dimensions to estimate. Default is 2.
#' @param nsamp An integer specifying the number of samples to draw from the posterior distribution. Default is 2000.
#' @param burnin An integer specifying the number of burn-in iterations for the sampler. Default is 1000.
#' @param cred.level A numeric value specifying the credible interval level (between 0 and 1). Default is 0.9.
#' @param slice.starts A character string specifying the starting points for the slice sampler. 
#' Options are `"lbfgs"` (default) for using the L-BFGS result as the starting point, or `"random"` for random starts.
#' @param print.lbfgs A character string specifying the output method for L-BFGS results. 
#' Options are `"console"` (default) or a file path to save the results.
#' @param print.slice A character string specifying the output method for slice sampling results. 
#' Options are `"console"` (default) or a file path to save the results.
#' @param ... Additional arguments passed to the `smacofRect` function for SMACOF metric unfolding.
#'
#' @return A list of class `bayesunfold` containing:
#' \describe{
#'   \item{retained.obs}{A logical vector indicating which rows (respondents) were retained based on the cutoff criteria.}
#'   \item{smacof.result}{The result of SMACOF metric unfolding for comparison purposes.}
#'   \item{lbfgs.result}{A list with the L-BFGS estimated positions for stimuli and individuals.}
#'   \item{samples}{The raw samples from the slice sampler.}
#'   \item{result4}{The result object from the slice sampler.}
#'   \item{sigma_squared_hat}{The estimated variance from the slice sampler.}
#'   \item{sigma_squared_hat_sd}{The standard deviation of the estimated variance from the slice sampler.}
#'   \item{unrotated}{A list containing unrotated results for stimuli and individuals, including means and credible intervals.}
#'   \item{rotated}{A list containing Procrustes-rotated results for stimuli and individuals, including means and credible intervals.}
#' }
#'
#' @details
#' The `bayesunfold` function implements a Bayesian approach to unfolding analysis, 
#' where respondents' preferences are modeled in a latent space. The method combines L-BFGS optimization for initial estimates and slice sampling for posterior inference. 
#' The function handles missing data by excluding respondents with insufficient responses and imputing missing values for the remaining data.
#'
#' The L-BFGS and slice sampling steps can be configured to either print results to the console or save them to specified files.
#' @importFrom smacof smacofRect
#' @importFrom plyr aaply
#' @importFrom coda mcmc
#' @importFrom vegan procrustes
#'
#' @examples
#' \dontrun{
#' # Load the ANES1968 dataset
#' data(anes.input)
#' # Perform Bayesian unfolding
#' result <- bayesunfold(anes.input, dims = 2, 
#' nsamp = 20, burnin = 5, cred.level = 0.9, slice.starts = "lbfgs")
#' }
#' @name bayesunfold
#' @export
# Load the shared object containing the C functions
# dyn.load("src/lbfgs_bu3.so")
# bayesunfold <- function(input, dims = 2, nsamp = 2000, burnin = 1000, cred.level = 0.9,
#                         slice.starts = c("lbfgs", "random"), print.lbfgs="console",
#                         print.slice="console", ...) {
#   # Input validation
#   if (cred.level > 1 | cred.level < 0) {
#     stop("cred.level must be between 0 and 1\n")
#   }
#   tailprob <- (1 - cred.level) / 2
#   ll <- tailprob
#   ul <- 1 - tailprob
#   slice.starts <- match.arg(slice.starts)
#   
#   # Helper functions
#   set_globals <- function(nslice, nburn, nrowX, ncolX, NS, N, NDIM, UNFOLD, NMISSING, X, CONSTRAINTS) {
#     .C("copyFromR",
#        as.integer(nslice),
#        as.integer(nburn),
#        as.integer(nrowX),
#        as.integer(ncolX),
#        as.integer(NS),
#        as.integer(N),
#        as.integer(NDIM),
#        as.integer(UNFOLD),
#        as.integer(NMISSING),
#        as.double(X),
#        as.double(CONSTRAINTS))
#   }
# 
#   do_lbfgs <- function(kpnp, kpnq, yrotate, rmatrix) {
#     .C("mainlbfgs",
#        as.integer(kpnp),
#        as.integer(kpnq),
#        as.double(yrotate),
#        as.double(rmatrix))
#   }
# 
#   do_logposterior <- function(theta, XCOORDS, sumsquared, SIGMAPRIOR) {
#     .C("keithrules",
#        as.double(theta),
#        as.double(XCOORDS),
#        as.double(sumsquared),
#        as.double(SIGMAPRIOR))
#   }
# 
#   do_sliceu <- function(theta, thetanow2, theta1000, ssenow, XTRUE, thetaLeft, thetaRight, WW, PP, XCOORDS, SIGMAPRIOR) {
#     .C("sliceunfolding",
#        as.double(theta),
#        as.double(thetanow2),
#        as.double(theta1000),
#        as.double(ssenow),
#        as.double(XTRUE),
#        as.double(thetaLeft),
#        as.double(thetaRight),
#        as.double(WW),
#        as.integer(PP),
#        as.double(XCOORDS),
#        as.double(SIGMAPRIOR))
#   }
# 
#   # Data preprocessing
#   T <- as.matrix(input)
#   T[T < 0 | T > 100] <- NA
#   nstimuli <- ncol(T)
#   cutoff <- 7
#   keep <- rowSums(!is.na(T)) >= cutoff
#   T <- T[keep,]
#   T <- (100-T)/50
# 
#   # Setup parameters
#   nrowX <- nrow(T)
#   ncolX <- ncol(T)
#   nburn <- burnin
#   nslice <- nsamp
#   NS <- dims
#   N <- NS*(nrowX+ncolX) - ((NS*(NS+1))/2)
#   NDIM <- NS*(nrowX+ncolX) - (NS-1)
#   UNFOLD <- 1
#   NMISSING <- 7
# 
#   TT <- T
#   TT[which(is.na(TT), arr.ind=TRUE)] <- -999.0
#   X <- as.vector(t(TT))
# 
#   # Set up constraints
#   CONSTRAINTS <- rep(1, NS*(nrowX+ncolX))
#   if (NS==1) {
#     CONSTRAINTS[NS*(nrowX+ncolX)] <- 0
#   } else if (NS==2) {
#     CONSTRAINTS[(NS*(nrowX+ncolX)-NS):(NS*(nrowX+ncolX))] <- 0
#   } else if (NS==3) {
#     CONSTRAINTS[(NS*(nrowX+ncolX)-4):(NS*(nrowX+ncolX))] <- 0
#     CONSTRAINTS[(NS*(nrowX+ncolX)-6)] <- 0
#   }
# 
#   # SMACOF (metric unfolding)
#   weights <- matrix(1, nrow=nrowX, ncol=ncolX)
#   weights[is.na(T)] <- 0
#   SMACOF.result <- smacofRect(T, ndim=NS, circle = "none", weightmat=weights, itmax=10000, ...)
# 
#   # Prepare for L-BFGS
#   zmetric <- as.numeric(t(SMACOF.result$conf.col))
#   xmetric <- as.numeric(t(SMACOF.result$conf.row))
#   rmatrix <- c(zmetric, xmetric)
#   rmatrix[(NS*(nrowX+ncolX)-NS):(NS*(nrowX+ncolX))] <- 0
#   yrotate <- rep(0, (NS*(nrowX+ncolX)))
# 
#   # L-BFGS optimization
#   set_globals(nslice, nburn, nrowX, ncolX, NS, N, NDIM, UNFOLD, NMISSING, X, CONSTRAINTS)
#   if(print.lbfgs != "console") {
#     sink(print.lbfgs)
#     lbfgs.result <- do_lbfgs(nrowX, ncolX, yrotate, rmatrix)
#     sink()
#   } else {
#     lbfgs.result <- do_lbfgs(nrowX, ncolX, yrotate, rmatrix)
#   }
# 
#   lbfgs.coords <- lbfgs.result[[3]]
#   dim(lbfgs.coords) <- c(NS, (nrowX+ncolX))
#   X3 <- t(lbfgs.coords)
#   lbfgs.stimuli <- X3[1:ncolX,]
#   lbfgs.individuals <- X3[(ncolX+1):(nrowX+ncolX),]
# 
#   # Bayesian Unfolding (Slice Sampling)
#   if(slice.starts == "random") {
#     theta <- runif(NDIM, min=-.5, max=.5) # Random starts for Slice Sampler
#   } else if(slice.starts == "lbfgs") {
#     theta <- lbfgs.result[[3]] # Non-random starts for Slice Sampler
#   }
# 
#   theta2 <- theta
#   theta1000 <- rep(0, nslice*NDIM)
#   dim(theta1000) <- c(nslice*NDIM, 1)
#   ssenow <- rep(0, (2*(nslice+nburn)))
#   dim(ssenow) <- c((2*(nslice+nburn)), 1)
#   XTRUE <- lbfgs.result[[3]]
#   thetaL <- rep(-99.0, NDIM)
#   thetaR <- rep(99.0, NDIM)
#   dim(thetaL) <- dim(thetaR) <- c(NDIM, 1)
#   thetaL[NDIM] <- 0.10
#   thetaR[NDIM] <- 0.50
#   WW <- 1.0
#   PP <- 3.0
#   XCOORDS <- rep(0, (nrowX+ncolX)*NS)
#   SIGMAPRIOR <- 100.0
# 
#   if(print.slice != "console") {
#     sink(print.slice)
#     result4 <- do_sliceu(theta, theta2, theta1000, ssenow, XTRUE, thetaL, thetaR, WW, PP, XCOORDS, SIGMAPRIOR)
#     sink()
#   } else {
#     result4 <- do_sliceu(theta, theta2, theta1000, ssenow, XTRUE, thetaL, thetaR, WW, PP, XCOORDS, SIGMAPRIOR)
#   }
# 
#   # Post-processing
#   sigma_squared_hat <- mean(result4[[4]][(burnin+1):(burnin+nsamp)])
#   sigma_squared_hat_sd <- sd(result4[[4]][(burnin+1):(burnin+nsamp)])
# 
#   samples <- matrix(result4[[3]], ncol=NDIM, byrow=TRUE)
# 
#   # Process stimuli
#   stims <- vector("list", nsamp-1)
#   for(j in 1:(nsamp-1)) {
#     stims[[j]] <- matrix(samples[j,1:(ncolX*NS)], ncol=NS, byrow=TRUE)
#   }
#   stim.array <- array(as.numeric(unlist(stims)), dim=c(ncolX, NS, (nsamp-1)))
#   stim.mean <- aaply(stim.array, c(1,2), mean, na.rm=TRUE)
#   stim.lower <- aaply(stim.array, c(1,2), quantile, ll, na.rm=TRUE)
#   stim.upper <- aaply(stim.array, c(1,2), quantile, ul, na.rm=TRUE)
# 
#   stim.samples <- matrix(c(stim.array), ncol=ncolX*NS, byrow=TRUE)
#   stim.samples <- stim.samples[-nrow(stim.samples), ]
# 
#   # Process individuals
#   individuals <- vector("list", (nsamp-1))
#   for(j in 1:(nsamp-1)) {
#     individuals[[j]] <- matrix(c(samples[j,-(1:(ncolX*NS))], 0), ncol=NS, byrow=TRUE)
#   }
#   indiv.array <- array(as.numeric(unlist(individuals)), dim=c(nrowX, NS, (nsamp-1)))
#   indiv.mean <- aaply(indiv.array, c(1,2), mean, na.rm=TRUE)
#   indiv.lower <- aaply(indiv.array, c(1,2), quantile, ll, na.rm=TRUE)
#   indiv.upper <- aaply(indiv.array, c(1,2), quantile, ul, na.rm=TRUE)
# 
#   indiv.samples <- matrix(c(indiv.array), ncol=nrowX*NS, byrow=TRUE)
# 
#   # Set class and names
#   class(stim.samples) <- "mcmc"
#   class(indiv.samples) <- "mcmc"
#   if(!is.null(colnames(input))) {
#     rownames(stim.mean) <- rownames(stim.lower) <- rownames(stim.upper) <- colnames(input)
#   }
# 
#   # Prepare results
#   orig.res <- list(
#     stim.samples = stim.samples,
#     indiv.samples = indiv.samples,
#     stimuli = list(mean = stim.mean, lower=stim.lower, upper=stim.upper),
#     individuals = list(mean = indiv.mean, lower=indiv.lower, upper=indiv.upper)
#   )
# 
#   # Results with Procrustes Rotation
#   stim.rot <- matrix(NA, nrow=(nsamp-1), ncol=ncolX*2)
#   indiv.rot <- matrix(NA, nrow=(nsamp-1), ncol=length(c(individuals[[1]])))
#   for(i in 1:(nsamp-1)) {
#     p <- procrustes(stim.array[,,i], lbfgs.stimuli, dilation=TRUE, translation=TRUE)
#     stim.rot[i,] <- c(p$X.new)
#     indiv.rot[i,] <- c(with(p, s * individuals[[i]] %*% R + matrix(tt, nrow(individuals[[i]]), ncol(individuals[[i]]), byrow = TRUE)))
#   }
#   stim.mean <- matrix(colMeans(stim.rot, na.rm=TRUE), ncol=dims)
#   stim.lower <- matrix(apply(stim.rot, 2, quantile, ll, na.rm=TRUE), ncol=dims)
#   stim.upper <- matrix(apply(stim.rot, 2, quantile, ul, na.rm=TRUE), ncol=dims)
#   indiv.mean <- matrix(colMeans(indiv.rot), ncol=dims)
#   indiv.lower <- matrix(apply(indiv.rot, 2, quantile, ll, na.rm=TRUE), ncol=dims)
#   indiv.upper <- matrix(apply(indiv.rot, 2, quantile, ul, na.rm=TRUE), ncol=dims)
# 
#   if(!is.null(colnames(input))) {
#     rownames(stim.mean) <- rownames(stim.lower) <- rownames(stim.upper) <- colnames(input)
#   }
#   stim.samples <- stim.rot
#   stim.samples <- stim.samples[-nrow(stim.samples), ]
#   class(stim.samples) <- "mcmc"
# 
#   indiv.samples <- indiv.rot
#   class(indiv.samples) <- "mcmc"
# 
#   rotated.res <- list(
#     stim.samples = stim.samples,
#     indiv.samples = indiv.samples,
#     stimuli = list(mean = stim.mean, lower=stim.lower, upper=stim.upper),
#     individuals = list(mean = indiv.mean, lower=indiv.lower, upper=indiv.upper)
#   )
# 
#   # Create final object
#   BUobject <- list(
#     retained.obs = keep,
#     smacof.result = SMACOF.result,
#     lbfgs.result = list(stimuli=lbfgs.stimuli, individuals=lbfgs.individuals),
#     samples = samples,
#     result4 = result4,
#     sigma_squared_hat = sigma_squared_hat,
#     sigma_squared_hat_sd = sigma_squared_hat_sd,
#     unrotated = orig.res,
#     rotated = rotated.res
#   )
# 
#   class(BUobject) <- "bayesunfold"
#   return(BUobject)
# }
bayesunfold <- function(input, dims = 2, nsamp = 2000, burnin = 1000, cred.level = 0.9,
                        slice.starts = c("lbfgs", "random"), print.lbfgs="console",
                        print.slice="console", ...) {
# # Check if the shared object file 'lbfgs_bu3.so' is loaded
# if (!"lbfgs_bu3.so" %in% getLoadedDLLs()) {
#   message("The shared object file 'src/lbfgs_bu3.so' is not loaded.")
#   message("Please manually load it using the following command:")
#   message("dyn.load('src/lbfgs_bu3.so')")
# } else {
#   message("'src/lbfgs_bu3.so' is already loaded.")
# }

  # Check if the shared object file 'lbfgs_bu3.so' is loaded
  if (!"lbfgs_bu3.so" %in% names(getLoadedDLLs())) {
    stop("Error: The shared object file 'src/lbfgs_bu3.so' is not loaded.\n",
         "Please load it manually using the following command:\n",
         "dyn.load('src/lbfgs_bu3.so')")
  } else {
    message("'lbfgs_bu3.so' is successfully loaded.")
  }
  
  # Input validation
  if (cred.level > 1 | cred.level < 0) {
    stop("cred.level must be between 0 and 1\n")
  }
  tailprob <- (1 - cred.level) / 2
  ll <- tailprob
  ul <- 1 - tailprob
  slice.starts <- match.arg(slice.starts)
  
  # Helper functions
  set_globals <- function(nslice, nburn, nrowX, ncolX, NS, N, NDIM, UNFOLD, NMISSING, X, CONSTRAINTS) {
    .C("copyFromR",
       as.integer(nslice),
       as.integer(nburn),
       as.integer(nrowX),
       as.integer(ncolX),
       as.integer(NS),
       as.integer(N),
       as.integer(NDIM),
       as.integer(UNFOLD),
       as.integer(NMISSING),
       as.double(X),
       as.double(CONSTRAINTS))
  }
  
  do_lbfgs <- function(kpnp, kpnq, yrotate, rmatrix) {
    .C("mainlbfgs",
       as.integer(kpnp),
       as.integer(kpnq),
       as.double(yrotate),
       as.double(rmatrix))
  }
  
  do_logposterior <- function(theta, XCOORDS, sumsquared, SIGMAPRIOR) {
    .C("keithrules",
       as.double(theta),
       as.double(XCOORDS),
       as.double(sumsquared),
       as.double(SIGMAPRIOR))
  }
  
  do_sliceu <- function(theta, thetanow2, theta1000, ssenow, XTRUE, thetaLeft, thetaRight, WW, PP, XCOORDS, SIGMAPRIOR) {
    .C("sliceunfolding",
       as.double(theta),
       as.double(thetanow2),
       as.double(theta1000),
       as.double(ssenow),
       as.double(XTRUE),
       as.double(thetaLeft),
       as.double(thetaRight),
       as.double(WW),
       as.integer(PP),
       as.double(XCOORDS),
       as.double(SIGMAPRIOR))
  }
  
  # Data preprocessing
  T <- as.matrix(input)
  T[T < 0 | T > 100] <- NA
  nstimuli <- ncol(T)
  cutoff <- 7
  keep <- rowSums(!is.na(T)) >= cutoff
  T <- T[keep,]
  T <- (100-T)/50
  
  # Setup parameters
  nrowX <- nrow(T)
  ncolX <- ncol(T)
  nburn <- burnin
  nslice <- nsamp
  NS <- dims
  N <- NS*(nrowX+ncolX) - ((NS*(NS+1))/2)
  NDIM <- NS*(nrowX+ncolX) - (NS-1)
  UNFOLD <- 1
  NMISSING <- 7
  
  TT <- T
  TT[which(is.na(TT), arr.ind=TRUE)] <- -999.0
  X <- as.vector(t(TT))
  
  # Set up constraints
  CONSTRAINTS <- rep(1, NS*(nrowX+ncolX))
  if (NS==1) {
    CONSTRAINTS[NS*(nrowX+ncolX)] <- 0
  } else if (NS==2) {
    CONSTRAINTS[(NS*(nrowX+ncolX)-NS):(NS*(nrowX+ncolX))] <- 0
  } else if (NS==3) {
    CONSTRAINTS[(NS*(nrowX+ncolX)-4):(NS*(nrowX+ncolX))] <- 0
    CONSTRAINTS[(NS*(nrowX+ncolX)-6)] <- 0
  }
  
  # SMACOF (metric unfolding)
  weights <- matrix(1, nrow=nrowX, ncol=ncolX)
  weights[is.na(T)] <- 0
  SMACOF.result <- smacofRect(T, ndim=NS, circle = "none", weightmat=weights, itmax=10000, ...)
  
  # Prepare for L-BFGS
  zmetric <- as.numeric(t(SMACOF.result$conf.col))
  xmetric <- as.numeric(t(SMACOF.result$conf.row))
  rmatrix <- c(zmetric, xmetric)
  rmatrix[(NS*(nrowX+ncolX)-NS):(NS*(nrowX+ncolX))] <- 0
  yrotate <- rep(0, (NS*(nrowX+ncolX)))
  
  # L-BFGS optimization
  set_globals(nslice, nburn, nrowX, ncolX, NS, N, NDIM, UNFOLD, NMISSING, X, CONSTRAINTS)
  if(print.lbfgs != "console") {
    sink(print.lbfgs)
    lbfgs.result <- do_lbfgs(nrowX, ncolX, yrotate, rmatrix)
    sink()
  } else {
    lbfgs.result <- do_lbfgs(nrowX, ncolX, yrotate, rmatrix)
  }
  
  lbfgs.coords <- lbfgs.result[[3]]
  dim(lbfgs.coords) <- c(NS, (nrowX+ncolX))
  X3 <- t(lbfgs.coords)
  lbfgs.stimuli <- X3[1:ncolX,]
  lbfgs.individuals <- X3[(ncolX+1):(nrowX+ncolX),]
  
  # Bayesian Unfolding (Slice Sampling)
  if(slice.starts == "random") {
    theta <- runif(NDIM, min=-.5, max=.5) # Random starts for Slice Sampler
  } else if(slice.starts == "lbfgs") {
    theta <- lbfgs.result[[3]] # Non-random starts for Slice Sampler
  }
  
  theta2 <- theta
  theta1000 <- rep(0, nslice*NDIM)
  dim(theta1000) <- c(nslice*NDIM, 1)
  ssenow <- rep(0, (2*(nslice+nburn)))
  dim(ssenow) <- c((2*(nslice+nburn)), 1)
  XTRUE <- lbfgs.result[[3]]
  thetaL <- rep(-99.0, NDIM)
  thetaR <- rep(99.0, NDIM)
  dim(thetaL) <- dim(thetaR) <- c(NDIM, 1)
  thetaL[NDIM] <- 0.10
  thetaR[NDIM] <- 0.50
  WW <- 1.0
  PP <- 3.0
  XCOORDS <- rep(0, (nrowX+ncolX)*NS)
  SIGMAPRIOR <- 100.0
  
  if(print.slice != "console") {
    sink(print.slice)
    result4 <- do_sliceu(theta, theta2, theta1000, ssenow, XTRUE, thetaL, thetaR, WW, PP, XCOORDS, SIGMAPRIOR)
    sink()
  } else {
    result4 <- do_sliceu(theta, theta2, theta1000, ssenow, XTRUE, thetaL, thetaR, WW, PP, XCOORDS, SIGMAPRIOR)
  }
  
  # Post-processing
  sigma_squared_hat <- mean(result4[[4]][(burnin+1):(burnin+nsamp)])
  sigma_squared_hat_sd <- sd(result4[[4]][(burnin+1):(burnin+nsamp)])
  
  samples <- matrix(result4[[3]], ncol=NDIM, byrow=TRUE)
  
  # Process stimuli
  stims <- vector("list", nsamp-1)
  for(j in 1:(nsamp-1)) {
    stims[[j]] <- matrix(samples[j,1:(ncolX*NS)], ncol=NS, byrow=TRUE)
  }
  stim.array <- array(as.numeric(unlist(stims)), dim=c(ncolX, NS, (nsamp-1)))
  stim.mean <- aaply(stim.array, c(1,2), mean, na.rm=TRUE)
  stim.lower <- aaply(stim.array, c(1,2), quantile, ll, na.rm=TRUE)
  stim.upper <- aaply(stim.array, c(1,2), quantile, ul, na.rm=TRUE)
  
  stim.samples <- matrix(c(stim.array), ncol=ncolX*NS, byrow=TRUE)
  stim.samples <- stim.samples[-nrow(stim.samples), ]
  
  # Process individuals
  individuals <- vector("list", (nsamp-1))
  for(j in 1:(nsamp-1)) {
    individuals[[j]] <- matrix(c(samples[j,-(1:(ncolX*NS))], 0), ncol=NS, byrow=TRUE)
  }
  indiv.array <- array(as.numeric(unlist(individuals)), dim=c(nrowX, NS, (nsamp-1)))
  indiv.mean <- aaply(indiv.array, c(1,2), mean, na.rm=TRUE)
  indiv.lower <- aaply(indiv.array, c(1,2), quantile, ll, na.rm=TRUE)
  indiv.upper <- aaply(indiv.array, c(1,2), quantile, ul, na.rm=TRUE)
  
  indiv.samples <- matrix(c(indiv.array), ncol=nrowX*NS, byrow=TRUE)
  
  # Set class and names
  class(stim.samples) <- "mcmc"
  class(indiv.samples) <- "mcmc"
  if(!is.null(colnames(input))) {
    rownames(stim.mean) <- rownames(stim.lower) <- rownames(stim.upper) <- colnames(input)
  }
  
  # Prepare results
  orig.res <- list(
    stim.samples = stim.samples,
    indiv.samples = indiv.samples,
    stimuli = list(mean = stim.mean, lower=stim.lower, upper=stim.upper),
    individuals = list(mean = indiv.mean, lower=indiv.lower, upper=indiv.upper)
  )
  
  # Results with Procrustes Rotation
  stim.rot <- matrix(NA, nrow=(nsamp-1), ncol=ncolX*2)
  indiv.rot <- matrix(NA, nrow=(nsamp-1), ncol=length(c(individuals[[1]])))
  for(i in 1:(nsamp-1)) {
    p <- procrustes(stim.array[,,i], lbfgs.stimuli, dilation=TRUE, translation=TRUE)
    stim.rot[i,] <- c(p$X.new)
    indiv.rot[i,] <- c(with(p, s * individuals[[i]] %*% R + matrix(tt, nrow(individuals[[i]]), ncol(individuals[[i]]), byrow = TRUE)))
  }
  stim.mean <- matrix(colMeans(stim.rot, na.rm=TRUE), ncol=dims)
  stim.lower <- matrix(apply(stim.rot, 2, quantile, ll, na.rm=TRUE), ncol=dims)
  stim.upper <- matrix(apply(stim.rot, 2, quantile, ul, na.rm=TRUE), ncol=dims)
  indiv.mean <- matrix(colMeans(indiv.rot), ncol=dims)
  indiv.lower <- matrix(apply(indiv.rot, 2, quantile, ll, na.rm=TRUE), ncol=dims)
  indiv.upper <- matrix(apply(indiv.rot, 2, quantile, ul, na.rm=TRUE), ncol=dims)
  
  if(!is.null(colnames(input))) {
    rownames(stim.mean) <- rownames(stim.lower) <- rownames(stim.upper) <- colnames(input)
  }
  stim.samples <- stim.rot
  stim.samples <- stim.samples[-nrow(stim.samples), ]
  class(stim.samples) <- "mcmc"
  
  indiv.samples <- indiv.rot
  class(indiv.samples) <- "mcmc"
  
  rotated.res <- list(
    stim.samples = stim.samples,
    indiv.samples = indiv.samples,
    stimuli = list(mean = stim.mean, lower=stim.lower, upper=stim.upper),
    individuals = list(mean = indiv.mean, lower=indiv.lower, upper=indiv.upper)
  )
  
  # Create final object
  BUobject <- list(
    retained.obs = keep,
    smacof.result = SMACOF.result,
    lbfgs.result = list(stimuli=lbfgs.stimuli, individuals=lbfgs.individuals),
    samples = samples,
    result4 = result4,
    sigma_squared_hat = sigma_squared_hat,
    sigma_squared_hat_sd = sigma_squared_hat_sd,
    unrotated = orig.res,
    rotated = rotated.res
  )
  
  class(BUobject) <- "bayesunfold"
  return(BUobject)
}
