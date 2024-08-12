#' Bayesian Multidimensional Scaling (BMDS)
#'
#' This function performs Bayesian Multidimensional Scaling (BMDS) using a specified distance matrix and predefined positions for stimuli. It uses the JAGS (Just Another Gibbs Sampler) framework to estimate the positions of stimuli in a lower-dimensional space.
#'
#' @param data A numeric matrix representing the observed distance matrix. The matrix should be symmetric with non-negative entries.
#' @param posStims A numeric vector of length 2 specifying the indices of the stimuli that should be constrained to have positive coordinates on the first and second dimensions, respectively.
#' @param negStims A numeric vector of length 2 specifying the indices of the stimuli that should be constrained to have negative coordinates on the first and second dimensions, respectively.
#' @param z A numeric matrix with the initial positions for the stimuli. Missing values (NA) indicate positions to be estimated.
#' @param fname A character string specifying the file name to which the JAGS model code should be written. This argument is required.
#' @param n.sample An integer specifying the number of MCMC samples to draw. The default is 2500.
#' @param ... Additional arguments passed to the JAGS model, such as `n.chains` and `n.adapt`.
#'
#' @return An invisible list with two components:
#' \item{zhat}{An object containing the MCMC samples for the estimated positions.}
#' \item{zhat.ci}{A data frame containing the summary statistics (mean, standard deviation, and credible intervals) for the estimated positions.}
#'
#' @details This function requires the `rjags` package and JAGS software to run. It constructs a JAGS model that imposes constraints on the positions of certain stimuli (specified by `posStims` and `negStims`) and estimates the remaining positions using MCMC sampling. The model is written to a file specified by `fname`, and the function then runs the JAGS model using the specified number of samples.
#' @importFrom stats runif
#' @examples
#' \dontrun{
#' # Load the nations dataset (assuming it's already loaded in your environment)
#' # If not, load or create a similar dataset with distance values.
#' data(nations)
#' # Initialize the z matrix with NA values and set specific starting points
#' z <- matrix(NA, nrow = nrow(nations), ncol = 2)
#' z[10, ] <- 0        # Fix both dimensions of the 10th stimulus to 0
#' z[11, 2] <- 0       # Fix the second dimension of the 11th stimulus to 0
#' 
#' # Perform Bayesian Multidimensional Scaling (BMDS)
#' nations.bmds <- BMDS(nations, 
#'                      posStims = c(7, 2),  # Stimuli to constrain to positive coordinates
#'                      negStims = c(9, 8),  # Stimuli to constrain to negative coordinates
#'                      z = z,               # Initial positions matrix
#'                      fname = "nations_jags.txt",  # File to write the JAGS model
#'                      n.sample = 10000,     # Number of MCMC samples
#'                      n.adapt = 10000)      # Number of adaptation steps
#' 
#' # Inspect the result
#' print(nations.bmds$zhat.ci)  # Print the summary statistics for the estimated positions
#' }
#' @seealso \code{\link[rjags]{jags.model}}, \code{\link[rjags]{coda.samples}}, \code{\link[coda]{summary.mcmc.list}}
#'
#' @export
BMDS <- function(data, posStims, negStims, z, fname=NULL, n.sample = 2500, ...){
  if(requireNamespace("rjags")){  
    jmod <- get("jags.model", asNamespace("rjags"))
    jcs <- get("coda.samples", asNamespace("rjags"))
    jsum <- get("summary.mcarray", asNamespace("rjags"))
    csum <- get("summary.mcmc.list", asNamespace("coda"))
    sumsamp <- function(x){
      if(!is.null(dim(x))){
        jsum(x)  
      }
      else{
        csum(x)
      }
    }
    args <- as.list(match.call(expand.dots = FALSE)$`...`)
    if(!("n.chains" %in% names(args)))args$n.chains = 2
    if(!("n.adapt" %in% names(args)))args$n.adapt = 10000
    if(!("inits" %in% names(args))){
      z.init <- array(runif(2*nrow(z), -5, 5), dim=dim(z))
      z.init[which(!is.na(z))] <- NA
      z.init[posStims[1], 1] <- abs(z.init[posStims[1], 1] )
      z.init[posStims[2], 2] <- abs(z.init[posStims[2], 2] )
      z.init[negStims[1], 1] <- -abs(z.init[negStims[1], 1] )
      z.init[negStims[2], 2] <- -abs(z.init[negStims[2], 2] )
      args$inits <- function(){list(z=z.init)}
    }
    
    m1 <- "model{
        for (i in 1:(N-1)){
            for (j in (i+1):N){
                dstar[i,j] ~ dlnorm(mu[i,j],tau)
                mu[i,j] <- log(sqrt((z[i,1]-z[j,1])*(z[i,1]-z[j,1])+(z[i,2]-z[j,2])*(z[i,2]-z[j,2])))
            }
        }
        tau ~ dunif(0,10)"
    nZ <- nrow(data)
    d1const <- c(posStims[1], negStims[1])
    d2const <- c(posStims[2], negStims[2])
    d1const2 <- d2const2 <- c("T(0, )", "T(,0)")
    d1const2 <- d1const2[order(d1const)]
    d1const <- sort(d1const)
    d2const2 <- d2const2[order(d2const)]
    d2const <- sort(d2const)
    for(i in 1:nZ){
      if(is.na(z[i,1])){
        m1 <- paste(m1, "\nz[", i,", 1] ~ dnorm(0,.01)", ifelse(i %in% d1const, d1const2[which(d1const == i)], ""), sep="")
      }
      if(is.na(z[i,2])){
        m1 <- paste(m1, "\nz[", i,", 2] ~ dnorm(0,.01)", ifelse(i %in% d2const, d2const2[which(d2const == i)], ""), sep="")
      }
    }
    m1 <- paste(m1, "\n}",sep="")
    if(is.null(fname))stop("Must specify a file name to write the code to")
    cat(m1, file=fname)
    data <- as.matrix(data)
    args$file <- fname
    args$data <- list('N'=nrow(data), dstar = as.matrix(max(data)-data), z=z)
    mod.sim <- do.call("jmod", args)
    samples <- jcs(mod.sim,'z',  n.sample, thin=1)
    zhat <- samples
    zhat.sum <- sumsamp(zhat)
    zhat.ci <- data.frame("stimulus" = c(outer(colnames(data), c(" D1", " D2"), paste0)),
                          "idealpt" = zhat.sum$statistics[,1],
                          "sd" = zhat.sum$statistics[,2],
                          "lower" = zhat.sum$quantiles[,1],
                          "upper" = zhat.sum$quantiles[,5])
    rownames(zhat.ci) <- NULL
    res.list = list(zhat=zhat, zhat.ci = zhat.ci)
    invisible(res.list)
  } else{
    stop("You must install JAGS (https://sourceforge.net/projects/mcmc-jags/) and the rjags package for this function to work.\n")
  }
}
