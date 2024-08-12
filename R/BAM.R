#' Bayesian Aldrich-McKelvey (BAM) Scaling
#'
#' This function performs Bayesian Aldrich-McKelvey (BAM) scaling using JAGS to estimate ideal points and optionally save intermediate parameters.
#'
#' @param data A dataset prepared using the `bamPrep` function.
#' @param polarity Integer specifying the column number of the stimulus used to set the polarity.
#' @param zhatSave Logical, whether to save the posterior samples of `zhat` (default is TRUE).
#' @param abSave Logical, whether to save the posterior samples of `a` and `b` parameters (default is FALSE).
#' @param resp.idealpts Logical, whether to compute and save respondent ideal points (default is FALSE).
#' @param n.sample Integer, the number of MCMC samples to draw (default is 2500).
#' @param ... Additional arguments passed to `jags.model`.
#' @return A list containing posterior samples and summaries for the estimated parameters.
#' @importFrom rjags jags.model 
#' @importFrom rjags coda.samples
#' @importFrom basicspace aldmck
#' @importFrom stats na.omit rnorm quantile
#' @examples
#' \dontrun{
#' data(bamdata)
#' bam.france <- BAM(bamdata, polarity=2, n.adapt=2500, n.sample=5000,
#'                   zhat=TRUE, ab=TRUE, resp.idealpts=TRUE)
#' }
#' @export

BAM <- function(data, polarity, zhatSave = TRUE, abSave = FALSE, resp.idealpts = FALSE, n.sample = 2500, ...) {
  if (requireNamespace("rjags")) {
    # Internal JAGS and coda functions
    jmod <- get("jags.model", asNamespace("rjags"))
    jcs <- get("coda.samples", asNamespace("rjags"))
    jsum <- get("summary.mcarray", asNamespace("rjags"))
    csum <- get("summary.mcmc.list", asNamespace("coda"))

    sumsamp <- function(x) {
      if (!is.null(dim(x))) {
        jsum(x)
      } else {
        csum(x)
      }
    }

    if (!("bamPrep" %in% class(data))) stop("Data should be output from the bamPrep function")
    args <- as.list(match.call(expand.dots = FALSE)$`...`)
    if (!("n.chains" %in% names(args))) args$n.chains = 2
    if (!("n.adapt" %in% names(args))) args$n.adapt = 10000
    if (!("inits" %in% names(args))) {
      orig <- aldmck(na.omit(data$stims), respondent = 0, polarity = polarity, verbose = FALSE)
      args$inits <- list()
      for (i in 1:args$n.chains) {
        zhs <- orig$stimuli + rnorm(length(orig$stimuli), 0, 1)
        zhs[polarity] <- -abs(zhs[polarity])
        args$inits[[i]] <- list(zhatstar = zhs)
      }
    }

    # Model specification
    args$file <- system.file("templates/BAM_JAGScode.bug", package = "asmcjr")
    lower <- rep(-100, ncol(data$stims))
    upper <- rep(100, ncol(data$stims))
    upper[polarity] <- 0
    args$data <- list('z' = data$stims, q = ncol(data$stims), N = nrow(data$stims), lower = lower, upper = upper)
    mod.sim <- do.call("jmod", args)

    # Save zhat or ab samples or both
    if (zhatSave & !abSave) {
      samples <- jcs(mod.sim, 'zhat', n.sample, thin = 1)
      zhat <- samples
      for (i in 1:length(zhat)) colnames(zhat[[i]]) <- colnames(data$stims)
      zhat.sum <- sumsamp(zhat)
      zhat.ci <- data.frame(
        "stimulus" = factor(colnames(data$stims), levels = colnames(data$stims)[order(zhat.sum$statistics)]),
        "idealpt" = zhat.sum$statistics[, 1],
        "sd" = zhat.sum$statistics[, 2],
        "lower" = zhat.sum$quantiles[, 1],
        "upper" = zhat.sum$quantiles[, 5]
      )
      rownames(zhat.ci) <- NULL
      class(zhat.ci) <- c("aldmck_ci", "data.frame")
      res.list = list(zhat = zhat, zhat.ci = zhat.ci)
    }

    if (abSave & !zhatSave) {
      samples <- jcs(mod.sim, c('a', 'b'), n.sample, thin = 1)
      a <- samples[, grep("^a", colnames(samples[[1]]))]
      b <- samples[, grep("^b", colnames(samples[[1]]))]
      res.list = list(a = a, b = b)
    }

    if (abSave & zhatSave) {
      samples <- jcs(mod.sim, c('zhat', 'a', 'b'), n.sample, thin = 1)
      zhat <- samples[, grep("^z", colnames(samples[[1]]))]
      for (i in 1:length(zhat)) colnames(zhat[[i]]) <- colnames(data$stims)
      zhat.sum <- sumsamp(zhat)
      zhat.ci <- data.frame(
        "stimulus" = factor(colnames(data$stims), levels = colnames(data$stims)[order(zhat.sum$statistics)]),
        "idealpt" = zhat.sum$statistics[, 1],
        "sd" = zhat.sum$statistics[, 2],
        "lower" = zhat.sum$quantiles[, 1],
        "upper" = zhat.sum$quantiles[, 5]
      )
      rownames(zhat.ci) <- NULL
      class(zhat.ci) <- c("aldmck_ci", "data.frame")
      a <- samples[, grep("^a", colnames(samples[[1]]))]
      b <- samples[, grep("^b", colnames(samples[[1]]))]
      res.list = list(zhat = zhat, zhat.ci = zhat.ci, a = a, b = b)
    }

    if (resp.idealpts) {
      amat <- do.call(rbind, res.list$a)
      bmat <- do.call(rbind, res.list$b)
      diffs <- t(apply(amat, 1, function(x) data$self - x))
      resp.ideals <- diffs / bmat
      resp.ideal.summary <- t(apply(resp.ideals, 2, quantile, c(.025, .5, .975), na.rm = TRUE))
      resp.ideal.summary <- as.data.frame(resp.ideal.summary)
      names(resp.ideal.summary) <- c("lower", "median", "upper")
      res.list$resp.samples = resp.ideals
      res.list$resp.summary = resp.ideal.summary
    }

    invisible(res.list)
  } else {
    stop("You must install JAGS (https://sourceforge.net/projects/mcmc-jags/) and the rjags package to use this function.\n")
  }
}
