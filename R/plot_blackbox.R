#' Plot Blackbox Results with Optional Grouping and Issue Vectors
#'
#' This function creates a scatter plot of the results from a blackbox analysis. It allows for optional grouping of data points and the inclusion of issue vectors based on ordinal regression.
#'
#' @param result A list object containing the results of the blackbox analysis. This object should include individual-level data for each dimension.
#' @param dims A numeric vector of length 2 indicating the dimensions to plot.
#' @param whichRes An optional integer specifying which result set to use from the `result$individuals` list. If `NULL`, the highest dimension is used.
#' @param groupVar An optional vector for grouping the data points. If `NULL`, no grouping is applied.
#' @param issueVector A vector indicating which issues to include as vectors in the plot. This can be a character vector of column names or a numeric vector of column indices.
#' @param data A data frame containing the original data used in the analysis. Required if `issueVector` is specified.
#' @param missing A vector of values to be treated as missing in the data.
#' @param rug Logical, if `TRUE`, adds marginal rugs to the plot.
#' @param xlab A character string for the x-axis label. If `NULL`, the label is generated automatically.
#' @param ylab A character string for the y-axis label. If `NULL`, the label is generated automatically.
#' @param main A character string for the plot title. If `NULL`, no title is added.
#' @param nudgeX A numeric vector for nudging the x positions of the issue vector labels. Defaults to `NULL`.
#' @param nudgeY A numeric vector for nudging the y positions of the issue vector labels. Defaults to `NULL`.
#' @param ... Additional arguments passed to `ggplot2` functions.
#' @return A `ggplot2` object representing the blackbox analysis plot.
#' @import ggplot2
#' @importFrom MASS polr
#' @examples
#' \dontrun{
#' # CDS2000 is a dataset containing U.S. Congressional data,
#' # where the first column represents party affiliation: 1 for Democrats, 2 for Republicans.
#' 
#' # Recode the party variable: 1 = 'Democrat', 2 = 'Republican', others as NA
#' party <- car::recode(CDS2000[,1],
#'                      "1='Democrat'; 2='Republican'; else=NA",
#'                      as.factor=TRUE)
#' 
#' # Create a Blackbox plot to compare the distribution of different parties in a two-dimensional space
#' plot_blackbox(result.repdem, dims=c(1, 2), groupVar=party,
#'               xlab="First Dimension\n(Left-Right)",
#'               ylab="Second Dimension") +
#'   theme(legend.position="bottom", aspect.ratio=1) +
#'   guides(shape=guide_legend(override.aes=list(size=4))) +
#'   labs(colour="Party")
#' }
#' @export
plot_blackbox <- function(result, dims, whichRes=NULL, groupVar=NULL, issueVector=NULL,
                          data=NULL, missing=NULL, rug=FALSE, xlab=NULL, main = NULL, ylab=NULL, nudgeX=NULL, nudgeY=NULL,...){
  wres <- ifelse(is.null(whichRes), max(dims), whichRes)
  dimdat <- result$individuals[[wres]][,dims]
  names(dimdat) <- c("x", "y")
  if(is.null(groupVar)){
    g <- ggplot(dimdat, aes(x=x, y=y)) + geom_point(shape=1, col="gray65") + theme_bw()
  }
  else{
    dimdat$group = groupVar
    dimdat$pch = substr(as.character(dimdat$group), 1, 1)
    ng <- length(unique(na.omit(groupVar)))
    g <- ggplot(dimdat, aes(x=x, y=y)) +
      geom_point(aes(group=group, col=group), alpha=0) +
      geom_text(aes(group=group, col=group, label=pch), show.legend=FALSE) +
      scale_color_manual(values=gray.palette(ng)) +
      guides(colour = guide_legend("Grouping", override.aes = list(size = 2, alpha = 1))) +
      theme_bw()
  }
  if(rug){
    g <- g+geom_rug(show.legend=FALSE)
  }
  if(is.null(xlab)){
    xlab <- paste0("Dimension ", dims[1])
  }
  if(is.null(ylab)){
    ylab <- paste0("Dimension ", dims[2])
  }
  g <- g+ylab(ylab) + xlab(xlab)
  if(!is.null(issueVector)){
    if(is.null(data))stop("If you want to plot issue vectors, you need to specify the data\n")
    if(is.character(issueVector)){
      iss <- match(issueVector, colnames(data))
      if(length(iss) != length(issueVector))stop("At least one of the specified issues didn't match names in the data\n")
    }
    else{
      iss <- issueVector
    }
    dv <- data[,iss, drop=FALSE]
    dv[which(dv %in% missing, arr.ind=TRUE)] <- NA
    dv <- do.call("data.frame", lapply(1:ncol(dv), function(x)as.factor(dv[,x])))
    names(dv) <- colnames(data)[iss]
    op <- list()
    for(i in 1:ncol(dv)){
      op[[i]] <- polr(dv[,i] ~ x + y, data=dimdat, method="probit")
    }
    b <- sapply(op, function(x)x$coef)
    Nvals <- apply(b, 2, function(x)x/sqrt(sum(x^2)))
    scale.fac <- min(apply(dimdat[,1:2], 2, function(x)diff(range(x, na.rm=TRUE)))/2)
    for(i in 1:ncol(Nvals)){
      tmp <- data.frame(x=c(0, scale.fac*Nvals[1,i]), y=c(0, scale.fac*Nvals[2,i]))
      g <- g + geom_line(data=tmp, arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"), size=1.1) +
        geom_line(data=-tmp, arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed"), size=1.1)
    }
    if(is.null(nudgeX)){nudgeX <- rep(0, ncol(Nvals))}
    if(is.null(nudgeY)){nudgeY <- rep(0, ncol(Nvals))}
    colnames(Nvals) <- colnames(data)[iss]
    nvals <- t(-Nvals*scale.fac)
    nvals <- as.data.frame(nvals)
    g <- g + geom_text(data=nvals, aes(x=x, y=y, label=rownames(nvals), size=2),
                       nudge_x = nudgeX, nudge_y=nudgeY, show.legend=FALSE)
  }
  return(g)
}
