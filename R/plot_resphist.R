#' Plot Respondent Ideal Points Histogram
#'
#' This function generates a histogram plot of respondent ideal points, optionally grouped or faceted by variables.
#'
#' @param result An object of class "aldmck" or "blackbox" containing respondent data.
#' @param groupVar Optional grouping variable.
#' @param facetVar Optional faceting variable.
#' @param addStim Logical, whether to add stimuli points to the plot.
#' @param scaleDensity Logical, whether to scale the density by group proportions.
#' @param weights Character, specifying which weights to use ("all", "positive", "negative").
#' @param xlab Character, label for the x-axis.
#' @param main Character, main title for the plot.
#' @param ylab Character, label for the y-axis.
#' @param whichRes Optional, which respondent data to use.
#' @param dim Optional, dimension for "blackbox" objects.
#' @param ... Additional arguments passed to ggplot2 functions.
#' @return A ggplot2 object
#' @importFrom MASS polr
#' @import ggplot2
#' @examples
#' \dontrun{
#' data(result.france)
#' plot_resphist(result.france, addStim = TRUE, xlab = "Left-Right") +
#'   theme(legend.position = "bottom", aspect.ratio = 1) +
#'   guides(shape = guide_legend(override.aes = list(size = 4), nrow = 3)) +
#'   labs(shape = "Party", colour = "Party")
#' }
#' @export
# plot_resphist <- function(result, groupVar = NULL, facetVar=NULL, addStim = FALSE,
#                           scaleDensity=TRUE, weights = c("all", "positive", "negative"),
#                           xlab = NULL, main = NULL, ylab = NULL, whichRes = NULL,
#                           dim = NULL, ...){
#   w <- match.arg(weights)
#   shapes <- c(15,16,18, 24, 25, 0,1,2,3,4,5,6,7)
#   if(class(result) == "aldmck"){
#     v <- result$respondents
#   }
#   if(class(result) == "blackbox"){
#     if(is.null(dim)){stop("For blackbox, 'dim' must be specified\n")}
#     if(is.null(whichRes)){wres <- dim}
#     else{wres <- whichRes}
#     v <- data.frame("idealpt" = result$individuals[[wres]][,dim], "weight" = 1)
#   }
#   gv <- NULL
#   if(!is.null(groupVar)){
#     v$stimulus <- groupVar
#     gv <- c(gv, "stimulus")
#   }
#   if(!is.null(facetVar)){
#     v$facet <- facetVar
#     gv <- c(gv, "facet")
#   }
#   v <- na.omit(v)
#   xl <- ifelse(is.null(xlab), "Ideal Points", xlab)
#   yl <- ifelse(is.null(ylab), "Density", ylab)
#   main <- ifelse(is.null(main), "", main)
#   if(is.null(groupVar)& is.null(facetVar)){
#     if(w == "all"){
#       g <- ggplot(v, aes(x=idealpt)) +  xlab(xl) + ylab(yl) + ggtitle(main) +stat_density(geom="line") + theme_bw()
#     }
#     if(w == "positive"){
#       posv <- v[which(v$weight > 0),]
#       xl <- paste0(xl, " (n=", nrow(posv), ")")
#       g <- ggplot(posv, aes(x=idealpt)) +  xlab(xl) + ylab(yl) + ggtitle(main) + stat_density(geom="line")  + theme_bw()
#     }
#     if(w == "negative"){
#       negv <- v[which(v$weight < 0),]
#       xl <- paste0(xl, " (n=", nrow(negv), ")")
#       g <- ggplot(negv, aes(x=idealpt)) +  xlab(xl) + ylab(yl) + ggtitle(main)+ stat_density(geom="line")  + theme_bw()
#     }
#   }
#   else{
#     v$newg <- apply(v[,gv,drop=FALSE], 1, paste, collapse=";")
#     ng <- length(table(v$newg))
#     if(w == "all"){
#       props <- table(v$newg)/sum(table(v$newg))
#       bd <- by(v$idealpt, list(v$newg), density)
#       lens <- sapply(bd, function(z)length(z$x))
#       w0 <- which(lens == 0)
#       if(length(w0) > 0){
#         for(j in length(w0):1){bd[[w0[j]]] <- NULL}
#       }
#       for(i in 1:length(bd)){
#         if(scaleDensity)bd[[i]]$y <- bd[[i]]$y*props[i]
#         bd[[i]]$newg <- names(bd)[i]
#       }
#       bd <- lapply(bd, function(z)data.frame("idealpt" = z$x, "Density"=z$y, "newg"=z$newg))
#       bd <- do.call(rbind, bd)
#       t1 <- do.call(rbind, strsplit(as.character(bd$newg), split=";", fixed=T))
#       t1 <- as.data.frame(t1)
#       names(t1) <- c(gv)
#       bd <- cbind(bd, t1)
#       g <- ggplot(bd, aes(x=idealpt, y=Density, group=stimulus, color=stimulus)) + geom_line() + scale_color_manual(values=gray.palette(ng)) +
#         xlab(xl) + ylab(yl) + ggtitle(main) + theme_bw()
#     }
#     if(w == "positive"){
#       posv <- v[which(v$weight > 0),]
#       xl <- paste0(xl, " (n=", nrow(posv), ")")
#       props <- table(posv$stimulus)/sum(table(posv$stimulus))
#       bd <- by(posv$idealpt, list(posv$stimulus), density)
#       lens <- sapply(bd, function(z)length(z$x))
#       w0 <- which(lens == 0)
#       if(length(w0) > 0){
#         for(j in length(w0):1){bd[[w0[j]]] <- NULL}
#       }
#       for(i in 1:length(bd)){
#         if(scaleDensity)bd[[i]]$y <- bd[[i]]$y*props[i]
#         bd[[i]]$stimulus <- factor(i, levels=1:length(bd), labels=names(bd))
#       }
#       bd <- lapply(bd, function(z)data.frame("idealpt" = z$x, "Density"=z$y, "stimulus"=z$stimulus))
#       bd <- do.call(rbind, bd)
#       g <- ggplot(bd, aes(x=idealpt, y=Density, group=stimulus, color=stimulus)) + geom_line() + scale_color_manual(values=gray.palette(ng)) +
#         xlab(xl) + ylab(yl) + ggtitle(main) + theme_bw()
#     }
#     if(w == "negative"){
#       negv <- v[which(v$weight < 0),]
#       xl <- paste0(xl, " (n=", nrow(negv), ")")
#       props <- table(negv$stimulus)/sum(table(negv$stimulus))
#       bd <- by(negv$idealpt, list(negv$stimulus), density)
#       lens <- sapply(bd, function(z)length(z$x))
#       w0 <- which(lens == 0)
#       if(length(w0) > 0){
#         for(j in length(w0):1){bd[[w0[j]]] <- NULL}
#       }
#       for(i in 1:length(bd)){
#         if(scaleDensity)bd[[i]]$y <- bd[[i]]$y*props[i]
#         bd[[i]]$stimulus <- factor(i, levels=1:length(bd), labels=names(bd))
#       }
#       bd <- lapply(bd, function(z)data.frame("idealpt" = z$x, "Density"=z$y, "stimulus"=z$stimulus))
#       bd <- do.call(rbind, bd)
#       g <- ggplot(bd, aes(x=idealpt, y=Density, group=stimulus, color=stimulus)) + geom_line() + scale_color_manual(values=gray.palette(ng)) +
#         xlab(xl) + ylab(yl) + ggtitle(main) + theme_bw()
#     }
#   }
#   if(addStim){
#     tmp <- na.omit(result$stimuli)
#     if(!is.null(groupVar)){
#       tmp <- tmp[which(names(tmp) %in% unique(groupVar))]
#       n <- names(tmp)
#       p <- data.frame("idealpt" = tmp, "stimulus" = factor(n, levels=n[order(tmp)]))
#       g <- g + geom_point(data=p, aes(y=0, group=stimulus, pch=stimulus, col=stimulus, size=2.5)) +
#         scale_shape_manual(values=shapes[1:nrow(p)]) + theme_bw() + scale_size(2.5, guide=FALSE)
#     }
#     else{
#       n <- names(tmp)
#       p <- data.frame("idealpt" = tmp, "stimulus" = factor(n, levels=n[order(tmp)]))
#       g <- g + geom_point(data=p, aes(y=0, group=stimulus, pch=stimulus, col=stimulus, size=2.5)) +
#         scale_shape_manual(values=shapes[1:nrow(p)]) + scale_color_manual(values=gray.palette(nrow(p))) +
#         theme_bw() + scale_size(2.5, guide=FALSE)
#     }
#   }
#   if(!is.null(facetVar)){
#     g <- g + facet_wrap(~facet)
#   }
#   return(g)
#   if(exists("bd")){
#     out <- list(data=bd)
#   }
#   else{
#     out <- list(data=v)
#   }
#   if(exists("p")){
#     out$stim_data <- p
#   }
#   return(out)
# }

plot_resphist <- function(result, groupVar = NULL, facetVar=NULL, addStim = FALSE,
                          scaleDensity=TRUE, weights = c("all", "positive", "negative"),
                          xlab = NULL, main = NULL, ylab = NULL, whichRes = NULL,
                          dim = NULL, ...){
  w <- match.arg(weights)
  shapes <- c(15,16,18, 24, 25, 0,1,2,3,4,5,6,7)
  if (inherits(result, "aldmck")) {
    v <- result$respondents
  }
  if (inherits(result, "blackbox")) {
    if (is.null(dim)) { stop("For blackbox, 'dim' must be specified\n") }
    if (is.null(whichRes)) { wres <- dim }
    else { wres <- whichRes }
    v <- data.frame("idealpt" = result$individuals[[wres]][,dim], "weight" = 1)
  }
  gv <- NULL
  if (!is.null(groupVar)) {
    v$stimulus <- groupVar
    gv <- c(gv, "stimulus")
  }
  if (!is.null(facetVar)) {
    v$facet <- facetVar
    gv <- c(gv, "facet")
  }
  v <- na.omit(v)
  xl <- ifelse(is.null(xlab), "Ideal Points", xlab)
  yl <- ifelse(is.null(ylab), "Density", ylab)
  main <- ifelse(is.null(main), "", main)
  if (is.null(groupVar) & is.null(facetVar)) {
    if (w == "all") {
      g <- ggplot(v, aes(x=idealpt)) +  xlab(xl) + ylab(yl) + ggtitle(main) +stat_density(geom="line") + theme_bw()
    }
    if (w == "positive") {
      posv <- v[which(v$weight > 0),]
      xl <- paste0(xl, " (n=", nrow(posv), ")")
      g <- ggplot(posv, aes(x=idealpt)) +  xlab(xl) + ylab(yl) + ggtitle(main) + stat_density(geom="line")  + theme_bw()
    }
    if (w == "negative") {
      negv <- v[which(v$weight < 0),]
      xl <- paste0(xl, " (n=", nrow(negv), ")")
      g <- ggplot(negv, aes(x=idealpt)) +  xlab(xl) + ylab(yl) + ggtitle(main)+ stat_density(geom="line")  + theme_bw()
    }
  } else {
    v$newg <- apply(v[,gv,drop=FALSE], 1, paste, collapse=";")
    ng <- length(table(v$newg))
    if (w == "all") {
      props <- table(v$newg)/sum(table(v$newg))
      bd <- by(v$idealpt, list(v$newg), density)
      lens <- sapply(bd, function(z)length(z$x))
      w0 <- which(lens == 0)
      if (length(w0) > 0) {
        for (j in length(w0):1) { bd[[w0[j]]] <- NULL }
      }
      for (i in 1:length(bd)) {
        if (scaleDensity) bd[[i]]$y <- bd[[i]]$y*props[i]
        bd[[i]]$newg <- names(bd)[i]
      }
      bd <- lapply(bd, function(z) data.frame("idealpt" = z$x, "Density"=z$y, "newg"=z$newg))
      bd <- do.call(rbind, bd)
      t1 <- do.call(rbind, strsplit(as.character(bd$newg), split=";", fixed=T))
      t1 <- as.data.frame(t1)
      names(t1) <- c(gv)
      bd <- cbind(bd, t1)
      g <- ggplot(bd, aes(x=idealpt, y=Density, group=stimulus, color=stimulus)) + geom_line() + scale_color_manual(values=gray.palette(ng)) +
        xlab(xl) + ylab(yl) + ggtitle(main) + theme_bw()
    }
    if (w == "positive") {
      posv <- v[which(v$weight > 0),]
      xl <- paste0(xl, " (n=", nrow(posv), ")")
      props <- table(posv$stimulus)/sum(table(posv$stimulus))
      bd <- by(posv$idealpt, list(posv$stimulus), density)
      lens <- sapply(bd, function(z)length(z$x))
      w0 <- which(lens == 0)
      if (length(w0) > 0) {
        for (j in length(w0):1) { bd[[w0[j]]] <- NULL }
      }
      for (i in 1:length(bd)) {
        if (scaleDensity) bd[[i]]$y <- bd[[i]]$y*props[i]
        bd[[i]]$stimulus <- factor(i, levels=1:length(bd), labels=names(bd))
      }
      bd <- lapply(bd, function(z)data.frame("idealpt" = z$x, "Density"=z$y, "stimulus"=z$stimulus))
      bd <- do.call(rbind, bd)
      g <- ggplot(bd, aes(x=idealpt, y=Density, group=stimulus, color=stimulus)) + geom_line() + scale_color_manual(values=gray.palette(ng)) +
        xlab(xl) + ylab(yl) + ggtitle(main) + theme_bw()
    }
    if (w == "negative") {
      negv <- v[which(v$weight < 0),]
      xl <- paste0(xl, " (n=", nrow(negv), ")")
      props <- table(negv$stimulus)/sum(table(negv$stimulus))
      bd <- by(negv$idealpt, list(negv$stimulus), density)
      lens <- sapply(bd, function(z)length(z$x))
      w0 <- which(lens == 0)
      if (length(w0) > 0) {
        for (j in length(w0):1) { bd[[w0[j]]] <- NULL }
      }
      for (i in 1:length(bd)) {
        if (scaleDensity) bd[[i]]$y <- bd[[i]]$y*props[i]
        bd[[i]]$stimulus <- factor(i, levels=1:length(bd), labels=names(bd))
      }
      bd <- lapply(bd, function(z)data.frame("idealpt" = z$x, "Density"=z$y, "stimulus"=z$stimulus))
      bd <- do.call(rbind, bd)
      g <- ggplot(bd, aes(x=idealpt, y=Density, group=stimulus, color=stimulus)) + geom_line() + scale_color_manual(values=gray.palette(ng)) +
        xlab(xl) + ylab(yl) + ggtitle(main) + theme_bw()
    }
  }
  if (addStim) {
    tmp <- na.omit(result$stimuli)
    if (!is.null(groupVar)) {
      tmp <- tmp[which(names(tmp) %in% unique(groupVar))]
      n <- names(tmp)
      p <- data.frame("idealpt" = tmp, "stimulus" = factor(n, levels=n[order(tmp)]))
      g <- g + geom_point(data=p, aes(y=0, group=stimulus, pch=stimulus, col=stimulus, size=2.5)) +
        scale_shape_manual(values=shapes[1:nrow(p)]) + theme_bw() + scale_size(2.5, guide=FALSE)
    } else {
      n <- names(tmp)
      p <- data.frame("idealpt" = tmp, "stimulus" = factor(n, levels=n[order(tmp)]))
      g <- g + geom_point(data=p, aes(y=0, group=stimulus, pch=stimulus, col=stimulus, size=2.5)) +
        scale_shape_manual(values=shapes[1:nrow(p)]) + scale_color_manual(values=gray.palette(nrow(p))) +
        theme_bw() + scale_size(2.5, guide=FALSE)
    }
  }
  if (!is.null(facetVar)) {
    g <- g + facet_wrap(~facet)
  }
  return(g)
  if (exists("bd")) {
    out <- list(data=bd)
  } else {
    out <- list(data=v)
  }
}
  