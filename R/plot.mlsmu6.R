#' Plot Method for `mlsmu6` Class Objects
#'
#' This function creates a 2D plot of the results from the `mlsmu6` function. It visualizes the estimated positions of both individuals and stimuli in a lower-dimensional space.
#'
#' @param x An object of class `mlsmu6`, typically the output from the `mlsmu6` function.
#' @param ... Additional arguments passed to the plot method (currently not used).
#' @param selected.stims A vector of stimulus names to highlight in the plot. If `NULL`, all stimuli are plotted.
#' @param ind.id.size Numeric value specifying the text size for individual labels in the plot. Default is 3.
#' @param stim.id.size Numeric value specifying the text size for stimulus labels in the plot. Default is 6.
#'
#' @return A `ggplot` object representing the 2D plot of the estimated positions.
#'
#' @details The plot shows individuals and stimuli in a 2D space defined by the first two dimensions (`Dim1` and `Dim2`). Individuals are represented with smaller labels, while stimuli are larger, making them more prominent. If `selected.stims` is specified, only the selected stimuli will be highlighted, while the rest will be omitted from the plot.
#' @import ggplot2
#' @examples
#' \dontrun{
#' # Assuming `result` is the output from `mlsmu6` function:
#' data(mlsmu6_out)
#' plot.mlsmu6(mlsmu6_out) 
#' }
#'
#' @export
plot.mlsmu6 <- function(x, ..., selected.stims = NULL, ind.id.size = 3, stim.id.size = 6) {
  dfi <- x$inds
  if (is.null(dfi$id)) {
    dfi$id <- "x"
  }
  dfi$id <- as.character(dfi$id)
  dfii <- sort(unique(dfi$id))
  dfs <- x$stims
  dfs$id <- rownames(dfs)
  if (!is.null(selected.stims)) {
    dfs <- dfs[rownames(dfs) %in% selected.stims, ]
  }
  df <- rbind(dfi, dfs)
  levs1 <- dfii
  levs2 <- dfs$id
  df$id <- factor(df$id, levels = c(sort(levs1), sort(levs2)))
  
  g <- ggplot(df, aes(x = .data$Dim1, y = .data$Dim2, label = .data$id, colour = .data$id, size = .data$id)) +
    geom_text(show.legend = FALSE) +
    scale_colour_manual(values = c("gray25", "gray50", rep("black", nrow(dfs))), guide = "none") +
    scale_size_manual(values = c(rep(ind.id.size, length(dfii)), rep(stim.id.size, nrow(dfs)))) +
    scale_shape_discrete(name = "ID", breaks = dfii, labels = dfii) +
    theme_bw() +
    theme(aspect.ratio = 1)
  return(g)
}

# plot.mlsmu6 <- function(x, ..., selected.stims = NULL, ind.id.size = 3, stim.id.size = 6) {
#   dfi <- x$inds
#   if (is.null(dfi$id)) {
#     dfi$id <- "x"
#   }
#   dfi$id <- as.character(dfi$id)
#   dfii <- sort(unique(dfi$id))
#   dfs <- x$stims
#   dfs$id <- rownames(dfs)
#   if (!is.null(selected.stims)) {
#     dfs <- dfs[which(rownames(dfs) %in% selected.stims), ]
#   }
#   df <- rbind(dfi, dfs)
#   levs <- unique(df$id)
#   levs1 <- dfii
#   levs2 <- dfs$id
#   df$id <- factor(df$id, c(sort(levs1), sort(levs2)))
#   g <- ggplot(df, aes_string(x = "Dim1", y = "Dim2", label = "id", colour = "id", size = "id")) +
#     geom_text(show.legend = FALSE) +
#     scale_colour_manual(values = c("gray25", "gray50", rep("black", nrow(dfs))), guide = "none") +
#     scale_size_manual(values = c(rep(ind.id.size, length(dfii)), rep(stim.id.size, nrow(dfs)))) +
#     scale_shape_discrete(name = "ID", breaks = dfii, labels = dfii) +
#     theme_bw() + theme(aspect.ratio = 1)
#   return(g)
# }


