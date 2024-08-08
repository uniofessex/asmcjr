# tests/testthat/test-plot_resphist.R

library(testthat)
library(ggplot2)
library(asmcjr) 

test_that("plot_resphist generates a ggplot object", {
  data(result.france)
  
  # Basic test
  p <- plot_resphist(result.france, addStim = TRUE, xlab = "Left-Right")
  expect_s3_class(p, "ggplot")
  
  # Test with theme and guides modifications
  p <- p + theme(legend.position = "bottom", aspect.ratio = 1) +
    guides(shape = guide_legend(override.aes = list(size = 4), nrow = 3)) +
    labs(shape = "Party", colour = "Party")
  expect_s3_class(p, "ggplot")
  
  # Test for correct labels
  expect_true("Left-Right" %in% p$labels$x)
  expect_true("Density" %in% p$labels$y)
  expect_true("Party" %in% p$labels$shape)
  expect_true("Party" %in% p$labels$colour)
  
  # Additional tests can be added here for further verification
})
