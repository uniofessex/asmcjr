library(testthat)
library(ggplot2)
library(car)

test_that("plot_blackbox handles missing groupVar correctly", {
  data(CDS2000)
  party <- recode(CDS2000[,1], "1='Democrat'; 2='Republican'; else=NA", as.factor=TRUE)
  g <-  plot_blackbox(result.repdem, dims=c(1, 2), groupVar=party, 
                      xlab="First Dimension\n(Left-Right)",
                      ylab="Second Dimension")
  expect_s3_class(g, "ggplot")
  expect_true("group" %in% names(g$data))
})


