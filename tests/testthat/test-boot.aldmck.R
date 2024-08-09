# tests/testthat/test-boot_aldmck.R

library(testthat)
library(boot)
library(asmcjr) 
library(basicspace)

test_that("boot.aldmck generates correct output structure", {
  data(franceEES2009)

  # Run the bootstrap analysis
  boot_result <- boot.aldmck(franceEES2009, polarity = 2, respondent = 1, missing = c(77, 88, 89), verbose = FALSE, boot.args = list(R = 100))

  # Check if the result is a list
  expect_type(boot_result, "list")

  # Check if the list has the correct components
  expect_true("sumstats" %in% names(boot_result))
  expect_true("bootres" %in% names(boot_result))

  # Check if sumstats is a dataframe
  expect_s3_class(boot_result$sumstats, "data.frame")

  # Check if bootres is of class "boot"
  expect_s3_class(boot_result$bootres, "boot")

  # Check the structure of sumstats
  sumstats <- boot_result$sumstats
  expect_true(all(c("stimulus", "idealpt", "sd", "lower", "upper") %in% names(sumstats)))
  expect_type(sumstats$stimulus, "integer")
  expect_type(sumstats$idealpt, "double")
  expect_type(sumstats$sd, "double")
  expect_type(sumstats$lower, "double")
  expect_type(sumstats$upper, "double")

  # Check the dimensions of sumstats
  expect_gt(nrow(sumstats), 0)
  expect_equal(ncol(sumstats), 5)
})

test_that("boot.aldmck handles different boot.args correctly", {
  data(franceEES2009)

  # Run the bootstrap analysis with different R values
  boot_result_50 <- boot.aldmck(franceEES2009, polarity = 2, respondent = 1, missing = c(77, 88, 89), verbose = FALSE, boot.args = list(R = 50))
  boot_result_200 <- boot.aldmck(franceEES2009, polarity = 2, respondent = 1, missing = c(77, 88, 89), verbose = FALSE, boot.args = list(R = 200))

  # Check if the results are lists
  expect_type(boot_result_50, "list")
  expect_type(boot_result_200, "list")

  # Check if bootres is of class "boot"
  expect_s3_class(boot_result_50$bootres, "boot")
  expect_s3_class(boot_result_200$bootres, "boot")

  # Check the number of bootstrap replicates
  expect_equal(boot_result_50$bootres$R, 50)
  expect_equal(boot_result_200$bootres$R, 200)
})
