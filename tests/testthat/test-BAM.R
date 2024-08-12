# # tests/testthat/test-BAM.R
# 
# library(testthat)
# library(asmcjr)
# 
# test_that("BAM function returns correct structure and data types", {
#   skip_if_not_installed("rjags")  # Skip test if rjags is not installed
# 
#   data(bamdata)  # Ensure bamdata is available
# 
#   # Run the BAM function with default settings
#   result <- BAM(bamdata, polarity = 1, zhatSave = TRUE, abSave = FALSE, resp.idealpts = TRUE, n.sample = 5)
# 
#   # Test if the result is a list
#   expect_type(result, "list")
# 
#   # Check that zhat and zhat.ci are present in the result if zhatSave is TRUE
#   expect_true("zhat" %in% names(result))
#   expect_true("zhat.ci" %in% names(result))
# 
#   # Check that the zhat.ci object has the correct structure
#   expect_s3_class(result$zhat.ci, "data.frame")
#   expect_true(all(c("stimulus", "idealpt", "sd", "lower", "upper") %in% names(result$zhat.ci)))
# 
#   # Check if the resp.summary exists when resp.idealpts is TRUE
#   expect_true("resp.summary" %in% names(result))
#   expect_s3_class(result$resp.summary, "data.frame")
#   expect_true(all(c("lower", "median", "upper") %in% names(result$resp.summary)))
# 
#   # Check that the length of samples corresponds to n.sample
#   expect_equal(length(result$zhat), 100)
# 
#   # Check the contents of zhat
#   expect_type(result$zhat, "list")
#   expect_equal(length(result$zhat[[1]]), ncol(bamdata$stims))
# })
# 
# test_that("BAM function handles different configurations", {
#   skip_if_not_installed("rjags")
# 
#   data(bamdata)
# 
#   # Test with abSave = TRUE and zhatSave = FALSE
#   result_ab <- BAM(bamdata, polarity = 1, zhatSave = FALSE, abSave = TRUE, n.sample = 5)
# 
#   # Check that 'a' and 'b' are present
#   expect_true("a" %in% names(result_ab))
#   expect_true("b" %in% names(result_ab))
# 
#   # Test with both abSave and zhatSave = TRUE
#   result_ab_zhat <- BAM(bamdata, polarity = 1, zhatSave = TRUE, abSave = TRUE, n.sample = 5)
# 
#   # Check that 'a', 'b', 'zhat', and 'zhat.ci' are present
#   expect_true("a" %in% names(result_ab_zhat))
#   expect_true("b" %in% names(result_ab_zhat))
#   expect_true("zhat" %in% names(result_ab_zhat))
#   expect_true("zhat.ci" %in% names(result_ab_zhat))
# })
