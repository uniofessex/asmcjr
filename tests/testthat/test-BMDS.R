# library(testthat)
# library(rjags)
# 
# test_that("BMDS handles missing file name", {
#   # Create a small mock dataset
#   nations <- matrix(runif(100, 0, 10), nrow = 10)
#   posStims <- c(1, 2)
#   negStims <- c(3, 4)
#   z <- matrix(NA, nrow = nrow(nations), ncol = 2)
#   
#   # Expect an error because fname is NULL
#   expect_error(BMDS(nations, posStims, negStims, z = z, fname = NULL, n.sample = 100),
#                "Must specify a file name to write the code to")
# })
