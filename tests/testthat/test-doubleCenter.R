library(testthat)

# Test for doubleCenter function
test_that("doubleCenter handles square matrices correctly", {
  # Create a simple square matrix
  x <- matrix(c(1, 2, 3, 4), nrow = 2)
  
  # Apply the doubleCenter function
  result <- doubleCenter(x)
  
  # Check that the output is a matrix
  expect_true(is.matrix(result))
  
  # Check the dimensions
  expect_equal(dim(result), dim(x))
  
  # Check that the row and column means are approximately zero
  expect_true(all(abs(rowMeans(result)) < 1e-10))
  expect_true(all(abs(colMeans(result)) < 1e-10))
})
