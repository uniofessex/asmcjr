library(testthat)

# Test for doubleCenterRect function
test_that("doubleCenterRect handles rectangular matrices correctly", {
  # Create a simple rectangular matrix
  x <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  
  # Apply the doubleCenterRect function
  result <- doubleCenterRect(x)
  
  # Check that the output is a matrix
  expect_true(is.matrix(result))
  
  # Check the dimensions
  expect_equal(dim(result), dim(x))
  
  # Check that the row and column means are approximately zero
  expect_true(all(abs(rowMeans(result)) < 1e-10))
  expect_true(all(abs(colMeans(result)) < 1e-10))
})

test_that("doubleCenterRect handles edge cases with one row or one column", {
  # Create a matrix with one row
  x <- matrix(c(1, 2, 3), nrow = 1)
  
  # Apply the doubleCenterRect function
  result <- doubleCenterRect(x)
  
  # Check that the output is a matrix and has the correct dimensions
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(x))
  
  # Create a matrix with one column
  y <- matrix(c(1, 2, 3), ncol = 1)
  
  # Apply the doubleCenterRect function
  result <- doubleCenterRect(y)
  
  # Check that the output is a matrix and has the correct dimensions
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(y))
})