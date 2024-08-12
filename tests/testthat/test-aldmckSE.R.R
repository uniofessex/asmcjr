
obj <- list(
  respondents = matrix(c(1, 2, 3, 4, 5, 6), ncol = 2),
  stimuli = c(0.5, 1.5)
)

data <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)

expected_output <- c(2.53311403, 4.66368953)  

test_that("aldmckSE returns correct standard errors", {
  result <- aldmckSE(obj, data)
  expect_equal(result, expected_output)
})

test_that("aldmckSE handles missing data correctly", {
  # 测试函数是否正确处理缺失数据
  data_with_na <- data
  data_with_na[1, 1] <- NA
  
  # 计算不带 NA 的临时数据
  tmp <- na.omit(cbind(obj$respondents[, 1:2], data_with_na))
  alpha <- tmp[, 1]
  beta <- tmp[, 2]
  z <- tmp[, 3:ncol(tmp)]
  zhat <- obj$stimuli
  sigmaj <- rep(0, length(zhat))
  for (j in 1:length(zhat)) {
    for (i in 1:length(alpha)) {
      sigmaj[j] <- sigmaj[j] + ((alpha[i] + beta[i] * zhat[j]) - z[i, j])^2
    }
  }
  for (i in 1:length(zhat)) {
    sigmaj[i] <- sqrt(sigmaj[i] / length(alpha))
  }
  
  result <- aldmckSE(obj, data_with_na)
  expect_equal(result, sigmaj)
})

test_that("aldmckSE returns correct length of standard errors", {
  # 测试函数返回的标准误差向量的长度是否正确
  result <- aldmckSE(obj, data)
  expect_length(result, length(obj$stimuli))
})