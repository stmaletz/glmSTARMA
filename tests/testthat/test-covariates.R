# Test covariate functions TimeConstant, SpatialConstant

test_that("TimeConstant covariate works correctly", {
  expect_error(TimeConstant(1))
  expect_error(TimeConstant(c("a", "b")))
  expect_error(TimeConstant(c(1, NA)))
  expect_error(TimeConstant(matrix(1:10, nrow=2)))  
  result <- TimeConstant(c(1,2,3))
  expect_true(is.numeric(result))
  expect_true(!is.null(attr(result, "const")))
  expect_equal(attr(result, "const"), "time")
  expect_s3_class(result, "time_constant")  
})

test_that("SpatialConstant covariate works correctly", {
  expect_error(SpatialConstant(1))
  expect_error(SpatialConstant(c("a", "b")))
  expect_error(SpatialConstant(c(1, NA)))
  expect_error(SpatialConstant(matrix(1:10, nrow=2)))  
  result <- SpatialConstant(c(1,2,3))
  expect_true(is.numeric(result))
  expect_true(!is.null(attr(result, "const")))
  expect_equal(attr(result, "const"), "space")
  expect_s3_class(result, "spatial_constant")  
})



