# dimension of matrices is correct
test_that("dimension of W matrices is correct", {
  methods <- c("rectangle", "line", "circle", "full", "independent")
  dim <- 30
  maxOrder <- 3
  width <- 5
  for(method in methods){
    Wlist <- generateW(method = method, dim = dim, maxOrder = maxOrder, width = width)
    for(W in Wlist){
      expect_equal(dim(W), c(dim, dim))
    }
  }
})


# length of all list is correct
test_that("length of Wlist is correct", {
  methods <- c("rectangle", "line", "circle")
  dim <- 50
  maxOrder <- 4
  width <- 5
  for(method in methods){
    Wlist <- generateW(method = method, dim = dim, maxOrder = maxOrder, width = width)
    expect_equal(length(Wlist), maxOrder + 1)
  }
  expect_equal(length(generateW(method = "full", dim = dim)), dim^2)
  expect_equal(length(generateW(method = "independent", dim = dim)), dim)
})



# test functions: rectangle, line, circle -> row sums = 1
test_that("matrices are row-normalized", {
  methods <- c("rectangle", "line", "circle")
  dim <- 100
  maxOrder <- 5
  width <- 10
  for(method in methods){
    Wlist <- generateW(method = method, dim = dim, maxOrder = maxOrder, width = width)
    for(W in Wlist){
      row_sums <- rowSums(W)
      expect_true(all(abs(row_sums - 1) < 1e-10))
    }
  }
})

# test entries of W matrices
test_that("entries of W matrices are correct", {
  # independent
  dim <- 10
  Wlist_indep <- generateW(method = "independent", dim = dim)
  for(i in 1:dim){
    W <- Wlist_indep[[i]]
    expect_equal(W[i, i], 1)
    expect_equal(sum(W > 0) , 1)
  }
  # full
  Wlist_full <- generateW(method = "full", dim = dim)
  for(i in 1:(dim^2)){
    W <- Wlist_full[[i]]
    expect_equal(sum(W), 1)
    expect_equal(sum(W != 0), 1)
  }
})

test_that("maxOrder is valid",{
  dim <- 20
  expect_warning(generateW(method = "rectangle", dim = dim, maxOrder = -2, width = 4))
  expect_warning(generateW(method = "line", dim = dim, maxOrder = -5))
  expect_warning(generateW(method = "circle", dim = dim, maxOrder = -2))
  expect_no_message(generateW(method = "rectangle", dim = dim, maxOrder = 15, width = 4))
  expect_error(generateW(method = "line", dim = dim, maxOrder = 11))
  expect_warning(generateW(method = "circle", dim = dim, maxOrder = 11))
})

test_that("width is valid for rectangle method",{
  dim <- 30
  maxOrder <- 4
  expect_error(generateW(method = "rectangle", dim = dim, maxOrder = maxOrder, width = 0))
  expect_error(generateW(method = "rectangle", dim = dim, maxOrder = maxOrder, width = -3))
  expect_error(generateW(method = "rectangle", dim = dim, maxOrder = maxOrder, width = 7))
})
test_that("method argument is valid",{
  dim <- 25
  maxOrder <- 3
  width <- 5
  expect_error(generateW(method = "invalid_method", dim = dim, maxOrder = maxOrder, width = width))
})
