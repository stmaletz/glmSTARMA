testthat::skip_on_cran()
test_that("Data can be loaded and unloaded", {
    # clean up if any data is already cached
    delete_glmSTARMA_data(c("sst", "rota", "chickenpox"))

    dat <- load_data("sst", directory = tempdir())
    expect_true(is.list(dat))
    expect_named(dat, c("SST", "W_directed", "locations"), ignore.order = TRUE)
    
    dat2 <- load_data("rota", directory = tempdir())
    expect_true(is.list(dat2))
    expect_named(dat2, c("rota", "gdr_feature", "population_germany", "W_germany"), ignore.order = TRUE)
    
    dat3 <- load_data("chickenpox", directory = tempdir())
    expect_true(is.list(dat3))
    expect_named(dat3, c("chickenpox", "population_hungary", "W_hungary"), ignore.order = TRUE)
    # invalid name
    expect_error(load_data("invalid_dataset"), "name must be in 'rota', 'chickenpox', or 'sst'")
    
    # false if data already deleted
    expect_message(x <- delete_glmSTARMA_data("sst"), "There is no dataset to delete.")
    expect_false(x)

    # new path can be created
    new_path <- file.path(tempdir(), "new_data_path")
    dat3 <- load_data("chickenpox", directory = tempdir())

})

testthat::skip_on_cran()
test_that("data can be deleted", {
    expect_message(x <- delete_glmSTARMA_data("sst"), "There is no dataset to delete.")
    expect_false(x)

    load_data("chickenpox")
    expect_message(x <- delete_glmSTARMA_data("rota"), "is not in the cache.")
    expect_false(x)
    expect_message(x <- delete_glmSTARMA_data("chickenpox"), "Deleted: chickenpox")
    expect_true(x) 
})


testthat::skip_on_cran()
test_that("existing file is loaded without download", {
  tmp_dir <- tempdir()

  x <- data.frame(a = 1:3)
  save(x, file = file.path(tmp_dir, "rota.rda"))

  res <- load_data("rota", directory = tmp_dir)

  expect_equal(res$x, x)
})


testthat::skip_on_cran()
test_that("refresh downloads new file", {

  tmp_dir <- tempdir()

  with_mocked_bindings(
    {
      res <- load_data(
        "rota",
        directory = tmp_dir,
        refresh = TRUE
      )

      expect_equal(res$x, data.frame(a = 1:3))
    },
    "download.file" = function(url, destfile, ...) {
      x <- data.frame(a = 1:3)
      save(x, file = destfile)
      0
    },
    .package = "utils"
  )
})


testthat::skip_on_cran()
test_that("load_data throws error when download fails and no cache exists", {
  tmp_dir <- withr::local_tempdir()

  local_mocked_bindings(
    download.file = function(url, destfile, ...) stop("Simulated network error"),
    .package = "utils"
  )

  expect_error(
    load_data("rota", directory = tmp_dir),
    "Data could not be downloaded"
  )
})


testthat::skip_on_cran()
test_that("delete_glmSTARMA_data messages when directory has no .rda files", {
  tmp_dir <- withr::local_tempdir()

  local_mocked_bindings(
    R_user_dir = function(...) tmp_dir,
    .package = "tools"
  )

  expect_message(
    result <- delete_glmSTARMA_data("rota"),
    "There is no dataset to delete"
  )
  expect_false(result)
})

