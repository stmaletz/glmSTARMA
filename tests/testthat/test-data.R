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

