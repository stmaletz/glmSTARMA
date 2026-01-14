
test_that("Data can be loaded and unloaded", {
    dat <- load_data("SST")
    expect_true(is.list(dat))
    expect_named(dat, c("SST", "W_directed", "locations"))
    
    dat2 <- load_data("rota")
    expect_true(is.list(dat2))
    expect_named(dat2, c("rota", "gdr_feature", "population_germany", "W_germany"))
    
    dat3 <- load_data("chickenpox")
    expect_true(is.list(dat3))
    expect_named(dat3, c("chickenpox", "population_hungary", "W_hungary"))

    # invalid name
    expect_error(load_data("invalid_dataset"), "name must be in 'rota', 'chickenpox', or 'SST'")

    # clean up
    expect_true(delete_glmSTARMA_data(c("SST", "rota", "chickenpox")))
    # false if data already deleted
    expect_message(x <- delete_glmSTARMA_data("SST"), "There is no dataset to delete.")
    expect_false(x)
})
