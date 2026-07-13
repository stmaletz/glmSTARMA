
testthat::skip_on_cran()
test_that("predict works for glmstarma objects", {
    dat <- load_data("chickenpox", directory = tempdir())
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary
    covariates <- list(population = population_hungary, 
                       season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                       season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        family = vpoisson("log"))
    # test predict function
    predictions <- predict(result, n.ahead = 10)
    expect_is(predictions, "matrix")
    expect_equal(ncol(predictions), 10)
    expect_equal(nrow(predictions), nrow(chickenpox))

    # inhomogeneous intercept
    result <- glmstarma(chickenpox, list(intercept = "inhomogeneous", past_obs = 1), wlist = W_hungary, 
                        family = vpoisson("log"))
    # test predict function
    predictions <- predict(result, n.ahead = 10)
    expect_is(predictions, "matrix")
    expect_equal(ncol(predictions), 10)
    expect_equal(nrow(predictions), nrow(chickenpox))

    # different copulas
    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        family = vpoisson("log", copula = "frank", copula_param = 2))
    predictions <- predict(result)
    expect_is(predictions, "matrix")

    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        family = vpoisson("log", copula = "normal", copula_param = 0.5))
    predictions <- predict(result)
    expect_is(predictions, "matrix")

    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        family = vpoisson("log"))
    # test predict function
    predictions <- predict(result, n.ahead = 10, type = "response")
    expect_is(predictions, "matrix")

    predictions <- predict(result, n.ahead = 10, type = "link")
    expect_is(predictions, "matrix")

    predictions <- predict(result, n.ahead = 10, type = "sample")
    expect_is(predictions, "matrix")
    expect_true(all(predictions >= 0))
    expect_true(all(predictions == round(predictions)))

    # TODO: With covariates, new observations, new covariate values

})






testthat::skip_on_cran()
test_that("predict works for dglmstarma objects", {
    dat <- load_data("chickenpox", directory = tempdir())
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary
    covariates <- list(population = population_hungary, 
                       season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                       season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    result <- dglmstarma(chickenpox, list(past_obs = 1), list(past_obs = 1), wlist = W_hungary, 
                        mean_family = vnegative.binomial("log"))
    # test predict function
    predictions <- predict(result, n.ahead = 10)
    expect_is(predictions, "list")
    expect_equal(length(predictions), 2)
    expect_named(predictions, c("mean", "dispersion"), ignore.order = TRUE)
    expect_is(predictions$mean, "matrix")
    expect_is(predictions$dispersion, "matrix")
    expect_equal(ncol(predictions$mean), 10)
    expect_equal(nrow(predictions$mean), nrow(chickenpox))
    expect_equal(ncol(predictions$dispersion), 10)
    expect_equal(nrow(predictions$dispersion), nrow(chickenpox))

    # inhomogeneous intercept
    result <- dglmstarma(chickenpox, list(intercept = "inhomogeneous", past_obs = 1), list(past_obs = 1), wlist = W_hungary, 
                        mean_family = vnegative.binomial("log"))
    # test predict function
    predictions <- predict(result, n.ahead = 10)
    expect_is(predictions, "list")
    expect_equal(length(predictions), 2)
    expect_named(predictions, c("mean", "dispersion"), ignore.order = TRUE)
    expect_is(predictions$mean, "matrix")
    expect_is(predictions$dispersion, "matrix")
    expect_equal(ncol(predictions$mean), 10)
    expect_equal(nrow(predictions$mean), nrow(chickenpox))
    expect_equal(ncol(predictions$dispersion), 10)
    expect_equal(nrow(predictions$dispersion), nrow(chickenpox))

    # different copulas
    result <- dglmstarma(chickenpox, list(past_obs = 1), list(past_obs = 1), wlist = W_hungary, 
                        mean_family = vquasipoisson("log", copula = "frank", copula_param = 2))
    predictions <- predict(result)
    expect_is(predictions, "list")
    expect_equal(length(predictions), 2)
    expect_named(predictions, c("mean", "dispersion"), ignore.order = TRUE)
    expect_is(predictions$mean, "matrix")
    expect_is(predictions$dispersion, "matrix")


    result <- dglmstarma(chickenpox, list(past_obs = 1), list(past_obs = 1), wlist = W_hungary, 
                        mean_family = vquasipoisson("log", copula = "normal", copula_param = 0.5))
    predictions <- predict(result)
    expect_is(predictions, "list")
    expect_equal(length(predictions), 2)
    expect_named(predictions, c("mean", "dispersion"), ignore.order = TRUE)
    expect_is(predictions$mean, "matrix")
    expect_is(predictions$dispersion, "matrix")

    result <- dglmstarma(chickenpox, list(past_obs = 1), list(past_obs = 1), wlist = W_hungary, 
                        mean_family = vquasipoisson("log"))
    # test predict function
    predictions <- predict(result, n.ahead = 10, type = "response")
    expect_is(predictions, "list")
    expect_equal(length(predictions), 2)
    expect_named(predictions, c("mean", "dispersion"), ignore.order = TRUE)
    expect_is(predictions$mean, "matrix")
    expect_is(predictions$dispersion, "matrix")

    predictions <- predict(result, n.ahead = 10, type = "link")
    expect_is(predictions, "list")
    expect_equal(length(predictions), 2)
    expect_named(predictions, c("mean", "dispersion"), ignore.order = TRUE)
    expect_is(predictions$mean, "matrix")
    expect_is(predictions$dispersion, "matrix")

    predictions <- predict(result, n.ahead = 10, type = "sample")
    expect_is(predictions, "list")
    expect_true(all(predictions$mean >= 0))
    expect_true(all(predictions$mean == round(predictions$mean)))

    # TODO: With covariates, new observations, new covariate values

})














