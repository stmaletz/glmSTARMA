
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

