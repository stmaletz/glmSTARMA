
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
    predictions <- predict(result, type = "sample")
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

    # with copula
    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        family = vpoisson("log", copula = "frank", copula_param = 2))
    expect_is(predictions, "matrix")
    expect_true(all(predictions >= 0))
    expect_true(all(predictions == round(predictions)))

    # With covariates, new observations, new covariate values
    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        covariates = covariates,
                        family = vpoisson("log"))
    expect_warning(predictions <- predict(result, n.ahead = 10, type = "response"))

    predictions <- predict(result, n.ahead = 10, type = "response", newxreg = list(population = population_hungary[, 1:10], 
                                                season_cos = SpatialConstant(cos(2 * pi / 52 * 523:532)),
                                                season_sin = SpatialConstant(sin(2 * pi / 52 * 523:532))))


    # too few new covariate values
    expect_warning(predictions <- predict(result, n.ahead = 10, type = "response", newxreg = list(population = population_hungary[, 1:5], 
                                                season_cos = SpatialConstant(cos(2 * pi / 52 * 523:532)),
                                                season_sin = SpatialConstant(sin(2 * pi / 52 * 523:525)))))                                            

    set.seed(123)
    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        covariates = covariates,
                        family = vquasipoisson("log", dispersion = matrix(runif(length(chickenpox), min = 1, max = 10), nrow = nrow(chickenpox), ncol = ncol(chickenpox))))
    expect_error(predict(result, n.ahead = 10, type = "sample", newxreg = list(population = population_hungary[1:10, ], 
                                                season_cos = SpatialConstant(cos(2 * pi / 52 * 523:532)),
                                                season_sin = SpatialConstant(sin(2 * pi / 52 * 523:532)))))


    # new observations provided                                            
    set.seed(123)
    predictions <- predict(result, n.ahead = 10, type = "response", newobs = matrix(rpois(200, 10), nrow = 20, ncol = 10), 
                       newxreg = list(population = population_hungary[, 1:10], 
                                      season_cos = SpatialConstant(cos(2 * pi / 52 * 523:532)),
                                      season_sin = SpatialConstant(sin(2 * pi / 52 * 523:532))))

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
                        mean_family = vnegative.binomial("log"),
                        control = list(maxit = 100, print_progress = FALSE))
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
                        mean_family = vnegative.binomial("log"),
                        control = list(maxit = 100, print_progress = FALSE))
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
                        mean_family = vquasipoisson("log", copula = "frank", copula_param = 2),
                        control = list(maxit = 100, print_progress = FALSE))
    predictions <- predict(result)
    expect_is(predictions, "list")
    expect_equal(length(predictions), 2)
    expect_named(predictions, c("mean", "dispersion"), ignore.order = TRUE)
    expect_is(predictions$mean, "matrix")
    expect_is(predictions$dispersion, "matrix")


    result <- dglmstarma(chickenpox, list(past_obs = 1), list(past_obs = 1), wlist = W_hungary, 
                        mean_family = vquasipoisson("log", copula = "normal", copula_param = 0.5),
                        control = list(maxit = 100, print_progress = FALSE))
    predictions <- predict(result)
    expect_is(predictions, "list")
    expect_equal(length(predictions), 2)
    expect_named(predictions, c("mean", "dispersion"), ignore.order = TRUE)
    expect_is(predictions$mean, "matrix")
    expect_is(predictions$dispersion, "matrix")

    result <- dglmstarma(chickenpox, list(past_obs = 1), list(past_obs = 1), wlist = W_hungary, 
                        mean_family = vquasipoisson("log"),
                        control = list(maxit = 100, print_progress = FALSE))
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

    set.seed(42)
    predictions <- predict(result, n.ahead = 10, type = "sample")
    expect_is(predictions, "list")
    expect_true(all(predictions$mean >= 0))
    expect_true(all(predictions$mean == round(predictions$mean)))

    result <- dglmstarma(chickenpox, list(past_obs = 1), list(past_obs = 1), wlist = W_hungary, 
                        mean_family = vquasipoisson("log", copula = "frank", copula_param = 2),
                        control = list(maxit = 100, print_progress = FALSE))
    predictions <- predict(result, n.ahead = 10, type = "sample")
    expect_is(predictions, "list")
    expect_true(all(predictions$mean >= 0))
    expect_true(all(predictions$mean == round(predictions$mean)))                    

    # With covariates, new observations, new covariate values

    result <- dglmstarma(chickenpox, list(past_obs = 1), list(past_obs = 1), wlist = W_hungary, 
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        mean_family = vquasipoisson("log"),
                        control = list(maxit = 100, print_progress = FALSE))

    # warning because new covariates are not provided
    expect_warning(predictions <- predict(result, n.ahead = 10, type = "response"))

    # with new covariate values 
    predictions <- predict(result, n.ahead = 10, type = "response", 
                            newxreg_mean = list(population = population_hungary[, 1:10], 
                                                season_cos = SpatialConstant(cos(2 * pi / 52 * 523:532)),
                                                season_sin = SpatialConstant(sin(2 * pi / 52 * 523:532))),
                            newxreg_dispersion = list(population = population_hungary[, 1:10],
                                                season_cos = SpatialConstant(cos(2 * pi / 52 * 523:532)),
                                                season_sin = SpatialConstant(sin(2 * pi / 52 * 523:532))))
    expect_is(predictions, "list")
    expect_equal(length(predictions), 2)
    expect_named(predictions, c("mean", "dispersion"), ignore.order = TRUE)

    # with new observations
    set.seed(123)
    predictions <- predict(result, n.ahead = 10, type = "response", newobs = matrix(rpois(200, 10), nrow = 20, ncol = 10),
                            newxreg_mean = list(population = population_hungary[, 1:10], 
                                                season_cos = SpatialConstant(cos(2 * pi / 52 * 523:532)),
                                                season_sin = SpatialConstant(sin(2 * pi / 52 * 523:532))),
                            newxreg_dispersion = list(population = population_hungary[, 1:10],
                                                season_cos = SpatialConstant(cos(2 * pi / 52 * 523:532)),
                                                season_sin = SpatialConstant(sin(2 * pi / 52 * 523:532))))
    expect_is(predictions, "list")
    expect_equal(length(predictions), 2)
    expect_named(predictions, c("mean", "dispersion"), ignore.order = TRUE)

})














