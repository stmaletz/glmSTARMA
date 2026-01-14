# test fitting of dglmstarma models
testthat::skip_on_cran()

test_that("ts is validated", {
    dat <- load_data("chickenpox")
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary
    covariates <- list(population = population_hungary, 
                       season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                       season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    # valid case
    result <- dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    # not allowed: vpoisson family with dispersion model
    expect_error(dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vpoisson("log")))
    # not allowed: vbinomial family with dispersion model
    expect_error(dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vbinomial("softclipping", size = max(chickenpox))))
    # invalid case: non-matrix input
    chickenpox_vec <- as.vector(chickenpox)
    expect_error(dglmstarma(chickenpox_vec, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log")))
    # invalid case: non-numeric input
    chickenpox_char <- matrix(as.character(chickenpox), nrow = nrow(chickenpox), ncol = ncol(chickenpox))
    expect_error(dglmstarma(chickenpox_char, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vpoisson("log")))
    # invalid case: family is missing
    expect_error(dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                          mean_covariates = covariates, dispersion_covariates = covariates,
                          wlist = W_hungary))

    # invalid case: missing values
    chickenpox_na <- chickenpox
    chickenpox_na[1,10] <- NA
    expect_error(dglmstarma(chickenpox_na, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                          mean_covariates = covariates, dispersion_covariates = covariates,
                          wlist = W_hungary, mean_family = vpoisson("log")))
    # invalid case: Inf values
    chickenpox_inf <- chickenpox
    chickenpox_inf[1,10] <- Inf
    expect_error(dglmstarma(chickenpox_inf, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                          mean_covariates = covariates, dispersion_covariates = covariates,
                          wlist = W_hungary, mean_family = vpoisson("log")))
    # invalid case: negative values for count data distribution
    chickenpox_neg <- chickenpox
    chickenpox_neg[1,10] <- -5
    expect_error(dglmstarma(chickenpox_neg, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                          mean_covariates = covariates, dispersion_covariates = covariates,
                          wlist = W_hungary, mean_family = vpoisson("log")))
    expect_error(dglmstarma(chickenpox_neg, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                          mean_covariates = covariates, dispersion_covariates = covariates,
                          wlist = W_hungary, mean_family = vquasipoisson("log")))
    expect_error(dglmstarma(chickenpox_neg, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                          mean_covariates = covariates, dispersion_covariates = covariates,
                          wlist = W_hungary, mean_family = vnegative.binomial("log")))
    expect_error(dglmstarma(chickenpox_neg, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                          mean_covariates = covariates, dispersion_covariates = covariates,
                          wlist = W_hungary, mean_family = vquasibinomial("softclipping", size = max(chickenpox))))
    # invalid case: non-integer values for count data distribution
    chickenpox_nonint <- chickenpox
    chickenpox_nonint[1,10] <- 3.7
    expect_error(dglmstarma(chickenpox_nonint, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                          mean_covariates = covariates, dispersion_covariates = covariates,
                          wlist = W_hungary, mean_family = vpoisson("log")))
    expect_error(dglmstarma(chickenpox_nonint, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                          mean_covariates = covariates, dispersion_covariates = covariates,
                          wlist = W_hungary, mean_family = vquasipoisson("log")))
    expect_error(dglmstarma(chickenpox_nonint, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                          mean_covariates = covariates, dispersion_covariates = covariates,
                          wlist = W_hungary, mean_family = vnegative.binomial("log")))
    expect_error(dglmstarma(chickenpox_nonint, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                          mean_covariates = covariates, dispersion_covariates = covariates,
                          wlist = W_hungary, mean_family = vquasibinomial("softclipping", size = max(chickenpox))))
    # invalid case: values above n for binomial distribution
    chickenpox_binom <- chickenpox
    chickenpox_binom[1,10] <- max(chickenpox) + 10
    expect_error(dglmstarma(chickenpox_binom, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                          mean_covariates = covariates, dispersion_covariates = covariates,
                          wlist = W_hungary, mean_family = vquasibinomial("softclipping", size = max(chickenpox))))
    # invalid case: zero values for positive continuous distribution
    chickenpox_cont <- chickenpox
    chickenpox_cont[1,10] <- 0
    expect_error(dglmstarma(chickenpox_cont, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                          mean_covariates = covariates, dispersion_covariates = covariates,
                          wlist = W_hungary, mean_family = vgamma("log")))
    expect_error(dglmstarma(chickenpox_cont, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                          mean_covariates = covariates, dispersion_covariates = covariates,
                          wlist = W_hungary, mean_family = vinverse.gaussian("log")))
    # normal distr. case:
    result <- dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vnormal("identity"))
    expect_s3_class(result, "dglmstarma")
    delete_glmSTARMA_data("chickenpox")  # Clean up cached data
})


test_that("model is handled correctly", {
    dat <- load_data("chickenpox")
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary
    covariates <- list(population = population_hungary, 
                       season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                       season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    # valid model: no covariate orders specified (spatial orders default to 0)
    model_orders <- list(intercept = "homogeneous", past_obs = 1)
    result <- dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    expect_true(is.matrix(result$mean$model$covariates))
    expect_equal(nrow(result$mean$model$covariates), 1)
    expect_equal(ncol(result$mean$model$covariates), 3)
    expect_true(all(result$mean$model$covariates == 1L))
    expect_true(is.matrix(result$dispersion$model$covariates))
    expect_equal(nrow(result$dispersion$model$covariates), 1)
    expect_equal(ncol(result$dispersion$model$covariates), 3)
    expect_true(all(result$dispersion$model$covariates == 1L))

    # invalid model: negative time lags or zero or NA, Inf
    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1,
                         past_obs_time_lags = -1)
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log")))

    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1,
                         past_mean_time_lags = -1)
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log")))

    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1,
                         past_obs_time_lags = NA)
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log")))
    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1,
                         past_mean_time_lags = NA)
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log")))

    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1,
                         past_obs_time_lags = Inf)
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log")))

    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1,
                         past_mean_time_lags = Inf)
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log")))

    # invalid model: time lags too large >= ncol(ts)
    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1,
                         past_obs_time_lags = 600)
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log")))

    # invalid model: no regression on past observations when past_mean is specified
    model_orders <- list(intercept = "homogeneous", past_mean = 1)
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log")))
    # invalid model: intercept is not "homogeneous" or "inhomogeneous"
    model_orders <- list(intercept = "invalid_value", past_obs = 1, past_mean = 1)
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log")))
    # valid model: no regression either on past observations or past mean
    model_orders <- list(intercept = "homogeneous")
    result <- dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    # valid model: only regression on past observations
    model_orders <- list(intercept = "homogeneous", past_obs = 1)
    result <- dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    # valid model: pure intercept model
    model_orders <- list(intercept = "homogeneous")
    result <- dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        wlist = W_hungary, mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    model_orders <- list(intercept = "inhomogeneous")
    result <- dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        wlist = W_hungary, mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    
    # valid model: spatial order can be given as binary matrix
    model_orders <- list(intercept = "homogeneous", past_obs = matrix(1, nrow = 2, ncol = 1), 
                    past_mean = matrix(c(0, 1), nrow = 2, ncol = 1))
    result <- dglmstarma(chickenpox, model_orders, dispersion_model = model_orders,
                        mean_covariates = covariates, dispersion_covariates = covariates,
                        wlist = W_hungary, mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    delete_glmSTARMA_data("chickenpox")  # Clean up cached data
})


test_that("covariates dimensions are validated (mean model)", {
    dat <- load_data("chickenpox")
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary
    # valid covariates
    covariates <- list(population = population_hungary, 
                       season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                       season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    model_orders <- list(intercept = "homogeneous", past_obs = 1)
    result <- dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        mean_covariates = covariates, mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    # invalid covariates dimension (TimeConstant length)
    covariates_invalid1 <- list(population = TimeConstant(population_hungary[1:10, 1]), 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        mean_covariates = covariates_invalid1, mean_family = vquasipoisson("log")))
    # incorrect covariates dimensions (SpatialConstant length)
    covariates_invalid2 <- list(population = population_hungary, 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:52)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:52)))
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        mean_covariates = covariates_invalid2, mean_family = vquasipoisson("log")))
    # incorrect covariate dimensions (SpatialConstant is too long, warning case)
    covariates_invalid3 <- list(population = population_hungary, 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:600)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:600)))
    expect_warning(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        mean_covariates = covariates_invalid3, mean_family = vquasipoisson("log")))
    # incorrect covariate dimensions (matrix covariates with wrong nrow)
    covariates_invalid4 <- list(population = population_hungary[1:10,], 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        mean_covariates = covariates_invalid4, mean_family = vquasipoisson("log")))
    # incorrect covariate dimensions (matrix covariates with wrong ncol, too few columns)
    covariates_invalid5 <- list(population = population_hungary[, 1:100], 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        mean_covariates = covariates_invalid5, mean_family = vquasipoisson("log")))
    # incorrect covariate dimensions (matrix covariates with wrong ncol, too many columns, warning case)
    covariates_invalid6 <- list(population = cbind(population_hungary, 0), 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    expect_warning(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        mean_covariates = covariates_invalid6, mean_family = vquasipoisson("log")))
    # covariates can be named or unnamed
    covariates_unnamed <- list(population_hungary, 
                               SpatialConstant(cos(2 * pi / 52 * 1:522)),
                               SpatialConstant(sin(2 * pi / 52 * 1:522)))
    result <- dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        mean_covariates = covariates_unnamed, mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    covariates2 <- list(population = population_hungary, 
                        season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                        SpatialConstant(sin(2 * pi / 52 * 1:522)))
    result <- dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        mean_covariates = covariates2, mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    delete_glmSTARMA_data("chickenpox")  # Clean up cached data
})


test_that("covariates dimensions are validated (dispersion model)", {
    dat <- load_data("chickenpox")
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary
    # valid covariates
    covariates <- list(population = population_hungary, 
                       season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                       season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    model_orders <- list(intercept = "homogeneous", past_obs = 1)
    result <- dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        dispersion_covariates = covariates, mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    # invalid covariates dimension (TimeConstant length)
    covariates_invalid1 <- list(population = TimeConstant(population_hungary[1:10, 1]), 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        dispersion_covariates = covariates_invalid1, mean_family = vquasipoisson("log")))
    # incorrect covariates dimensions (SpatialConstant length)
    covariates_invalid2 <- list(population = population_hungary, 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:52)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:52)))
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        dispersion_covariates = covariates_invalid2, mean_family = vquasipoisson("log")))
    # incorrect covariate dimensions (SpatialConstant is too long, warning case)
    covariates_invalid3 <- list(population = population_hungary, 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:600)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:600)))
    expect_warning(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        dispersion_covariates = covariates_invalid3, mean_family = vquasipoisson("log")))
    # incorrect covariate dimensions (matrix covariates with wrong nrow)
    covariates_invalid4 <- list(population = population_hungary[1:10,], 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        dispersion_covariates = covariates_invalid4, mean_family = vquasipoisson("log")))
    # incorrect covariate dimensions (matrix covariates with wrong ncol, too few columns)
    covariates_invalid5 <- list(population = population_hungary[, 1:100], 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    expect_error(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        dispersion_covariates = covariates_invalid5, mean_family = vquasipoisson("log")))
    # incorrect covariate dimensions (matrix covariates with wrong ncol, too many columns, warning case)
    covariates_invalid6 <- list(population = cbind(population_hungary, 0), 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    expect_warning(dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        dispersion_covariates = covariates_invalid6, mean_family = vquasipoisson("log")))
    # covariates can be named or unnamed
    covariates_unnamed <- list(population_hungary, 
                               SpatialConstant(cos(2 * pi / 52 * 1:522)),
                               SpatialConstant(sin(2 * pi / 52 * 1:522)))
    result <- dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, wlist = W_hungary, 
                        dispersion_covariates = covariates_unnamed, mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    covariates2 <- list(population = population_hungary, 
                        season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                        SpatialConstant(sin(2 * pi / 52 * 1:522)))
    result <- dglmstarma(chickenpox, model_orders, dispersion_model = model_orders, 
                         wlist = W_hungary, dispersion_covariates = covariates2, 
                         mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    delete_glmSTARMA_data("chickenpox")  # Clean up cached data
})


test_that("wlist arguments are validated", {
    dat <- load_data("chickenpox")
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary
    # valid covariates
    covariates <- list(population = population_hungary, 
                       season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                       season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    # too few matrices wlist
    expect_error(dglmstarma(chickenpox, list(past_obs = 2), dispersion_model = list(past_obs = 1), 
                            wlist = W_hungary, mean_covariates = covariates,
                        dispersion_covariates = covariates, mean_family = vquasipoisson("log")))
    expect_error(dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 2), 
                            wlist = W_hungary, mean_covariates = covariates,
                        dispersion_covariates = covariates, mean_family = vquasipoisson("log")))
    # wlist is not a list
    expect_error(dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1), 
                            wlist = "not_a_list", mean_covariates = covariates,
                            dispersion_covariates = covariates, mean_family = vquasipoisson("log")))
    expect_error(dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1), 
                            wlist = W_hungary[[1]], mean_covariates = covariates,
                            dispersion_covariates = covariates, mean_family = vquasipoisson("log")))
    # wlist not specified
    expect_error(dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                            dispersion_covariates = covariates, mean_covariates = covariates,
                            wlist_past_mean = W_hungary, mean_family = vquasipoisson("log")))
    expect_error(dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                            dispersion_covariates = covariates, mean_covariates = covariates, 
                            wlist_covariates = W_hungary, mean_family = vquasipoisson("log")))
    # wlist_past_mean has not enough matrices
    expect_error(dglmstarma(chickenpox, list(past_obs = 1, past_mean = 2), dispersion_model = list(past_obs = 1),
                            wlist = W_hungary, mean_covariates = covariates,
                            dispersion_covariates = covariates, wlist_past_mean = W_hungary, mean_family = vquasipoisson("log")))
    # wlist_covariates has not enough matrices
    expect_error(dglmstarma(chickenpox, list(past_obs = 1, covariates = c(2, 2, 2)), dispersion_model = list(past_obs = 1), 
                            wlist = W_hungary, mean_covariates = covariates,
                            dispersion_covariates = covariates, wlist_covariates = W_hungary, mean_family = vquasipoisson("log")))
    # wlist_pseudo_obs has not enough matrices
    expect_error(dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 2),
                            wlist = W_hungary, mean_covariates = covariates,
                            dispersion_covariates = covariates, wlist_pseudo_obs = W_hungary, mean_family = vquasipoisson("log")))
    # wlist_past_dispersion has not enough matrices
    expect_error(dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1, past_mean = 2),
                            wlist = W_hungary, mean_covariates = covariates,
                            dispersion_covariates = covariates, wlist_past_dispersion = W_hungary, mean_family = vquasipoisson("log")))
    # wlist_covariates_dispersion has not enough matrices
    expect_error(dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1, covariates = c(2, 2, 2)),
                            wlist = W_hungary, mean_covariates = covariates,
                            dispersion_covariates = covariates, wlist_covariates_dispersion = W_hungary, mean_family = vquasipoisson("log")))

    # elements of wlist are not numeric matrices
    W_invalid <- W_hungary
    W_invalid[[1]] <- "not_a_matrix"
    expect_error(dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                            wlist = W_invalid, mean_covariates = covariates,
                            dispersion_covariates = covariates, mean_family = vquasipoisson("log")))

    # elements of wlist are not of the same dimension
    W_invalid2 <- W_hungary
    W_invalid2[[1]] <- diag(50)
    expect_error(dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                            wlist = W_invalid2, mean_covariates = covariates,
                            dispersion_covariates = covariates, mean_family = vquasipoisson("log")))
    # correct case: sparse matrices
    W_sparse <- lapply(W_hungary, function(mat) as(mat, "dgCMatrix"))
    result <- dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                            wlist = W_sparse, mean_covariates = covariates,
                            dispersion_covariates = covariates, mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    # correct case: mixed dense and sparse matrices
    W_mixed <- W_hungary
    W_mixed[[1]] <- as(W_mixed[[1]], "dgCMatrix")
    result <- dglmstarma(chickenpox, list(past_obs = 1), dispersion_model = list(past_obs = 1),
                            wlist = W_mixed, mean_covariates = covariates,
                            dispersion_covariates = covariates, mean_family = vquasipoisson("log"))
    expect_s3_class(result, "dglmstarma")
    delete_glmSTARMA_data("chickenpox")  # Clean up cached data
})


