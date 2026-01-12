# glmstarma:

testthat::skip_on_cran()

test_that("ts is validated", {
    data("chickenpox")
    covariates <- list(population = population_hungary, 
                       season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                       season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    # valid case
    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log"))
    expect_s3_class(result, "glmstarma")
    # invalid case: non-matrix input
    chickenpox_vec <- as.vector(chickenpox)
    expect_error(glmstarma(chickenpox_vec, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vpoisson("log")))
    # invalid case: non-numeric input
    chickenpox_char <- matrix(as.character(chickenpox), nrow = nrow(chickenpox), ncol = ncol(chickenpox))
    expect_error(glmstarma(chickenpox_char, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vpoisson("log")))
    # invalid case: family is missing
    expect_error(glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates))

    # invalid case: missing values
    chickenpox_na <- chickenpox
    chickenpox_na[1,10] <- NA
    expect_error(glmstarma(chickenpox_na, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vpoisson("log")))
    # invalid case: Inf values
    chickenpox_inf <- chickenpox
    chickenpox_inf[1,10] <- Inf
    expect_error(glmstarma(chickenpox_inf, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vpoisson("log")))
    # invalid case: negative values for count data distribution
    chickenpox_neg <- chickenpox
    chickenpox_neg[1,10] <- -5
    expect_error(glmstarma(chickenpox_neg, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vpoisson("log")))
    expect_error(glmstarma(chickenpox_neg, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vquasipoisson("log")))
    expect_error(glmstarma(chickenpox_neg, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vnegative.binomial("log")))
    expect_error(glmstarma(chickenpox_neg, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vbinomial("softclipping", size = max(chickenpox))))
    expect_error(glmstarma(chickenpox_neg, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vquasibinomial("softclipping", size = max(chickenpox))))
    # invalid case: non-integer values for count data distribution
    chickenpox_nonint <- chickenpox
    chickenpox_nonint[1,10] <- 3.7
    expect_error(glmstarma(chickenpox_nonint, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vpoisson("log")))
    expect_error(glmstarma(chickenpox_nonint, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vquasipoisson("log")))
    expect_error(glmstarma(chickenpox_nonint, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vnegative.binomial("log")))
    expect_error(glmstarma(chickenpox_nonint, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vbinomial("softclipping", size = max(chickenpox))))
    expect_error(glmstarma(chickenpox_nonint, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vquasibinomial("softclipping", size = max(chickenpox))))
    # invalid case: values above n for binomial distribution
    chickenpox_binom <- chickenpox
    chickenpox_binom[1,10] <- max(chickenpox) + 10
    expect_error(glmstarma(chickenpox_binom, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vbinomial("softclipping", size = max(chickenpox))))
    expect_error(glmstarma(chickenpox_binom, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vquasibinomial("softclipping", size = max(chickenpox))))
    # invalid case: zero values for positive continuous distribution
    chickenpox_cont <- chickenpox
    chickenpox_cont[1,10] <- 0
    expect_error(glmstarma(chickenpox_cont, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vgamma("log")))
    expect_error(glmstarma(chickenpox_cont, list(past_obs = 1), wlist = W_hungary, 
                          covariates = covariates, family = vinverse.gaussian("log")))
    # normal distr. case:
    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        covariates = covariates, family = vnormal("identity"))
    expect_s3_class(result, "glmstarma")
})


test_that("model is handled correctly", {
    data("chickenpox")
    covariates <- list(population = population_hungary, 
                       season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                       season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    # valid model: no covariate orders specified (spatial orders default to 0)
    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1)
    result <- glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log"))
    expect_s3_class(result, "glmstarma")
    expect_true(is.matrix(result$model$covariates))
    expect_equal(nrow(result$model$covariates), 1)
    expect_equal(ncol(result$model$covariates), 3)
    expect_true(all(result$model$covariates == 1L))

    # invalid model: negative time lags or zero or NA, Inf
    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1,
                         past_obs_time_lags = -1)
    expect_error(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log")))

    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1,
                         past_mean_time_lags = -1)
    expect_error(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log")))

    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1,
                         past_obs_time_lags = NA)
    expect_error(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log")))
    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1,
                         past_mean_time_lags = NA)
    expect_error(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log")))

    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1,
                         past_obs_time_lags = Inf)
    expect_error(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log")))

    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1,
                         past_mean_time_lags = Inf)
    expect_error(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log")))

    # invalid model: time lags too large >= ncol(ts)
    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1,
                         past_obs_time_lags = 600)
    expect_error(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log")))

    # invalid model: no regression on past observations when past_mean is specified
    model_orders <- list(intercept = "homogeneous", past_mean = 1)
    expect_error(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log")))
    # invalid model: intercept is not "homogeneous" or "inhomogeneous"
    model_orders <- list(intercept = "invalid_value", past_obs = 1, past_mean = 1)
    expect_error(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log")))
    # valid model: no regression either on past observations or past mean
    model_orders <- list(intercept = "homogeneous")
    result <- glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log"))
    expect_s3_class(result, "glmstarma")
    # valid model: only regression on past observations
    model_orders <- list(intercept = "homogeneous", past_obs = 1)
    result <- glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log"))
    expect_s3_class(result, "glmstarma")
    # valid model: pure intercept model
    model_orders <- list(intercept = "homogeneous")
    result <- glmstarma(chickenpox, model_orders, wlist = W_hungary, family = vpoisson("log"))
    expect_s3_class(result, "glmstarma")
    model_orders <- list(intercept = "inhomogeneous")
    result <- glmstarma(chickenpox, model_orders, wlist = W_hungary, family = vpoisson("log"))
    expect_s3_class(result, "glmstarma")
    
    # valid model: spatial order can be given as binary matrix
    model_orders <- list(intercept = "homogeneous", past_obs = matrix(1, nrow = 2, ncol = 1), 
                    past_mean = matrix(c(0, 1), nrow = 2, ncol = 1))
    result <- glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log"))
    expect_s3_class(result, "glmstarma")
})


test_that("covariates dimensions are validated", {
    data("chickenpox")
    # valid covariates
    covariates <- list(population = population_hungary, 
                       season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                       season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    model_orders <- list(intercept = "homogeneous", past_obs = 1, past_mean = 1)
    result <- glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log"))
    expect_s3_class(result, "glmstarma")
    # invalid covariates dimension (TimeConstant length)
    covariates_invalid1 <- list(population = TimeConstant(population_hungary[1:10, 1]), 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    expect_error(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates_invalid1, family = vpoisson("log")))
    # incorrect covariates dimensions (SpatialConstant length)
    covariates_invalid2 <- list(population = population_hungary, 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:52)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:52)))
    expect_error(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates_invalid2, family = vpoisson("log")))
    # incorrect covariate dimensions (SpatialConstant is too long, warning case)
    covariates_invalid3 <- list(population = population_hungary, 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:600)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:600)))
    expect_warning(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates_invalid3, family = vpoisson("log")))
    # incorrect covariate dimensions (matrix covariates with wrong nrow)
    covariates_invalid4 <- list(population = population_hungary[1:10,], 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    expect_error(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates_invalid4, family = vpoisson("log")))
    # incorrect covariate dimensions (matrix covariates with wrong ncol, too few columns)
    covariates_invalid5 <- list(population = population_hungary[, 1:100], 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    expect_error(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates_invalid5, family = vpoisson("log")))
    # incorrect covariate dimensions (matrix covariates with wrong ncol, too many columns, warning case)
    covariates_invalid6 <- list(population = cbind(population_hungary, 0), 
                                season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                                season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    expect_warning(glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates_invalid6, family = vpoisson("log")))
    # covariates can be named or unnamed
    covariates_unnamed <- list(population_hungary, 
                               SpatialConstant(cos(2 * pi / 52 * 1:522)),
                               SpatialConstant(sin(2 * pi / 52 * 1:522)))
    result <- glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates_unnamed, family = vpoisson("log"))
    expect_s3_class(result, "glmstarma")
    covariates2 <- list(population = population_hungary, 
                        season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                        SpatialConstant(sin(2 * pi / 52 * 1:522)))
    result <- glmstarma(chickenpox, model_orders, wlist = W_hungary, 
                        covariates = covariates2, family = vpoisson("log"))
    expect_s3_class(result, "glmstarma")
})



test_that("wlist arguments are validated", {
    data("chickenpox")
    # valid covariates
    covariates <- list(population = population_hungary, 
                       season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                       season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    # too few matrices wlist
    expect_error(glmstarma(chickenpox, list(past_obs = 2), wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log")))
    # wlist is not a list
    expect_error(glmstarma(chickenpox, list(past_obs = 2), wlist = "not_a_list", 
                        covariates = covariates, family = vpoisson("log")))
    expect_error(glmstarma(chickenpox, list(past_obs = 2), wlist = W_hungary[[1]], 
                        covariates = covariates, family = vpoisson("log")))
    # wlist not specified
    expect_error(glmstarma(chickenpox, list(past_obs = 1), covariates = covariates, 
                            wlist_past_mean = W_hungary, family = vpoisson("log")))
    expect_error(glmstarma(chickenpox, list(past_obs = 1), covariates = covariates, 
                            wlist_covariates = W_hungary, family = vpoisson("log")))
    # wlist_past_mean has not enough matrices
    expect_error(glmstarma(chickenpox, list(past_obs = 1, past_mean = 2), wlist = W_hungary, 
                            covariates = covariates, wlist_past_mean = W_hungary, family = vpoisson("log")))
    # wlist_covariates has not enough matrices
    expect_error(glmstarma(chickenpox, list(past_obs = 1, covariates = c(2, 2, 2)), wlist = W_hungary, 
                            covariates = covariates, wlist_covariates = W_hungary, family = vpoisson("log")))
    # elements of wlist are not numeric matrices
    W_invalid <- W_hungary
    W_invalid[[1]] <- "not_a_matrix"
    expect_error(glmstarma(chickenpox, list(past_obs = 1), wlist = W_invalid, 
                        covariates = covariates, family = vpoisson("log")))
    # elements of wlist are not of the same dimension
    W_invalid2 <- W_hungary
    W_invalid2[[1]] <- diag(50)
    expect_error(glmstarma.sim(200, parameter, model_orders, W_invalid2, covariates, family = fam))
    # correct case: sparse matrices
    W_sparse <- lapply(W_hungary, function(mat) as(mat, "dgCMatrix"))
    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_sparse, 
                        covariates = covariates, family = vpoisson("log"))
    expect_s3_class(result, "glmstarma")
    # correct case: mixed dense and sparse matrices
    W_mixed <- W_hungary
    W_mixed[[1]] <- as(W_mixed[[1]], "dgCMatrix")
    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_mixed, 
                        covariates = covariates, family = vpoisson("log"))
    expect_s3_class(result, "glmstarma")
})
