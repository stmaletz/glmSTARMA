# dglmstarma_sim:

test_that("ntime is handled correctly", {
    set.seed(42)
    n_obs <- 200L
    W <- generateW("rectangle", 100, 2, 10)
    model_orders_mean <- list(intercept = "homogeneous", 
                              past_obs = 2, past_mean = 1, 
                              covariates = c(0, 0))
    model_orders_dispersion <- list(intercept = "homogeneous", 
                                    past_obs = 1, 
                                    covariates = c(0, 0))
    covariates_mean <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(n_obs))),
                            location = TimeConstant(rnorm(100, sd = 0.81)))

    covariates_dispersion <- list(season = SpatialConstant(sin(2 * pi / 24 * seq(n_obs))),
                                  location = TimeConstant(runif(100)))
    params_mean <- list(intercept = 0.6, 
                        past_mean = matrix(c(0.2, 0.1), nrow = 2), 
                        past_obs = matrix(c(0.2, 0.1, 0.05), nrow = 3), 
                        covariates = matrix(c(0.9, 0.2), ncol = 2))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = matrix(c(0.5, 0.2), nrow = 2), 
                              covariates = matrix(c(0.1, 0.75), ncol = 2))
    family <- vnormal(copula = "frank", copula_param = 2)
    

    # correct n_obs input
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))
    expect_type(result, "list")
    expect_named(
        result,
        c("observations", "link_values", "pseudo_observations", "dispersion_values", 
          "mean_model", "dispersion_model", "parameters_mean", "parameters_dispersion")
    )
    expect_length(result, 8)

    expect_type(result$observations, "double")
    expect_true(is.matrix(result$observations))
    expect_equal(dim(result$observations), c(100, 200))

    expect_type(result$link_values, "double")
    expect_true(is.matrix(result$link_values))
    expect_equal(dim(result$link_values), c(100, 200))

    # incorrect n_obs input
    expect_error(dglmstarma.sim(NA, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    expect_error(dglmstarma.sim(NULL, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    expect_error(dglmstarma.sim(Inf, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    expect_error(dglmstarma.sim(-10, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    expect_error(dglmstarma.sim(0, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    expect_error(dglmstarma.sim(10.5, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
})

test_that("n_start is handled correctly", {
    set.seed(42)
    n_obs <- 200L
    W <- generateW("rectangle", 100, 2, 10)
    model_orders_mean <- list(intercept = "homogeneous", 
                              past_obs = 2, past_mean = 1, 
                              covariates = c(0, 0))
    model_orders_dispersion <- list(intercept = "homogeneous", 
                                    past_obs = 1, 
                                    covariates = c(0, 0))
    covariates_mean <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(n_obs))),
                            location = TimeConstant(rnorm(100, sd = 0.81)))

    covariates_dispersion <- list(season = SpatialConstant(sin(2 * pi / 24 * seq(n_obs))),
                                  location = TimeConstant(runif(100)))
    params_mean <- list(intercept = 0.6, 
                        past_mean = matrix(c(0.2, 0.1), nrow = 2), 
                        past_obs = matrix(c(0.2, 0.1, 0.05), nrow = 3), 
                        covariates = matrix(c(0.9, 0.2), ncol = 2))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = matrix(c(0.5, 0.2), nrow = 2), 
                              covariates = matrix(c(0.1, 0.75), ncol = 2))
    family <- vnormal(copula = "frank", copula_param = 2)
    # correct n_obs input
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion,
                   control = list(return_burn_in = TRUE))
    expect_type(result$observations, "double")
    expect_true(is.matrix(result$observations))
    expect_equal(dim(result$observations), c(100, 300))

    expect_type(result$link_values, "double")
    expect_true(is.matrix(result$link_values))
    expect_equal(dim(result$link_values), c(100, 300))
    # incorrect nstart input
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion,
                   n_start = -10L))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion,
                               n_start = 0L))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion,
                               n_start = 10.5))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion,
                               n_start = NA))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion,
                               n_start = Inf))
    # n_start is too small
    model_orders_mean$past_obs_time_lags <- 10L
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion,
                               n_start = 5L))
})

test_that("control gets checked", {
    set.seed(42)
    n_obs <- 200L
    W <- generateW("rectangle", 100, 2, 10)
    model_orders_mean <- list(intercept = "homogeneous", 
                              past_obs = 2, past_mean = 1, 
                              covariates = c(0, 0))
    model_orders_dispersion <- list(intercept = "homogeneous", 
                                    past_obs = 1, 
                                    covariates = c(0, 0))
    covariates_mean <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(n_obs))),
                            location = TimeConstant(rnorm(100, sd = 0.81)))

    covariates_dispersion <- list(season = SpatialConstant(sin(2 * pi / 24 * seq(n_obs))),
                                  location = TimeConstant(runif(100)))
    params_mean <- list(intercept = 0.6, 
                        past_mean = matrix(c(0.2, 0.1), nrow = 2), 
                        past_obs = matrix(c(0.2, 0.1, 0.05), nrow = 3), 
                        covariates = matrix(c(0.9, 0.2), ncol = 2))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = matrix(c(0.5, 0.2), nrow = 2), 
                              covariates = matrix(c(0.1, 0.75), ncol = 2))
    family <- vnormal(copula = "frank", copula_param = 2)
    
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion,
                    control = "not_a_list"))
})



test_that("wlist arguments are validated", {
    set.seed(42)
    n_obs <- 200L
    W <- generateW("rectangle", 100, 2, 10)
    model_orders_mean <- list(intercept = "homogeneous", 
                              past_obs = 2, past_mean = 1, 
                              covariates = c(0, 0))
    model_orders_dispersion <- list(intercept = "homogeneous", 
                                    past_obs = 1, 
                                    covariates = c(0, 0))
    covariates_mean <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(n_obs))),
                            location = TimeConstant(rnorm(100, sd = 0.81)))

    covariates_dispersion <- list(season = SpatialConstant(sin(2 * pi / 24 * seq(n_obs))),
                                  location = TimeConstant(runif(100)))
    params_mean <- list(intercept = 0.6, 
                        past_mean = matrix(c(0.2, 0.1), nrow = 2), 
                        past_obs = matrix(c(0.2, 0.1, 0.05), nrow = 3), 
                        covariates = matrix(c(0.9, 0.2), ncol = 2))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = matrix(c(0.5, 0.2), nrow = 2), 
                              covariates = matrix(c(0.1, 0.75), ncol = 2))
    family <- vnormal(copula = "frank", copula_param = 2)
    # too few matrices wlist
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W[1:2], pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    # wlist is not a list
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = "not_a_list", pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W[[1]], pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    # wlist not specified
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist_past_mean = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist_covariates = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    # wlist_past_mean has not enough matrices
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion,
                   wlist_past_mean = W[1]))
    # wlist_covariates has not enough matrices
    model_orders2 <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1, covariates = c(1, 1))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders2, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion,
                   wlist_covariates = W[1]))
    # elements of wlist are not numeric matrices
    W_invalid <- W
    W_invalid[[1]] <- "not_a_matrix"
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W_invalid, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    # elements of wlist are not of the same dimension
    W_invalid2 <- W
    W_invalid2[[1]] <- diag(50)
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W_invalid2, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    # correct case: sparse matrices
    W_sparse <- lapply(W, function(mat) as(mat, "dgCMatrix"))
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W_sparse, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))
    # correct case: mixed dense and sparse matrices
    W_mixed <- W
    W_mixed[[1]] <- as(W_mixed[[1]], "dgCMatrix")
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W_mixed, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))
})


test_that("model is handled correctly", {
    set.seed(42)
    n_obs <- 200L
    W <- generateW("rectangle", 100, 2, 10)
    model_orders_mean <- list(intercept = "homogeneous", 
                              past_obs = 2, past_mean = 1, 
                              covariates = c(0, 0))
    model_orders_dispersion <- list(intercept = "homogeneous", 
                                    past_obs = 1, 
                                    covariates = c(0, 0))
    covariates_mean <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(n_obs))),
                            location = TimeConstant(rnorm(100, sd = 0.81)))

    covariates_dispersion <- list(season = SpatialConstant(sin(2 * pi / 24 * seq(n_obs))),
                                  location = TimeConstant(runif(100)))
    params_mean <- list(intercept = 0.6, 
                        past_mean = matrix(c(0.2, 0.1), nrow = 2), 
                        past_obs = matrix(c(0.2, 0.1, 0.05), nrow = 3), 
                        covariates = matrix(c(0.9, 0.2), ncol = 2))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = matrix(c(0.5, 0.2), nrow = 2), 
                              covariates = matrix(c(0.1, 0.75), ncol = 2))
    family <- vnormal(copula = "frank", copula_param = 2)
    # invalid model: negative time lags or zero or NA, Inf
    model_orders_mean <- list(intercept = "homogeneous", 
                              past_obs = 2, past_mean = 1, 
                              past_obs_time_lags = -1,
                              covariates = c(0, 0))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))

    model_orders_mean <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1,
                         past_mean_time_lags = -1)
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))

    model_orders_mean <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1,
                         past_obs_time_lags = NA)
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))

    model_orders_mean <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1,
                         past_mean_time_lags = NA)
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))

    model_orders_mean <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1,
                         past_obs_time_lags = Inf)
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))

    model_orders_mean <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1,
                         past_mean_time_lags = Inf)
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))

    # invalid model: no regression on past observations when past_mean is specified
    model_orders_mean <- list(intercept = "homogeneous", past_mean = 1)
    parameter2 <- params_mean
    parameter2$past_obs <- NULL
    expect_error(dglmstarma.sim(n_obs, parameter2, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))

    # invalid model: intercept is not "homogeneous" or "inhomogeneous"
    model_orders_mean <- list(intercept = "invalid_value", past_obs = 2, past_mean = 1)
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))

    # valid model: no regression either on past observations or past mean
    model_orders_mean <- list(intercept = "homogeneous")
    parameter2 <- params_mean
    parameter2$past_obs <- NULL
    parameter2$past_mean <- NULL
    result <- dglmstarma.sim(n_obs, parameter2, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))

    # valid model: only regression on past observations
    model_orders_mean <- list(intercept = "homogeneous", past_obs = 2)
    parameter2 <- params_mean
    parameter2$past_mean <- NULL
    result <- dglmstarma.sim(n_obs, parameter2, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))

    # valid model: no covariate orders specified (spatial orders default to 0)
    model_orders_mean <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1)
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))
    expect_true(is.matrix(result$mean_model$covariates))
    expect_equal(nrow(result$mean_model$covariates), 1)
    expect_equal(ncol(result$mean_model$covariates), 2)
    expect_true(all(result$mean_model$covariates == 1L))

    # valid model: spatial order can be given as binary matrix
    model_orders_mean <- list(intercept = "homogeneous", past_obs = matrix(1, nrow = 3, ncol = 1), 
                    past_mean = matrix(1, nrow = 2, ncol = 1))
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))
    model_orders_mean <- list(intercept = "homogeneous", past_obs = matrix(c(0, 1), nrow = 2, ncol = 1), 
                    past_mean = matrix(1, nrow = 2, ncol = 1))
    parameter2 <- params_mean
    parameter2$past_obs <- matrix(c(0, 0.3), nrow = 2)
    result <- dglmstarma.sim(n_obs, parameter2, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))
})



test_that("covariates dimensions are validated", {
    set.seed(42)
    n_obs <- 200L
    W <- generateW("rectangle", 100, 2, 10)
    model_orders_mean <- list(intercept = "homogeneous", 
                              past_obs = 2, past_mean = 1, 
                              covariates = c(0, 0))
    model_orders_dispersion <- list(intercept = "homogeneous", 
                                    past_obs = 1, 
                                    covariates = c(0, 0))
    covariates_mean <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(n_obs))),
                            location = TimeConstant(rnorm(100, sd = 0.81)))

    covariates_dispersion <- list(season = SpatialConstant(sin(2 * pi / 24 * seq(n_obs))),
                                  location = TimeConstant(runif(100)))

    params_mean <- list(intercept = 0.6, 
                        past_mean = matrix(c(0.2, 0.1), nrow = 2), 
                        past_obs = matrix(c(0.2, 0.1, 0.05), nrow = 3), 
                        covariates = matrix(c(0.9, 0.2), ncol = 2))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = matrix(c(0.5, 0.2), nrow = 2), 
                              covariates = matrix(c(0.1, 0.75), ncol = 2))
    family <- vnormal(copula = "frank", copula_param = 2)
    # correct covariates dimensions
    covariates_invalid <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                               location = TimeConstant(rnorm(100, sd = 0.81)))
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))
    # incorrect covariates dimensions (TimeConstant length)
    covariates_invalid2 <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                                location = TimeConstant(rnorm(50, sd = 0.81)))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_invalid2, 
                   dispersion_covariates = covariates_dispersion))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_invalid2))
    # incorrect covariates dimensions (SpatialConstant length)
    covariates_invalid3 <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(150))),
                                location = TimeConstant(rnorm(100, sd = 0.81)))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_invalid3, 
                   dispersion_covariates = covariates_dispersion))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_invalid3))
    # incorrect covariate dimensions (SpatialConstant is too long, warning case)
    covariates_invalid4 <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(250))),
                                location = TimeConstant(rnorm(100, sd = 0.81)))
    expect_warning(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_invalid4, 
                   dispersion_covariates = covariates_dispersion))
    expect_warning(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_invalid4))
    # incorrect covariate dimensions (matrix covariates with wrong nrow)
    covariates_invalid5 <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                                matr = matrix(rnorm(50 * 200), ncol = 200))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_invalid5, 
                   dispersion_covariates = covariates_dispersion))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_invalid5))
    # incorrect covariate dimensions (matrix covariates with wrong ncol, too few columns)
    covariates_invalid6 <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                                matr = matrix(rnorm(100 * 50), ncol = 50))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_invalid6, 
                   dispersion_covariates = covariates_dispersion))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_invalid6))
    # incorrect covariate dimensions (matrix covariates with wrong ncol, too many columns, warning case)
    covariates_invalid7 <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                                matr = matrix(rnorm(100 * 250), ncol = 250))
    expect_warning(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_invalid7, 
                   dispersion_covariates = covariates_dispersion))
    expect_warning(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_invalid7))
    # covariates can be named or unnamed
    covariates_unnamed <- list(SpatialConstant(sin(2 * pi / 12 * seq(200))),
                               TimeConstant(rnorm(100, sd = 0.81)))
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_unnamed, 
                   dispersion_covariates = covariates_unnamed)
    expect_true(is.list(result))
    covariates2 <- list(SpatialConstant(sin(2 * pi / 12 * seq(200))),
                        location = TimeConstant(rnorm(100, sd = 0.81)))
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates2, 
                   dispersion_covariates = covariates2)
    expect_true(is.list(result))
})


test_that("family argument is validated", {
    set.seed(42)
    n_obs <- 200L
    W <- generateW("rectangle", 100, 2, 10)
    model_orders_mean <- list(intercept = "homogeneous", 
                              past_obs = 2, past_mean = 1, 
                              covariates = c(0, 0))
    model_orders_dispersion <- list(intercept = "homogeneous", 
                                    past_obs = 1, 
                                    covariates = c(0, 0))
    covariates_mean <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(n_obs))),
                            location = TimeConstant(rnorm(100, sd = 0.81)))

    covariates_dispersion <- list(season = SpatialConstant(sin(2 * pi / 24 * seq(n_obs))),
                                  location = TimeConstant(runif(100)))

    params_mean <- list(intercept = 0.6, 
                        past_mean = matrix(c(0.2, 0.1), nrow = 2), 
                        past_obs = matrix(c(0.2, 0.1, 0.05), nrow = 3), 
                        covariates = matrix(c(0.9, 0.2), ncol = 2))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = matrix(c(0.5, 0.2), nrow = 2), 
                              covariates = matrix(c(0.1, 0.75), ncol = 2))
    family <- vnormal(copula = "frank", copula_param = 2)
    # family not specified
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = NULL, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    # family is not of class "stfamily"
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = "not_a_family", 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = poisson("log"), 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = list(), 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
})



test_that("parameter matches model orders", {
    set.seed(42)
    n_obs <- 200L
    W <- generateW("rectangle", 100, 2, 10)
    model_orders_mean <- list(intercept = "homogeneous", 
                              past_obs = 2, past_mean = 1, 
                              covariates = c(0, 0))
    model_orders_dispersion <- list(intercept = "homogeneous", 
                                    past_obs = 1, 
                                    covariates = c(0, 0))
    covariates_mean <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(n_obs))),
                            location = TimeConstant(rnorm(100, sd = 0.81)))

    covariates_dispersion <- list(season = SpatialConstant(sin(2 * pi / 24 * seq(n_obs))),
                                  location = TimeConstant(runif(100)))

    
    family <- vnormal(copula = "frank", copula_param = 2)

    # parameters can be submitted in matrix (always) in the list, intercept as vector
    params_mean <- list(intercept = 0.6, 
                        past_mean = matrix(c(0.2, 0.1), nrow = 2), 
                        past_obs = matrix(c(0.2, 0.1, 0.05), nrow = 3), 
                        covariates = matrix(c(0.9, 0.2), ncol = 2))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = matrix(c(0.5, 0.2), nrow = 2), 
                              covariates = matrix(c(0.1, 0.75), ncol = 2))
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))
    # parameters can be given as numeric vectors (if all spatial orders are 0 or only one time lag)
    params_mean <- list(intercept = 0.5, 
                      past_obs = c(0.3, 0.2, 0.05),
                      past_mean = c(0.1, 0.05),
                      covariates = c(0.75, 0.5))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = c(0.5, 0.2),
                              covariates = c(0.1, 0.75))
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))

    model_orders2 <- list(intercept = "homogeneous", past_obs = 2, past_mean = c(1, 0))
    params_mean <- list(intercept = 0.5, 
                      past_obs = c(0.3, 0.2, 0.05),
                      past_mean = c(0.1, 0.05, 0.1),
                      covariates = c(0.75, 0.5))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))

    # error if parameters for some model element are missing
    params_mean <- list(intercept = 0.5, 
                      past_mean = c(0.1, 0.05),
                      covariates = c(0.75, 0.5))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    params_dispersion <- list(intercept = 0.5, 
                      past_obs = c(0.3, 0.2, 0.05),
                      covariates = c(0.75, 0.5))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    params_mean <- list(intercept = 0.5, 
                      past_obs = c(0.3, 0.2, 0.05),
                      past_mean = c(0.1, 0.05))
    model_orders2 <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1, covariates = c(1, 1))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders2, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))

    # error if big mismatch between parameter and model order
    parameter <- list(intercept = c(0.5, 0.3), 
                      past_obs = c(0.3, 0.2, 0.05), 
                      past_mean = c(0.1, 0.05),
                      covariates = c(0.75, 0.5))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))

    params_mean <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2), nrow = 2),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    params_dispersion <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2, 0.05, 0.1), nrow = 4),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    params_mean <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2, 0.05, 0.1, 0.1, 0.1), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))

    params_mean <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1), nrow = 1),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    params_mean <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5, 0.1), ncol = 3))
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    
    # warning if non-zero parameters for spatial order excluded in model (only if spatial order is given as binary matrix)
    model_orders3 <- list(intercept = "homogeneous", past_obs = matrix(c(0, 1, 1), nrow = 3, ncol = 1), 
                    past_mean = matrix(1, nrow = 2, ncol = 1),
                    covariates = c(0, 0))
    parameter <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = matrix(c(0.5, 0.2), nrow = 2), 
                              covariates = matrix(c(0.1, 0.75), ncol = 2))
    expect_warning(dglmstarma.sim(n_obs, parameter, params_dispersion, model_orders3, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
})



test_that("errors and warnings are thrown if invalid parameters are provided", {
    set.seed(42)
    n_obs <- 200L
    W <- generateW("rectangle", 100, 2, 10)
    model_orders_mean <- list(intercept = "homogeneous", 
                              past_obs = 2, past_mean = 1, 
                              covariates = c(0, 0))
    model_orders_dispersion <- list(intercept = "homogeneous", 
                                    past_obs = 1, 
                                    covariates = c(0, 0))
    covariates_mean_positive <- list(season = SpatialConstant(abs(sin(2 * pi / 12 * seq(n_obs)))),
                            location = TimeConstant(abs(rnorm(100, sd = 0.81))))
    covariates_mean_negative <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(n_obs))),
                            location = TimeConstant(rnorm(100, sd = 0.81)))

    covariates_dispersion_positive <- list(season = SpatialConstant(abs(sin(2 * pi / 24 * seq(n_obs)))),
                                  location = TimeConstant(abs(runif(100))))
    covariates_dispersion_negative <- list(season = SpatialConstant(sin(2 * pi / 24 * seq(n_obs))),
                                  location = TimeConstant(runif(100)))
    family <- vnormal(copula = "frank", copula_param = 2)


    # negative parameters are provided for family that allows only non-negative parameters
    params_mean <- list(intercept = 0.6, 
                        past_mean = matrix(c(0.2, 0.1), nrow = 2), 
                        past_obs = matrix(c(0.2, -0.1, 0.05), nrow = 3), 
                        covariates = matrix(c(0.9, 0.2), ncol = 2))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = matrix(c(0.5, -0.2), nrow = 2), 
                              covariates = matrix(c(0.1, 0.75), ncol = 2))
    family <- vpoisson("identity", copula = "frank", copula_param = 2)

    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean_positive, 
                   dispersion_covariates = covariates_dispersion_positive))
    family <- vpoisson("log", copula = "frank", copula_param = 2)               
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, dispersion_link = "identity",
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean_positive, 
                   dispersion_covariates = covariates_dispersion_positive))

    # negative covariate values when family requires non-negative parameters
    params_mean <- list(intercept = 0.6, 
                        past_mean = matrix(c(0.2, 0.1), nrow = 2), 
                        past_obs = matrix(c(0.2, 0.1, 0.05), nrow = 3), 
                        covariates = matrix(c(0.9, 0.2), ncol = 2))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = matrix(c(0.5, 0.2), nrow = 2), 
                              covariates = matrix(c(0.1, 0.75), ncol = 2))

    family <- vpoisson("identity", copula = "frank", copula_param = 2)

    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean_negative, 
                   dispersion_covariates = covariates_dispersion_positive))
    family <- vpoisson("log", copula = "frank", copula_param = 2)               
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, dispersion_link = "identity",
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean_positive, 
                   dispersion_covariates = covariates_dispersion_negative))

    # NA parameter
    params_mean <- list(intercept = 0.6, 
                        past_mean = matrix(c(0.2, 0.1), nrow = 2), 
                        past_obs = matrix(c(0.2, NA, 0.05), nrow = 3), 
                        covariates = matrix(c(0.9, 0.2), ncol = 2))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = matrix(c(0.5, NA), nrow = 2), 
                              covariates = matrix(c(0.1, 0.75), ncol = 2))

    family <- vpoisson("log", copula = "frank", copula_param = 2)

    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean_negative, 
                   dispersion_covariates = covariates_dispersion_negative))            

    # warning if parameters might lead to exploding values (sum of ar coefficients >= 1)
    params_mean <- list(intercept = 0.6, 
                        past_mean = matrix(c(0.2, 0.1), nrow = 2), 
                        past_obs = matrix(c(0.2, 0.4, 0.1), nrow = 3), 
                        covariates = matrix(c(0.9, 0.2), ncol = 2))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = matrix(c(0.5, 0.2), nrow = 2), 
                              covariates = matrix(c(0.1, 0.75), ncol = 2))

    family <- vpoisson("log", copula = "frank", copula_param = 2)

    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean_negative, 
                   dispersion_covariates = covariates_dispersion_negative))  
})






















test_that("generated values match distribution", {
    set.seed(42)
    n_obs <- 200L
    W <- generateW("rectangle", 100, 2, 10)
    model_orders_mean <- list(intercept = "homogeneous", 
                              past_obs = 2, past_mean = 1, 
                              covariates = c(0, 0))
    model_orders_dispersion <- list(intercept = "homogeneous", 
                                    past_obs = 1, 
                                    covariates = c(0, 0))
    covariates_mean <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(n_obs))),
                            location = TimeConstant(rnorm(100, sd = 0.81)))

    covariates_dispersion <- list(season = SpatialConstant(sin(2 * pi / 24 * seq(n_obs))),
                                  location = TimeConstant(runif(100)))
    params_mean <- list(intercept = 0.6, 
                        past_mean = matrix(c(0.2, 0.1), nrow = 2), 
                        past_obs = matrix(c(0.2, 0.1, 0.05), nrow = 3), 
                        covariates = matrix(c(0.9, 0.2), ncol = 2))
    params_dispersion <- list(intercept = 0.5, 
                              past_obs = matrix(c(0.5, 0.2), nrow = 2), 
                              covariates = matrix(c(0.1, 0.75), ncol = 2))
    

    # poisson
    fam <- vpoisson("log")
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = fam, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    # quasipoisson
    fam <- vquasipoisson("log")
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = fam, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))
    expect_true(any(result$observations > 0))
    expect_true(all(result$observations >= 0))
    expect_true(all(result$observations == floor(result$observations)))
    # negative binomial
    fam <- vnegative.binomial("log")
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = fam, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))
    expect_true(any(result$observations > 0))
    expect_true(all(result$observations >= 0))
    expect_true(all(result$observations == floor(result$observations)))
    # binomial
    fam <- vbinomial("softclipping", size = 10)
    expect_error(dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = fam, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion))
    # quasibinomial
    fam <- vquasibinomial("softclipping", size = 10)
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = fam, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))
    expect_true(any(result$observations > 0))
    expect_true(all(result$observations >= 0))
    expect_true(all(result$observations <= 10))
    expect_true(all(result$observations == floor(result$observations)))
    # gamma
    fam <- vgamma("log")
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = fam, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))
    expect_true(all(result$observations > 0))
    # inverse gaussian
    fam <- vinverse.gaussian("log")
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = fam, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))
    expect_true(all(result$observations > 0))
    # normal
    fam <- vnormal("identity")
    result <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = fam, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)
    expect_true(is.list(result))
    expect_true(all(result$observations != 0))
})

