# glmstarma_sim:
test_that("ntime is handled correctly", {
    set.seed(42)
    W <- generateW("rectangle", 100, 2, 10)
    fam <- vpoisson("log")
    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1)
    parameter <- list(intercept = 0.5, past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    covariates <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                       location = TimeConstant(rnorm(100, sd = 0.81)))
    # correct n_obs input
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))
    expect_type(result, "list")
    expect_named(
        result,
        c("observations", "link_values", "model", "parameters")
    )
    expect_length(result, 4)

    expect_type(result$observations, "double")
    expect_true(is.matrix(result$observations))
    expect_equal(dim(result$observations), c(100, 200))

    expect_type(result$link_values, "double")
    expect_true(is.matrix(result$link_values))
    expect_equal(dim(result$link_values), c(100, 200))

    # incorrect n_obs input
    expect_error(glmstarma.sim(NA, parameter, model_orders, W, covariates, family = fam))
    expect_error(glmstarma.sim(NULL, parameter, model_orders, W, covariates, family = fam))
    expect_error(glmstarma.sim(Inf, parameter, model_orders, W, covariates, family = fam))
    expect_error(glmstarma.sim(-10, parameter, model_orders, W, covariates, family = fam))
    expect_error(glmstarma.sim(0, parameter, model_orders, W, covariates, family = fam))
    expect_error(glmstarma.sim(10.5, parameter, model_orders, W, covariates, family = fam))
})

test_that("n_start is handled correctly", {
    set.seed(42)
    W <- generateW("rectangle", 100, 2, 10)
    fam <- vpoisson("log")
    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1)
    parameter <- list(intercept = 0.5, past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    covariates <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                       location = TimeConstant(rnorm(100, sd = 0.81)))
    # correct n_obs input
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam,
                            n_start = 100L, control = list(return_burn_in = TRUE))
    expect_type(result$observations, "double")
    expect_true(is.matrix(result$observations))
    expect_equal(dim(result$observations), c(100, 300))

    expect_type(result$link_values, "double")
    expect_true(is.matrix(result$link_values))
    expect_equal(dim(result$link_values), c(100, 300))
    # incorrect nstart input
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam,
                               n_start = -10L))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam,
                               n_start = 0L))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam,
                               n_start = 10.5))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam,
                               n_start = NA))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam,
                               n_start = Inf))
    # n_start is too small
    model_orders$past_obs_time_lags <- 10L
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam,
                               n_start = 5L))
})

test_that("control gets checked", {
    set.seed(42)
    W <- generateW("rectangle", 100, 2, 10)
    fam <- vpoisson("log")
    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1)
    parameter <- list(intercept = 0.5, past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    covariates <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                       location = TimeConstant(rnorm(100, sd = 0.81)))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam,
                               control = "not_a_list"))
})


test_that("wlist arguments are validated", {
    set.seed(42)
    W <- generateW("rectangle", 100, 2, 10)
    fam <- vpoisson("log")
    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1)
    parameter <- list(intercept = 0.5, past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    covariates <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                       location = TimeConstant(rnorm(100, sd = 0.81)))
    # too few matrices wlist
    expect_error(glmstarma.sim(200, parameter, model_orders, W[1:2], covariates, family = fam))
    # wlist is not a list
    expect_error(glmstarma.sim(200, parameter, model_orders, "not_a_list", covariates, family = fam,
                            n_start = 100L, control = list(return_burn_in = TRUE)))
    expect_error(glmstarma.sim(200, parameter, model_orders, W[[1]], covariates, family = fam))
    # wlist not specified
    expect_error(glmstarma.sim(200, parameter, model_orders, wlist_past_mean = W, covariates = covariates, family = fam))
    expect_error(glmstarma.sim(200, parameter, model_orders, wlist_covariates = W, covariates = covariates, family = fam))
    # wlist_past_mean has not enough matrices
    expect_error(glmstarma.sim(200, parameter, model_orders, wlist = W, covariates, wlist_past_mean = W[1], family = fam))
    # wlist_covariates has not enough matrices
    model_orders2 <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1, covariates = c(1, 1))
    parameter2 <- list(intercept = 0.5, past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5, 0.1, 0.1), ncol = 2))
    expect_error(glmstarma.sim(200, parameter2, model_orders2, wlist = W, covariates, wlist_covariates = W[1], family = fam))
    # elements of wlist are not numeric matrices
    W_invalid <- W
    W_invalid[[1]] <- "not_a_matrix"
    expect_error(glmstarma.sim(200, parameter, model_orders, W_invalid, covariates, family = fam))
    # elements of wlist are not of the same dimension
    W_invalid2 <- W
    W_invalid2[[1]] <- diag(50)
    expect_error(glmstarma.sim(200, parameter, model_orders, W_invalid2, covariates, family = fam))
    # correct case: sparse matrices
    W_sparse <- lapply(W, function(mat) as(mat, "dgCMatrix"))
    result <- glmstarma.sim(200, parameter, model_orders, W_sparse, covariates, family = fam)
    expect_true(is.list(result))
    # correct case: mixed dense and sparse matrices
    W_mixed <- W
    W_mixed[[1]] <- as(W_mixed[[1]], "dgCMatrix")
    result <- glmstarma.sim(200, parameter, model_orders, W_mixed, covariates, family = fam)
    expect_true(is.list(result))
})


test_that("model is handled correctly", {
    set.seed(42)
    W <- generateW("rectangle", 100, 2, 10)
    fam <- vpoisson("log")
    parameter <- list(intercept = 0.5, past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    covariates <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                       location = TimeConstant(rnorm(100, sd = 0.81)))
    # invalid model: negative time lags or zero or NA, Inf
    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1,
                         past_obs_time_lags = -1)
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))

    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1,
                         past_mean_time_lags = -1)
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))

    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1,
                         past_obs_time_lags = NA)
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))

    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1,
                         past_mean_time_lags = NA)
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))

    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1,
                         past_obs_time_lags = Inf)
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))

    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1,
                         past_mean_time_lags = Inf)
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))

    # invalid model: no regression on past observations when past_mean is specified
    model_orders <- list(intercept = "homogeneous", past_mean = 1)
    parameter2 <- parameter
    parameter2$past_obs <- NULL
    expect_error(glmstarma.sim(200, parameter2, model_orders, W, covariates, family = fam))

    # invalid model: intercept is not "homogeneous" or "inhomogeneous"
    model_orders <- list(intercept = "invalid_value", past_obs = 2, past_mean = 1)
    expect_error(glmstarma.sim(200, parameter2, model_orders, W, covariates, family = fam))

    # valid model: no regression either on past observations or past mean
    model_orders <- list(intercept = "homogeneous")
    parameter2 <- parameter
    parameter2$past_obs <- NULL
    parameter2$past_mean <- NULL
    result <- glmstarma.sim(200, parameter2, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))

    # valid model: only regression on past observations
    model_orders <- list(intercept = "homogeneous", past_obs = 2)
    parameter2 <- parameter
    parameter2$past_mean <- NULL
    result <- glmstarma.sim(200, parameter2, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))

    # valid model: no covariate orders specified (spatial orders default to 0)
    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1)
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))
    expect_true(is.matrix(result$model$covariates))
    expect_equal(nrow(result$model$covariates), 1)
    expect_equal(ncol(result$model$covariates), 2)
    expect_true(all(result$model$covariates == 1L))

    # valid model: spatial order can be given as binary matrix
    model_orders <- list(intercept = "homogeneous", past_obs = matrix(1, nrow = 3, ncol = 1), 
                    past_mean = matrix(1, nrow = 2, ncol = 1))
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))
    model_orders <- list(intercept = "homogeneous", past_obs = matrix(c(0, 1), nrow = 2, ncol = 1), 
                    past_mean = matrix(1, nrow = 2, ncol = 1))
    parameter2 <- parameter
    parameter2$past_obs <- matrix(c(0, 0.3), nrow = 2)
    result <- glmstarma.sim(200, parameter2, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))
})


test_that("covariates dimensions are validated", {
    set.seed(42)
    W <- generateW("rectangle", 100, 2, 10)
    fam <- vpoisson("log")
    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1)
    parameter <- list(intercept = 0.5, past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    # correct covariates dimensions
    covariates_invalid <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                               location = TimeConstant(rnorm(100, sd = 0.81)))
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates_invalid, family = fam)
    expect_true(is.list(result))
    # incorrect covariates dimensions (TimeConstant length)
    covariates_invalid2 <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                                location = TimeConstant(rnorm(50, sd = 0.81)))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates_invalid2, family = fam))
    # incorrect covariates dimensions (SpatialConstant length)
    covariates_invalid3 <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(150))),
                                location = TimeConstant(rnorm(100, sd = 0.81)))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates_invalid3, family = fam))
    # incorrect covariate dimensions (SpatialConstant is too long, warning case)
    covariates_invalid4 <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(250))),
                                location = TimeConstant(rnorm(100, sd = 0.81)))
    expect_warning(glmstarma.sim(200, parameter, model_orders, W, covariates_invalid4, family = fam))

    # incorrect covariate dimensions (matrix covariates with wrong nrow)
    covariates_invalid5 <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                                matr = matrix(rnorm(50 * 200), ncol = 200))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates_invalid5, family = fam))
    # incorrect covariate dimensions (matrix covariates with wrong ncol, too few columns)
    covariates_invalid6 <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                                matr = matrix(rnorm(100 * 50), ncol = 50))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates_invalid6, family = fam))
    # incorrect covariate dimensions (matrix covariates with wrong ncol, too many columns, warning case)
    covariates_invalid7 <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                                matr = matrix(rnorm(100 * 250), ncol = 250))
    expect_warning(glmstarma.sim(200, parameter, model_orders, W, covariates_invalid7, family = fam))
    # covariates can be named or unnamed
    covariates_unnamed <- list(SpatialConstant(sin(2 * pi / 12 * seq(200))),
                               TimeConstant(rnorm(100, sd = 0.81)))
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates_unnamed, family = fam)
    expect_true(is.list(result))
    covariates2 <- list(SpatialConstant(sin(2 * pi / 12 * seq(200))),
                        location = TimeConstant(rnorm(100, sd = 0.81)))
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates2, family = fam)
    expect_true(is.list(result))
})


test_that("family argument is validated", {
    set.seed(42)
    W <- generateW("rectangle", 100, 2, 10)
    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1)
    parameter <- list(intercept = 0.5, past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    covariates <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                       location = TimeConstant(rnorm(100, sd = 0.81)))
    # family is not specified
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates))
    # family is not of class "stfamily"
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = "not_a_family"))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = poisson("log")))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family =  list()))
    # correct family
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = vpoisson("log"))
    expect_true(is.list(result))
})


test_that("parameter matches model orders", {
    set.seed(42)
    W <- generateW("rectangle", 100, 2, 10)
    fam <- vpoisson("log")
    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1)
    covariates <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                       location = TimeConstant(rnorm(100, sd = 0.81)))
    # parameters can be submitted in matrix (always) in the list, intercept as vector
    parameter <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))
    # parameters can be given as numeric vectors (if all spatial orders are 0 or only one time lag)
    parameter <- list(intercept = 0.5, 
                      past_obs = c(0.3, 0.2, 0.05),
                      past_mean = c(0.1, 0.05),
                      covariates = c(0.75, 0.5))
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))

    # model_orders2 <- list(intercept = "homogeneous", past_obs = 2, past_mean = c(1, 0))
    # parameter <- list(intercept = 0.5, 
    #                   past_obs = c(0.3, 0.2, 0.05),
    #                   past_mean = c(0.1, 0.05, 0.1),
    #                   covariates = c(0.75, 0.5))
    # expect_error(glmstarma.sim(200, parameter, model_orders2, W, covariates, family = fam))

    # error if parameters for some model element are missing
    parameter <- list(intercept = 0.5, 
                      past_mean = c(0.1, 0.05),
                      covariates = c(0.75, 0.5))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))
    parameter <- list(intercept = 0.5, 
                      past_obs = c(0.3, 0.2, 0.05),
                      covariates = c(0.75, 0.5))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))
    parameter <- list(intercept = 0.5, 
                      past_obs = c(0.3, 0.2, 0.05),
                      past_mean = c(0.1, 0.05))
    model_orders2 <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1, covariates = c(1, 1))
    expect_error(glmstarma.sim(200, parameter, model_orders2, W, covariates, family = fam))

    # error if big mismatch between parameter and model order
    parameter <- list(intercept = c(0.5, 0.3), 
                      past_obs = c(0.3, 0.2, 0.05), 
                      past_mean = c(0.1, 0.05),
                      covariates = c(0.75, 0.5))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))

    parameter <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2), nrow = 2),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))
    parameter <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2, 0.05, 0.1), nrow = 4),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))

    parameter <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2, 0.05, 0.1, 0.1, 0.1), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))

    parameter <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1), nrow = 1),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))

    parameter <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5, 0.1), ncol = 3))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))
    
    # warning if non-zero parameters for spatial order excluded in model (only if spatial order is given as binary matrix)
    model_orders3 <- list(intercept = "homogeneous", past_obs = matrix(c(0, 1, 1), nrow = 3, ncol = 1), 
                    past_mean = matrix(1, nrow = 2, ncol = 1))
    parameter <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    expect_warning(glmstarma.sim(200, parameter, model_orders3, W, covariates, family = fam))
})

test_that("errors and warnings are thrown if invalid parameters are provided", {
    set.seed(42)
    W <- generateW("rectangle", 100, 2, 10)
    fam <- vpoisson("identity")
    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1)
    covariates <- list(season = SpatialConstant(abs(sin(2 * pi / 12 * seq(200)))),
                       location = TimeConstant(abs(rnorm(100, sd = 0.81))))
    covariates_negative <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                                location = TimeConstant(runif(100, min = -1, max = 0)))
    # negative parameters are provided for family that requires non-negative parameters
    parameter <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, -0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))
    # negative covariates values when family requires non-negative parameters
    parameter2 <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    expect_error(glmstarma.sim(200, parameter2, model_orders, W, covariates_negative, family = fam))
    # NA parameter
    parameter <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, NA, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    expect_error(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))
    # warning if parameters might lead to exploding values (sum of ar coefficients >= 1)
    parameter <- list(intercept = 0.5, 
                      past_obs = matrix(c(0.3, 0.2, 0.1), nrow = 3),
                      past_mean = matrix(c(0.2, 0.2), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    expect_warning(glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam))
})


test_that("generated values match distribution", {
    set.seed(42)
    W <- generateW("rectangle", 100, 2, 10)
    model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1)
    parameter <- list(intercept = 0.5, past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3),
                      past_mean = matrix(c(0.1, 0.05), nrow = 2),
                      covariates = matrix(c(0.75, 0.5), ncol = 2))
    covariates <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(200))),
                       location = TimeConstant(rnorm(100, sd = 0.81)))                   
    # poisson
    fam <- vpoisson("log")
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))
    expect_true(any(result$observations > 0))
    expect_true(all(result$observations >= 0))
    expect_true(all(result$observations == floor(result$observations)))
    expect_equal(mean((exp(result$link_values) - result$observations)^2 / exp(result$link_values)), 1, tolerance = 0.1)
    # quasipoisson
    fam <- vquasipoisson("log", dispersion = 2)
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))
    expect_true(any(result$observations > 0))
    expect_true(all(result$observations >= 0))
    expect_true(all(result$observations == floor(result$observations)))
    expect_equal(mean((exp(result$link_values) - result$observations)^2 / exp(result$link_values)), 2, tolerance = 0.2)
    # negative binomial
    fam <- vnegative.binomial("log", dispersion = 2)
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))
    expect_true(any(result$observations > 0))
    expect_true(all(result$observations >= 0))
    expect_true(all(result$observations == floor(result$observations)))
    expect_equal(mean((exp(result$link_values) - result$observations)^2 / (exp(result$link_values) + 2 * exp(result$link_values)^2)), 1, tolerance = 0.2)
    # binomial
    fam <- vbinomial("softclipping", size = 10)
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))
    expect_true(any(result$observations > 0))
    expect_true(all(result$observations >= 0))
    expect_true(all(result$observations <= 10))
    expect_true(all(result$observations == floor(result$observations)))
    expectation <- 10 * log((1 + exp(result$link_values)) / (1 + exp(result$link_values - 1)))
    expect_equal(mean((expectation - result$observations)^2 / (expectation * (1 - expectation / 10))), 1, tolerance = 0.2)
    # quasibinomial
    fam <- vquasibinomial("softclipping", size = 10, dispersion = 2)
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))
    expect_true(any(result$observations > 0))
    expect_true(all(result$observations >= 0))
    expect_true(all(result$observations <= 10))
    expect_true(all(result$observations == floor(result$observations)))
    expectation <- 10 * log((1 + exp(result$link_values)) / (1 + exp(result$link_values - 1)))
    expect_equal(mean((expectation - result$observations)^2 / (expectation * (1 - expectation / 10))), 2, tolerance = 0.2)
    # gamma
    fam <- vgamma("log", dispersion = 2)
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))
    expect_true(all(result$observations > 0))
    expectation <- exp(result$link_values)
    expect_equal(mean((expectation - result$observations)^2 / (expectation)^2), 2, tolerance = 0.2)
    # inverse gaussian
    fam <- vinverse.gaussian("log", dispersion = 2)
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))
    expect_true(all(result$observations > 0))
    expectation <- exp(result$link_values)
    expect_equal(mean((expectation - result$observations)^2 / (expectation)^3), 2, tolerance = 0.2)
    # normal
    fam <- vnormal("identity", dispersion = 2)
    result <- glmstarma.sim(200, parameter, model_orders, W, covariates, family = fam)
    expect_true(is.list(result))
    expect_true(all(result$observations != 0))
    expect_equal(mean((result$link_values - result$observations)^2), 2, tolerance = 0.2)
})
