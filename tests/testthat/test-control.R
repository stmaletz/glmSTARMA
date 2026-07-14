
testthat::skip_on_cran()
test_that("different control parameters work for glmstarma", {
    dat <- load_data("chickenpox", directory = tempdir())
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary
    covariates <- list(population = population_hungary, 
                       season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                       season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                    covariates = covariates, family = vpoisson("log"),
                    control = list(parameter_init = list(intercept = 1.5, 
                                                         past_obs = matrix(c(0.5, 0.1), nrow = 2),
                                                         covariates = matrix(c(0.1, 0.1, 0.5), nrow = 1))))

    set.seed(123)
    result <- glmstarma(chickenpox, list(past_obs = 1, past_mean = 1), wlist = W_hungary, 
                    covariates = covariates, family = vpoisson("log"),
                    control = list(init_link = matrix(runif(20, min = 1, max = 10), nrow = 20)))                                                     

    # Constrained optimization
    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log"),
                        control = list(constrained = FALSE))
    expect_s3_class(result, "glmstarma")

    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log"),
                        control = list(constrain_method = "sum_of_absolutes"))
    expect_s3_class(result, "glmstarma")

    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log"),
                        control = list(constrain_method = "absolute_sum"))
    expect_s3_class(result, "glmstarma")

    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("softplus"),
                        control = list(constrain_method = "soft"))
    expect_s3_class(result, "glmstarma")

    # Optimization algorithm
    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        covariates = covariates, family = vpoisson("log"),
                        control = list(method = "optim"))   
    expect_s3_class(result, "glmstarma")

    result <- glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
                        family = vpoisson("identity"),
                        control = list(method = "optim"))                       
    expect_s3_class(result, "glmstarma")

})




testthat::skip_on_cran()
test_that("different control parameters work for dglmstarma", {
    # Simulate data
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
    sim <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
                   model_orders_dispersion, mean_family = family, 
                   wlist = W, pseudo_observations = "deviance", 
                   mean_covariates = covariates_mean, 
                   dispersion_covariates = covariates_dispersion)

    # Parameter initialization
    result <- dglmstarma(sim$observations, model_orders_mean, model_orders_dispersion, mean_family = vnormal(), 
                         dispersion_link = "log", W, mean_covariates = covariates_mean, dispersion_covariates = covariates_dispersion,
                         control = list(parameter_init = params_mean,
                                        parameter_init_dispersion = params_dispersion,
                                        print_progress = FALSE))
    expect_s3_class(result, "dglmstarma")

    # Initialization of feedback term
    result <- dglmstarma(sim$observations, model_orders_mean, model_orders_dispersion, mean_family = vnormal(), 
                     dispersion_link = "log", W, mean_covariates = covariates_mean, dispersion_covariates = covariates_dispersion,
                     control = list(init_link = sim$link_values[, 1, drop = FALSE],
                                    maxit = 100,
                                    print_progress = FALSE))
    expect_s3_class(result, "dglmstarma")

    result <- dglmstarma(sim$observations, model_orders_mean, list(past_obs = 1, past_mean = 0), mean_family = vnormal(), 
                     dispersion_link = "log", W, mean_covariates = covariates_mean, dispersion_covariates = covariates_dispersion,
                     control = list(init_link = sim$link_values[, 1, drop = FALSE],
                                    maxit = 100,
                                    print_progress = FALSE))
    expect_s3_class(result, "dglmstarma")
})



















