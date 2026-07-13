
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
    dat <- load_data("chickenpox", directory = tempdir())
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary
    covariates <- list(population = population_hungary, 
                       season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
                       season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))

    # Parameter initialization
    set.seed(123)
    result <- dglmstarma(chickenpox, list(past_obs = 1), list(past_obs = 1), mean_family = vquasipoisson("log"), 
          dispersion_link = "log", W_hungary, 
          control = list(parameter_init = list(intercept = 1.5, 
                                               past_obs = matrix(c(0.5, 0.1), nrow = 2),
                                               covariates = matrix(c(0.1, 0.1, 0.5), nrow = 1)),
                        parameter_init_dispersion = list(intercept = 1.5, 
                                                        past_obs = matrix(c(0.5, 0.1), nrow = 2)),                      
                        print_progress = FALSE))
    

    result <- dglmstarma(chickenpox, list(past_obs = 1), list(past_obs = 1), mean_family = vquasipoisson("log"), 
          dispersion_link = "log", W_hungary, 
          control = list(parameter_init = list(intercept = 1.5, 
                                               past_obs = matrix(c(0.5, 0.1), nrow = 2),
                                               covariates = matrix(c(0.1, 0.1, 0.5), nrow = 1)),
                        parameter_init_dispersion = list(intercept = 1.5, 
                                                        past_obs = matrix(c(0.5, 0.1), nrow = 2)),                      
                        print_progress = FALSE))

    # initialization of feedback term
    result <- dglmstarma(chickenpox, list(past_obs = 1, past_mean = 1), list(past_obs = 1), mean_family = vquasipoisson("log"), 
          dispersion_link = "log", W_hungary, 
          control = list(init_link = matrix(runif(20, min = 1, max = 10), nrow = 20)))

    result <- dglmstarma(chickenpox, list(past_obs = 1), list(past_obs = 1, past_mean = 1), mean_family = vquasipoisson("log"), 
          dispersion_link = "log", W_hungary, 
          control = list(init_dispersion = matrix(runif(20, min = 1, max = 10), nrow = 20)))

})



















