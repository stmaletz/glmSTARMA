
test_that("coef function works", {
    dat <- load_data("chickenpox")
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary

    model_autoregressive <- list(past_obs = rep(1, 7))
    fit <- glmstarma(chickenpox, model_autoregressive, W_hungary, family = vpoisson("log"),
                    covariates = list(population = population_hungary))
    x <- coef(fit)
    expect_true(is.numeric(x))
    expect_equal(length(x), 16)
    x <- coef(fit, asList = TRUE)
    expect_true(is.list(x))
    expect_named(x, c("intercept", "past_obs", "past_mean", "covariates"))
    
    mean_model <- list(past_obs = rep(1, 7))
    dispersion_model <- list(past_obs = 1)
    fit2 <- dglmstarma(chickenpox, mean_model, dispersion_model, mean_family = vquasipoisson("log"),
                   dispersion_link = "log",
                   wlist = W_hungary, 
                   mean_covariates = list(population = population_hungary))
    x <- coef(fit2)
    expect_true(is.numeric(x))
    expect_equal(length(x), 19)
    x <- coef(fit2, asList = TRUE)
    expect_true(is.list(x))
    expect_named(x, c("mean", "dispersion"))
    expect_named(x$mean, c("intercept", "past_obs", "past_mean", "covariates"))
    expect_named(x$dispersion, c("intercept", "past_obs", "past_mean", "covariates"))
    delete_glmSTARMA_data("chickenpox")  # Clean up cached data
})  

test_that("vcov function works", {
    dat <- load_data("chickenpox")
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary

    model_autoregressive <- list(past_obs = rep(1, 7))
    fit <- glmstarma(chickenpox, model_autoregressive, W_hungary, family = vpoisson("log"),
                    covariates = list(population = population_hungary))
    x <- vcov(fit)
    expect_true(is.matrix(x))
    expect_true(is.numeric(x))
    expect_equal(dim(x), c(16, 16))
    expect_true(all(eigen(x)$values > 0))

    mean_model <- list(past_obs = rep(1, 7))
    dispersion_model <- list(past_obs = 1)
    fit2 <- dglmstarma(chickenpox, mean_model, dispersion_model, mean_family = vquasipoisson("log"),
                   dispersion_link = "log",
                   wlist = W_hungary, 
                   mean_covariates = list(population = population_hungary))
    x <- vcov(fit2, "mean")
    expect_true(is.matrix(x))
    expect_true(is.numeric(x))
    expect_equal(dim(x), c(16, 16))
    expect_true(all(eigen(x)$values > 0))
    y <- vcov(fit2, "dispersion")
    expect_true(is.matrix(y))
    expect_true(is.numeric(y))
    expect_equal(dim(y), c(3, 3))
    expect_true(all(eigen(y)$values > 0))

    z <- vcov(fit2, "both")
    expect_true(is.list(z))
    expect_named(z, c("mean", "dispersion"))
    expect_equal(z$mean, x)
    expect_equal(z$dispersion, y)
    delete_glmSTARMA_data("chickenpox")  # Clean up cached data
})

## rest will be added later