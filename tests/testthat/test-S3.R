testthat::skip_on_cran()
test_that("coef function works", {
    dat <- load_data("chickenpox", directory = tempdir())
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
})
 
testthat::skip_on_cran()
test_that("vcov function works", {
    dat <- load_data("chickenpox", directory = tempdir())
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

    dat <- load_data("chickenpox", directory = tempdir())
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary

    fit2 <- dglmstarma(chickenpox, list(past_obs = 1, past_mean = 1), list(intercept = "inhomogeneous"), mean_family = vquasipoisson("identity"),
                   dispersion_link = "log",
                   wlist = W_hungary, 
                   mean_covariates = list(population = population_hungary),
                   control = list(use_fast_if_const_dispersion = TRUE))
    z <- vcov(fit2, "dispersion")
    expect_true(is.matrix(z))
    expect_true(is.numeric(z))
    expect_equal(dim(z), c(20, 20))
    expect_true(all(eigen(z)$values > 0))

})

testthat::skip_on_cran()
test_that("information criteria work", {
    dat <- load_data("chickenpox", directory = tempdir())
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary

    model_autoregressive <- list(past_obs = rep(1, 7))
    fit <- glmstarma(chickenpox, model_autoregressive, W_hungary, family = vpoisson("log"),
                    covariates = list(population = population_hungary))
    # Test AIC, BIC, QIC, logLik with adjust = TRUE
    x <- AIC(fit)
    expect_true(is.numeric(x))
    x <- BIC(fit)
    expect_true(is.numeric(x))
    x <- QIC(fit)
    expect_true(is.numeric(x))
    x <- BIC(fit)
    expect_true(is.numeric(x))
    x <- logLik(fit)
    expect_true(is.numeric(x))

    # Test AIC, BIC, QIC, logLik with adjust = FALSE
    x <- AIC(fit, adjust = FALSE)
    expect_true(is.numeric(x))
    x <- BIC(fit, adjust = FALSE)
    expect_true(is.numeric(x))
    x <- QIC(fit, adjust = FALSE)
    expect_true(is.numeric(x))
    x <- BIC(fit, adjust = FALSE)
    expect_true(is.numeric(x))
    x <- logLik(fit, adjust = FALSE)
    expect_true(is.numeric(x))


    mean_model <- list(past_obs = rep(1, 7))
    dispersion_model <- list(past_obs = 1)
    fit2 <- dglmstarma(chickenpox, mean_model, dispersion_model, mean_family = vquasipoisson("log"),
                   dispersion_link = "log",
                   wlist = W_hungary, 
                   mean_covariates = list(population = population_hungary))

    # Test AIC, BIC, QIC, logLik with adjust = TRUE
    x <- AIC(fit2)
    expect_true(is.numeric(x))
    x <- BIC(fit2)
    expect_true(is.numeric(x))
    x <- QIC(fit2)
    expect_true(is.numeric(x))
    x <- logLik(fit2)
    expect_true(is.numeric(x))

    # Test AIC, BIC, QIC, logLik with adjust = FALSE
    x <- AIC(fit2, adjust = FALSE)
    expect_true(is.numeric(x))
    x <- BIC(fit2, adjust = FALSE)
    expect_true(is.numeric(x))
    x <- QIC(fit2, adjust = FALSE)
    expect_true(is.numeric(x))
    x <- BIC(fit2, adjust = FALSE)
    expect_true(is.numeric(x))
    x <- logLik(fit2, adjust = FALSE)
    expect_true(is.numeric(x))
})


testthat::skip_on_cran()
test_that("residuals function works", {
    dat <- load_data("chickenpox", directory = tempdir())
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary

    model_autoregressive <- list(past_obs = rep(1, 7))
    fit <- glmstarma(chickenpox, model_autoregressive, W_hungary, family = vpoisson("log"),
                    covariates = list(population = population_hungary))
    x <- residuals(fit)
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 515)
    x <- residuals(fit, drop_init = FALSE)
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 522)

    x <- residuals(fit, type = "pearson")
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 515)
    x <- residuals(fit, type = "pearson", drop_init = FALSE)
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 522)

    x <- residuals(fit, type = "deviance")
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 515)
    expect_true(all(x >= 0))
    x <- residuals(fit, type = "deviance", drop_init = FALSE)
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 522)
    expect_true(all(x >= 0))


    
    mean_model <- list(past_obs = rep(1, 7))
    dispersion_model <- list(past_obs = 1)
    fit2 <- dglmstarma(chickenpox, mean_model, dispersion_model, mean_family = vquasipoisson("log"),
                   dispersion_link = "log",
                   wlist = W_hungary, 
                   mean_covariates = list(population = population_hungary))
    x <- residuals(fit2)
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 515)
    x <- residuals(fit2, drop_init = FALSE)
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 522)

    x <- residuals(fit2, type = "pearson")
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 515)
    x <- residuals(fit2, type = "pearson", drop_init = FALSE)
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 522)

    x <- residuals(fit2, type = "deviance")
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 515)
    expect_true(all(x >= 0))
    x <- residuals(fit2, type = "deviance", drop_init = FALSE)
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 522)
    expect_true(all(x >= 0))
})



testthat::skip_on_cran()
test_that("fitted function works", {
    dat <- load_data("chickenpox", directory = tempdir())
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary

    model_autoregressive <- list(past_obs = rep(1, 7))
    fit <- glmstarma(chickenpox, model_autoregressive, W_hungary, family = vpoisson("log"),
                    covariates = list(population = population_hungary))
    x <- fitted(fit)
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 515)
    x <- fitted(fit, drop_init = FALSE)
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 522)


    
    mean_model <- list(past_obs = rep(1, 7))
    dispersion_model <- list(past_obs = 1)
    fit2 <- dglmstarma(chickenpox, mean_model, dispersion_model, mean_family = vquasipoisson("log"),
                   dispersion_link = "log",
                   wlist = W_hungary, 
                   mean_covariates = list(population = population_hungary))
    x <- fitted(fit2)
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 515)
    x <- fitted(fit2, drop_init = FALSE)
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 522)

    x <- fitted(fit2, return_value = "dispersion")
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 521)
    x <- fitted(fit2, return_value = "dispersion", drop_init = FALSE)
    expect_true(is.numeric(x))
    expect_equal(nrow(x), 20)
    expect_equal(ncol(x), 522)
})



testthat::skip_on_cran()
test_that("summary.glmstarma function works", {
    dat <- load_data("chickenpox", directory = tempdir())
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary

    fit <- glmstarma(chickenpox, list(past_obs = 1), W_hungary, family = vpoisson("log"),
                    covariates = list(population = population_hungary))
    x <- summary(fit)
    expect_s3_class(x, "summary.glmstarma")
    print(x)

    set.seed(123)
    fit <- glmstarma(chickenpox, list(past_obs = 1), W_hungary, family = vquasipoisson("log", dispersion = matrix(runif(length(chickenpox), min = 1, max = 10), nrow = nrow(chickenpox), ncol = ncol(chickenpox))),
                    covariates = list(population = population_hungary))
    x <- summary(fit)
    expect_s3_class(x, "summary.glmstarma")
    print(x)


    fit$coefficients <- NULL
    fit$coefficients_list <- NULL
    x <- summary(fit)
    expect_s3_class(x, "summary.glmstarma")
    expect_equal(length(x), 1)
    print(x)

    fit <- glmstarma(chickenpox, list(intercept = "inhomogeneous", past_obs = 1), W_hungary, family = vpoisson("log"),
                 covariates = list(population = population_hungary))
    x <- summary(fit)
    print(x)

    fit <- glmstarma(chickenpox, list(past_obs = cbind(c(1, 1), c(0, 1))), W_hungary, family = vpoisson("log"),
                 covariates = list(population = population_hungary))
    x <- summary(fit)
    print(x)

    fit <- glmstarma(chickenpox, list(past_obs = c(1, 1), past_mean = 1), W_hungary, family = vpoisson("log"),
                  covariates = list(population = population_hungary))
    x <- summary(fit)
    print(x)

    fit <- glmstarma(chickenpox, list(past_obs = c(1, 1)), W_hungary, family = vquasipoisson("identity"),
                  covariates = list(population = population_hungary))
    x <- summary(fit)
    print(x)
})




testthat::skip_on_cran()
test_that("summary.dglmstarma function works", {
    dat <- load_data("chickenpox", directory = tempdir())
    chickenpox <- dat$chickenpox
    population_hungary <- dat$population_hungary
    W_hungary <- dat$W_hungary

    fit2 <- dglmstarma(chickenpox, list(past_obs = 1, past_mean = 1), list(past_obs = 1), mean_family = vquasipoisson("identity"),
                   dispersion_link = "log",
                   wlist = W_hungary, 
                   mean_covariates = list(population = population_hungary))
    x <- summary(fit2)
    expect_s3_class(x, "summary.dglmstarma")
    print(x)

    fit2$mean$coefficients <- NULL
    fit2$mean$coefficients_list <- NULL
    fit2$dispersion$coefficients <- NULL
    fit2$dispersion$coefficients_list <- NULL
    x <- summary(fit2)
    expect_s3_class(x, "summary.dglmstarma")
    expect_equal(length(x), 1)


    fit2 <- dglmstarma(chickenpox, list(past_obs = 1, past_mean = 1), list(intercept = "inhomogeneous"), mean_family = vquasipoisson("identity"),
                   dispersion_link = "log",
                   wlist = W_hungary, 
                   mean_covariates = list(population = population_hungary),
                   control = list(use_fast_if_const_dispersion = TRUE))
    x <- summary(fit2)
    expect_s3_class(x, "summary.dglmstarma")
    print(x)

    x <- summary(fit2, alternative = "less")
    expect_s3_class(x, "summary.dglmstarma")
    print(x)

    x <- summary(fit2, alternative = "greater")
    expect_s3_class(x, "summary.dglmstarma")
    print(x)


    fit2 <- dglmstarma(chickenpox, list(past_obs = 1, intercept = "inhomogeneous"), list(past_obs = 1), mean_family = vquasipoisson("identity"),
                   dispersion_link = "log",
                   wlist = W_hungary)
    x <- summary(fit2)
    expect_s3_class(x, "summary.dglmstarma")
    print(x)

    fit2 <- dglmstarma(chickenpox, list(past_obs = 1), list(past_obs = 1), mean_family = vquasipoisson("identity"),
                   dispersion_link = "log",
                   wlist = W_hungary, 
                   dispersion_covariates = list(population = population_hungary))
    x <- summary(fit2)
    expect_s3_class(x, "summary.dglmstarma")
    print(x)

    fit2 <- dglmstarma(chickenpox, list(past_obs = rep(1, 3)), list(past_obs = 1, past_mean = 1), mean_family = vquasipoisson("log"),
                   dispersion_link = "identity",
                   wlist = W_hungary,
                   mean_covariates = list(population = population_hungary))
    x <- summary(fit2)
    expect_s3_class(x, "summary.dglmstarma")
    print(x)

})






## rest will be added later