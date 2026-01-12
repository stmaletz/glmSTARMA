# -----------------------------------------------------------------------------
# File: summary_glmstarma.R
# Purpose: Functions to summarize the results of a fitted (d)glmstarma model
# Author: Steffen Maletz
# Last modified: 2025-11-15
# -----------------------------------------------------------------------------


#' @rdname summary.glmstarma
#' @title Summarize the results of a glmstarma model
#'
#' @description This functions summarizes the model fit of a glmstarma model
#'
#' @param object An object of class \code{glmstarma}
#' @param ... Additional arguments passed to specific methods.
#' @return An object of class \code{summary.glmstarma} which contains the following elements
#'   * call: The function call to fit the model
#'   * coefficients: The estimated coefficients of the model with approximate standard errors, z- and p-values. See details.
#'   * distribution: The marginal distribution of the conditional observations.
#'   * link: The link-function used to connect the conditional mean with the linear process.
#'   * dispersion: The dispersion parameter of the conditional distribution
#'   * estimate_dispersion: logical value indicating whether \code{dispersion} was estimated (\code{TRUE}) or fixed by the distribution or user (\code{FALSE})
#'   * df: Number of estimated coefficients in the model
#'   * log_likelihood: The quasi-log-likelihood of the estimated model. See details.
#'   * aic: Akaike Information Criterion of the estimated model, see \link{information_criteria} with \code{adjust = TRUE}.
#'   * bic: Bayesian Information Criterion of the estimated model, see \link{information_criteria} with \code{adjust = TRUE}.
#'   * qic: Quasi Information Criterion of the estimated model, see \link{QIC} with \code{adjust = TRUE}.
#' @details  Standard errors, z-values and p-values are computed assuming asymptotic normality of the parameter estimation. The variance estimation is based on the sandwich estimator to adjust for quasi-maximum-likelihood estimation.
#'
#' If the model requires non-negative parameters, the p-values are adjusted accordingly. Note that this adjustment is only valid for testing single parameters against the null hypothesis of being zero.
#' If multiple parameters are tested simultaneously, or a linear combination of them, a different adjustment is necessary.
#' @seealso [glmstarma], [logLik], [AIC], [BIC], [QIC], [logLik.glmstarma], [AIC.glmstarma], [BIC.glmstarma]
#' @examples
#' data("chickenpox")
#'
#' model_autoregressive <- list(past_obs = rep(1, 7))
#' fit <- glmstarma(chickenpox, model_autoregressive, W_hungary, family = vpoisson("log"),
#'                  covariates = list(population = population_hungary))
#' summary(fit)
#' @exportS3Method base::summary
summary.glmstarma <- function(object, ...) {
  cl <- object$call
  coef_vector <- stats::coef(object)
  
  if (length(coef_vector) == 0) {
    return(structure(list(call = cl), class = "summary.glmstarma"))
  }
  
  # === Vorbereitung ===
  ll <- stats::logLik(object) 
  k <- length(coef_vector)
  n <- object$n_obs_effective
  model <- object$model
  coeffs <- object$coefficients_list
  fam <- object$family
  
  # === Intercept names ===
  param_names <- if (model$intercept == "homogeneous") {
    "(Intercept)"
  } else {
    paste0("(Intercept ", seq(object$target_dim), ")")
  }
  
  # === Helper: Create names ===
  get_param_names <- function(mat, rowlab, collab, prefix) {
    if (!is.null(mat) && nrow(mat) > 0) {
      idx <- which(mat == 1, arr.ind = TRUE)
      apply(idx, 1, function(i) {
        paste0(prefix, "_{", rowlab[i[1]], ", ", collab[i[2]], "}")
      })
    } else {
      character(0)
    }
  }
  
  # Past Means
  param_names <- c(param_names, get_param_names(model$past_mean, rownames(model$past_mean), colnames(model$past_mean), "past_mean"))
  
  # Past Observations
  param_names <- c(param_names, get_param_names(model$past_obs, rownames(model$past_obs), colnames(model$past_obs), "past_obs"))
  
  # Covariates (slightly different structure)
  if (!is.null(coeffs$covariates) && nrow(coeffs$covariates) > 0) {
    idx <- which(model$covariates == 1, arr.ind = TRUE)
    cov_names <- apply(idx, 1, function(i) {
      paste0(colnames(model$covariates)[i[2]], "_{", rownames(model$covariates)[i[1]], "}")
    })
    param_names <- c(param_names, cov_names)
  }

  # === results ===
  se <- sqrt(diag(object$variance_estimation) / n)
  if(!is.null(coeffs$past_mean) && (length(coeffs$past_mean) > 0)){
    k_mean <- sum(model$past_mean)
    k_intercept <- ifelse(model$intercept == "homogeneous", 1, object$target_dim)
    coef_vector <- c(coef_vector[(k_mean + 1):(k_mean + k_intercept)], coef_vector[1:k_mean], coef_vector[(k_mean + k_intercept + 1):length(coef_vector)])
    se <- c(se[(k_mean + 1):(k_mean + k_intercept)], se[1:k_mean], se[(k_mean + k_intercept + 1):length(se)])
  }


  z_values <- coef_vector / se
  p_vals <- 2 * stats::pnorm(abs(z_values), lower.tail = FALSE)
  if (fam$non_negative_parameters) {
    p_vals <- p_vals / 2
  }
  
  results <- data.frame(
    Estimate = coef_vector,
    `Std. Error` = se,
    `z value` = z_values,
    `Pr(>|z|)` = p_vals,
    row.names = param_names
  )
  names(results) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  # === Order of parameters ===
  n_intercept <- ifelse(model$intercept == "homogeneous", 1, object$target_dim)
  n_ma <- ifelse(is.null(model$past_mean), 0, sum(model$past_mean))

  if(is.null(object$dispersion)){
    dispersion_est <- fam$dispersion
  } else {
    dispersion_est <- object$dispersion
  }


  # === summary ===
  summary_list <- list(
    call = cl,
    coefficients = results,
    distribution = fam$distribution,
    link = fam$link,
    dispersion = dispersion_est,
    estimate_dispersion = fam$estimate_dispersion,
    df = length(coef_vector),
    log_likelihood = object$log_likelihood,
    aic = object$aic,
    bic = object$bic,
    qic = object$qic
  )
  
  class(summary_list) <- "summary.glmstarma"
  return(summary_list)
}


#' @rdname summary.dglmstarma
#' @title Summarize a dglmstarma Model
#'
#' @description This functions summarizes the model fit of a dglmstarma model.
#'
#' @param object An object of class \code{dglmstarma}
#' @param phi Numeric value indicating the null hypothesis value for the dispersion parameter test. Default is 1.
#' @param alternative Character string specifying the alternative hypothesis for the dispersion parameter test. Must be one of "two.sided" (default), "less" or "greater".
#' @param ... Additional arguments passed to specific methods.
#' @return An object of class \code{summary.dglmstarma} which contains the following elements
#'   * call: The function call to fit the model
#'   * coefficients_mean: The estimated coefficients of the mean model with approximate standard errors, z- and p-values. See details.
#'   * coefficients_dispersion: The estimated coefficients of the dispersion model with approximate standard errors, z- and p-values. See details.
#'   * distribution: The marginal distribution of the conditional observations.
#'   * link: The link-function used to connect the conditional mean with the linear process of the mean model.
#'   * dispersion_link: The link-function used to connect the dispersion with the linear process of the dispersion model.
#'   * dispersion_disp_parameter: The dispersion parameter of the dispersion family
#'   * df_mean: Number of estimated coefficients in the mean model
#'   * df_dispersion: Number of estimated coefficients in the dispersion model
#'   * log_likelihood: The quasi-log-likelihood of the estimated model. See details.
#'   * aic: Akaike Information Criterion of the estimated model, see \link{information_criteria} with \code{adjust = TRUE}.
#'   * bic: Bayesian Information Criterion of the estimated model, see \link{information_criteria} with \code{adjust = TRUE}.
#'   * qic: Quasi Information Criterion of the estimated model, see \link{QIC} with \code{adjust = TRUE}.
#' @details  Standard errors, z-values and p-values are computed assuming asymptotic normality of the parameter estimation. The variance estimation is based on the sandwich estimator to adjust for quasi-maximum-likelihood estimation.
#' If the model requires non-negative parameters, the p-values are adjusted accordingly. Note that this adjustment is only valid for testing single parameters against the null hypothesis of being zero.
#' If multiple parameters are tested simultaneously, or a linear combination of them, a different adjustment is necessary.
#' 
#' If the dispersion model is constant, i.e. it is only an intercept model, a test is performed to test whether the estimated dispersion parameter is significantly different from the null hypothesis value \code{phi}.
#' The alternative hypothesis can be specified via the \code{alternative} argument. This can be useful to test for overdispersion or underdispersion in the data.
#' @seealso [dglmstarma], [summary.glmstarma]
#' @examples
#' data("chickenpox")
#' mean_model <- list(past_obs = rep(1, 7))
#' dispersion_model <- list(past_obs = 1)
#' fit <- dglmstarma(chickenpox, mean_model, dispersion_model, mean_family = vquasipoisson("log"), 
#'                    dispersion_link = "log", wlist = W_hungary, 
#'                    mean_covariates = list(population = population_hungary))
#' summary(fit)
#' @exportS3Method base::summary
summary.dglmstarma <- function(object, phi = 1, alternative = c("two.sided", "less", "greater"), ...) {
  cl <- object$call
  coef_vector_mean <- object$mean$coefficients
  coef_vector_dispersion <- object$dispersion$coefficients
  
  
  if (length(coef_vector_mean) == 0 && length(coef_vector_dispersion) == 0) {
    return(structure(list(call = cl), class = "summary.dglmstarma"))
  }
  
  # === Vorbereitung ===
  ll <- object$log_likelihood_history[length(object$log_likelihood_history)]
  k_mean <- length(coef_vector_mean)
  k_dispersion <- length(coef_vector_dispersion)
  n_mean <- object$mean$n_obs_effective
  n_dispersion <- object$dispersion$n_obs_effective
  model_mean <- object$mean$model
  model_dispersion <- object$dispersion$model
  coeffs <- object$mean$coefficients_list
  coeffs_dispersion <- object$dispersion$coefficients_list
  fam <- object$mean$family
  fam_dispersion <- object$dispersion$family
  
  # === Intercept names ===
  param_names_mean <- if (model_mean$intercept == "homogeneous") {
    "(Intercept)"
  } else {
    paste0("(Intercept ", seq(object$target_dim), ")")
  }
  param_names_dispersion <- if (model_dispersion$intercept == "homogeneous") {
    "(Intercept)"
  } else {
    paste0("(Intercept ", seq(object$target_dim), ")")
  }

  # === Helper: Generate names ===
  get_param_names <- function(mat, rowlab, collab, prefix) {
    if (!is.null(mat) && nrow(mat) > 0) {
      idx <- which(mat == 1, arr.ind = TRUE)
      apply(idx, 1, function(i) {
        paste0(prefix, "_{", rowlab[i[1]], ", ", collab[i[2]], "}")
      })
    } else {
      character(0)
    }
  }
  
  # Past Means
  param_names_mean <- c(param_names_mean, get_param_names(model_mean$past_mean, rownames(model_mean$past_mean), colnames(model_mean$past_mean), "past_mean"))
  param_names_dispersion <- c(param_names_dispersion, get_param_names(model_dispersion$past_mean, rownames(model_dispersion$past_mean), colnames(model_dispersion$past_mean), "past_mean"))
  
  # Past Observations
  param_names_mean <- c(param_names_mean, get_param_names(model_mean$past_obs, rownames(model_mean$past_obs), colnames(model_mean$past_obs), "past_obs"))
  param_names_dispersion <- c(param_names_dispersion, get_param_names(model_dispersion$past_obs, rownames(model_dispersion$past_obs), colnames(model_dispersion$past_obs), "past_obs"))
  
  # covariates (slightly different structure)
  if (!is.null(coeffs$covariates) && nrow(coeffs$covariates) > 0) {
    idx <- which(model_mean$covariates == 1, arr.ind = TRUE)
    cov_names <- apply(idx, 1, function(i) {
      paste0(colnames(model_mean$covariates)[i[2]], "_{", rownames(model_mean$covariates)[i[1]], "}")
    })
    param_names_mean <- c(param_names_mean, cov_names)
  }
  if( !is.null(coeffs_dispersion$covariates) && nrow(coeffs_dispersion$covariates) > 0) {
    idx <- which(model_dispersion$covariates == 1, arr.ind = TRUE)
    cov_names <- apply(idx, 1, function(i) {
      paste0(colnames(model_dispersion$covariates)[i[2]], "_{", rownames(model_dispersion$covariates)[i[1]], "}")
    })
    param_names_dispersion <- c(param_names_dispersion, cov_names)
  }
  
  # === result table ===
  se_mean <- sqrt(diag(object$mean$variance_estimation) / n_mean)

  if(!is.null(coeffs$past_mean) && (length(coeffs$past_mean) > 0)){
    k_mean <- sum(model_mean$past_mean)
    k_intercept <- ifelse(model_mean$intercept == "homogeneous", 1, object$target_dim)
    coef_vector_mean <- c(coef_vector_mean[(k_mean + 1):(k_mean + k_intercept)], coef_vector_mean[1:k_mean], coef_vector_mean[(k_mean + k_intercept + 1):length(coef_vector_mean)])
    se_mean <- c(se_mean[(k_mean + 1):(k_mean + k_intercept)], se_mean[1:k_mean], se_mean[(k_mean + k_intercept + 1):length(se_mean)])
  }

  z_values_mean <- coef_vector_mean / se_mean
  p_vals_mean <- 2 * stats::pnorm(abs(z_values_mean), lower.tail = FALSE)
  if (fam$non_negative_parameters) {
    p_vals_mean <- p_vals_mean / 2
  }
  
  last_name <- "Pr(>|z|)"
  if(is.null(object$dispersion$variance_estimation) && object$control$use_fast_if_const_dispersion) {
    fam_dispersion$link <- "identity"
    fam_dispersion$dispersion <- NULL
    se_dispersion <- 2 * (drop(coef_vector_dispersion)^2) 
    nominator <- ifelse(length(coef_vector_dispersion) == 1, object$mean$n_obs_effective - length(coef_vector_mean), object$mean$n_obs_effective - length(coef_vector_mean) / object$target_dim)
    se_dispersion <- sqrt(se_dispersion / nominator)
    z_values_dispersion <- (coef_vector_dispersion - phi) / se_dispersion
    alternative <- match.arg(alternative)
    if (alternative == "two.sided") {
      p_vals_dispersion <- 2 * stats::pnorm(abs(z_values_dispersion), lower.tail = FALSE)
    } else if (alternative == "less") {
      p_vals_dispersion <- stats::pnorm(z_values_dispersion, lower.tail = TRUE)
      last_name <- "Pr(< z )"
    } else if (alternative == "greater") {
      p_vals_dispersion <- stats::pnorm(z_values_dispersion, lower.tail = FALSE)
      last_name <- "Pr(> z )"
    }
  } else {
    se_dispersion <- sqrt(diag(object$dispersion$variance_estimation) / n_dispersion)

    if(!is.null(coeffs_dispersion$past_mean) && (length(coeffs_dispersion$past_mean) > 0)){
      k_mean <- sum(model_dispersion$past_mean)
      k_intercept <- ifelse(model_dispersion$intercept == "homogeneous", 1, object$target_dim)
      coef_vector_dispersion <- c(coef_vector_dispersion[(k_mean + 1):(k_mean + k_intercept)], coef_vector_dispersion[1:k_mean], coef_vector_dispersion[(k_mean + k_intercept + 1):length(coef_vector_dispersion)])
      se_dispersion <- c(se_dispersion[(k_mean + 1):(k_mean + k_intercept)], se_dispersion[1:k_mean], se_dispersion[(k_mean + k_intercept + 1):length(se_dispersion)])
    }
    z_values_dispersion <- coef_vector_dispersion / se_dispersion
    p_vals_dispersion <- 2 * stats::pnorm(abs(z_values_dispersion), lower.tail = FALSE)
    if (fam_dispersion$non_negative_parameters) {
      p_vals_dispersion <- p_vals_dispersion / 2
    }
  }
  
  results_mean <- data.frame(
    Estimate = coef_vector_mean,
    `Std. Error` = se_mean,
    `z value` = z_values_mean,
    `Pr(>|z|)` = p_vals_mean,
    row.names = param_names_mean
  )
  
  results_dispersion <- data.frame(
    Estimate = coef_vector_dispersion,
    `Std. Error` = se_dispersion,
    `z value` = z_values_dispersion,
    `Pr(>|z|)` = p_vals_dispersion,
    row.names = param_names_dispersion
  )
  
  # === order of parameters ===
  n_intercept_mean <- ifelse(model_mean$intercept == "homogeneous", 1, object$target_dim)
  n_ma_mean <- ifelse(is.null(model_mean$past_mean), 0, sum(model_mean$past_mean))
  
  n_intercept_dispersion <- ifelse(model_dispersion$intercept == "homogeneous", 1, object$target_dim)
  n_ma_dispersion <- ifelse(is.null(model_dispersion$past_mean), 0, sum(model_dispersion$past_mean))
  
  names(results_mean) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  names(results_dispersion) <- c("Estimate", "Std. Error", "z value", last_name)
  
  # === summary ===
  summary_list <- list(
    call = cl,
    coefficients_mean = results_mean,
    coefficients_dispersion = results_dispersion,
    distribution = fam$distribution,
    link = fam$link,
    dispersion_link = fam_dispersion$link,
    dispersion_disp_parameter = fam_dispersion$dispersion,
    df_mean = length(coef_vector_mean),
    df_dispersion = length(coef_vector_dispersion),
    log_likelihood = object$total_log_likelihood,
    aic = object$aic,
    bic = object$bic,
    qic = object$qic
  )
  
  class(summary_list) <- "summary.dglmstarma"
  return(summary_list)
}  