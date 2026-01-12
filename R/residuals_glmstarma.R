# -----------------------------------------------------------------------------
# File: residuals_glmstarma.R
# Purpose: Implement S3 method 'residuals' for glmSTARMA models
# Author: Steffen Maletz
# Last modified: 2025-12-05
# -----------------------------------------------------------------------------

#' @rdname residuals-glmstarma
#' @name residuals.glmstarma
#' @aliases residuals.glmstarma residuals.dglmstarma
#' @title Residuals for glmstarma and dglmstarma Models
#' @description Compute residuals for fitted \code{glmstarma} and \code{dglmstarma} models.
#' @param object A fitted \code{glmstarma} or \code{dglmstarma} object.
#' @param type Type of residuals to compute. Options are \code{"response"} (raw residuals), \code{"pearson"}. See details.
#' @param drop_init Logical; if \code{TRUE}, initial first \code{max_time_lag} columns of residuals are dropped.
#' @param ignore_dispersion Logical; if \code{TRUE}, values are not scaled by the dispersion parameter
#' @return A matrix of residuals.
#' @details The \code{type} argument specifies the type of residuals to compute:
#' * \code{"response"}: Raw residuals, computed as the difference between observed and fitted values.
#' * \code{"pearson"}: Pearson residuals, defined as \deqn{r_i = \frac{y_i - \mu_i}{\sqrt{V(\mu_i)}}}, where \eqn{V(\mu_i)} is the variance function of the specified family.
#' * \code{"deviance"}: Deviance residuals, defined as \deqn{r_i = 2 \cdot (\ell(y_i; y_i) - \ell(y_i; \mu_i)),} i.e. the log-likelihood difference of a saturated model and the fitted model.
#' If \code{ignore_dispersion} is set to \code{FALSE}, pearson and deviance residuals are scaled by the dispersion parameter(s).
#' @examples
#' data("chickenpox")
#'
#' model_autoregressive <- list(past_obs = rep(1, 7))
#' fit <- glmstarma(chickenpox, model_autoregressive, W_hungary, family = vpoisson("log"),
#'                  covariates = list(population = population_hungary))
#' residuals(fit)
#'
#' mean_model <- list(past_obs = rep(1, 7))
#' dispersion_model <- list(past_obs = 1)
#' fit2 <- dglmstarma(chickenpox, mean_model, dispersion_model, mean_family = vquasipoisson("log"),
#'                    dispersion_link = "log",
#'                    wlist = W_hungary, 
#'                    mean_covariates = list(population = population_hungary))
#' residuals(fit2)
#' @exportS3Method stats::residuals
residuals.glmstarma <- function(object, type = c("response", "pearson", "deviance"), drop_init = TRUE, ignore_dispersion = TRUE){
    type <- match.arg(type)
    res <- object$ts - object$fitted.values
    if(type == "pearson"){
        res <- res / sqrt(object$family$variance(object$fitted.values, object$dispersion_est, ignore_dispersion))
    }
    if(type == "deviance"){
        res <- object$family$dev.resids(object$ts, object$fitted.values, object$dispersion_est, ignore_dispersion)
    }
    if(drop_init && object$max_time_lag > 0){
        res <- res[, -seq(object$max_time_lag), drop = FALSE]
    }
    return(res)
}


#' @rdname residuals-glmstarma
#' @exportS3Method stats::residuals
residuals.dglmstarma <- function(object, type = c("response", "pearson", "deviance"), drop_init = TRUE, ignore_dispersion = TRUE){
    type <- match.arg(type)
    res <- object$mean$ts - object$mean$fitted.values
    dispersion_est <- object$dispersion$dispersion_matrix
    if(ncol(dispersion_est) < ncol(object$mean$ts)){
        dispersion_est <- cbind(matrix(1, nrow(object$mean$ts), ncol(object$mean$ts) - ncol(dispersion_est)), dispersion_est)
    }
    if(type == "pearson"){
        res <- res / sqrt(object$mean$family$variance(object$mean$fitted.values, dispersion_est, ignore_dispersion))
    }
    if(type == "deviance"){
        res <- object$mean$family$dev.resids(object$mean$ts, object$mean$fitted.values, dispersion_est, ignore_dispersion)
    }
    if(drop_init && object$mean$max_time_lag > 0){
        res <- res[, -seq(object$mean$max_time_lag), drop = FALSE]
    }
    return(res)
}

