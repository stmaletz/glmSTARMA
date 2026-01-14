# -----------------------------------------------------------------------------
# File: AIC_glmstarma.R
# Purpose: Functions for computing the Akaike Information Criterion (AIC)
# Author: Steffen Maletz
# Last modified: 2025-12-14
# -----------------------------------------------------------------------------



#' @name information_criteria
#' @title Information Criteria for glmstarma and dglmstarma objects
#' @description Compute AIC, BIC, and QIC and (Quasi-)log-likelihood for \code{glmstarma} and \code{dglmstarma} objects.
#' @param object An object of class \code{glmstarma} or \code{dglmstarma}
#' @param k Numeric; penalty per parameter to be used. Default is \code{2} (standard AIC).
#' @param adjust Logical; if \code{TRUE} (default), the (quasi-)log-likelihood is adjusted for the effective sample size. See Details.
#' @details
#' During model fitting, the (quasi-)log-likelihood is computed only on the last \code{n_eff} time-points, where \code{n_eff = n - max_time_lag_mean - max_time_lag_dispersion}.
#' Here \code{n} is the total number of time-points, \code{max_time_lag_mean} the maximum temporal lag in the mean model, and \code{max_time_lag_dispersion} the maximum temporal lag in the dispersion model (for \code{dglmstarma} objects).
#' If no dispersion model is present (class \code{glmstarma}), \code{max_time_lag_dispersion} is zero.
#'
#' To be more specific the (quasi-)log-likelihood calculated during model estimation is given by 
#' \deqn{\ell(\mathbf{\theta}) = \sum_{t = \tau}^n \sum_{i = 1}^p \ell_{i, t}(\mathbf{\theta}),}
#' where \eqn{\ell_{i, t}(\mathbf{\theta})} denotes the (quasi-)log-likelihood of the observation at location \eqn{i} at time \eqn{t},
#' and \eqn{\tau = n - n_{\mathrm{eff}}}. 
#'
#' This calculation of the (quasi-)log-likelihood introduces bias when comparing models of different temporal orders. 
#' If \code{adjust = TRUE}, the (quasi-)log-likelihood is rescaled to \code{n} observations by multiplying with \eqn{n / n_{\mathrm{eff}}}, before calculating the AIC, BIC, QIC or (quasi)-log-likelihood.
#' @return A numeric value for the (possibly adjusted) AIC, BIC, QIC or (quasi-)log-likelihood.
#' @seealso \code{\link{AIC}}, \code{\link{BIC}}, \code{\link{logLik}}, \code{\link{QIC}} 
#' @examples
#' dat <- load_data("chickenpox")
#' chickenpox <- dat$chickenpox
#' population_hungary <- dat$population_hungary
#' W_hungary <- dat$W_hungary
#'
#' model_autoregressive <- list(past_obs = rep(1, 7))
#' fit <- glmstarma(chickenpox, model_autoregressive, W_hungary, family = vpoisson("log"),
#'                  covariates = list(population = population_hungary))
#' AIC(fit)
#' BIC(fit)
#' logLik(fit)
#' QIC(fit)
#'
#' mean_model <- list(past_obs = rep(1, 7))
#' dispersion_model <- list(past_obs = 1)
#' fit2 <- dglmstarma(chickenpox, mean_model, dispersion_model, mean_family = vquasipoisson("log"),
#'                    dispersion_link = "log",
#'                    wlist = W_hungary, 
#'                    mean_covariates = list(population = population_hungary))
#' AIC(fit2)
#' BIC(fit2)
#' logLik(fit2)
#' QIC(fit2)
#' delete_glmSTARMA_data("chickenpox")  # Clean up cached data
NULL


#' @rdname information_criteria
#' @exportS3Method stats::AIC
AIC.glmstarma <- function(object, k = 2, adjust = TRUE){
    ll_weight <- ifelse(adjust, 1, object$n_obs_effective / ncol(object$ts))
    return(-2 * ll_weight * object$log_likelihood + k * length(object$coefficients))
}

#' @rdname information_criteria
#' @exportS3Method stats::AIC
AIC.dglmstarma <- function(object, k = 2, adjust = TRUE){
    if(adjust){
        return(object$aic)
    } else {
        n_obs_effective <- ncol(object$mean$ts) - (object$mean$max_time_lag + object$dispersion$max_time_lag)
        ll_weight <- n_obs_effective / ncol(object$mean$ts)
        return(object$aic + 2 * object$total_log_likelihood - 2 * ll_weight * object$total_log_likelihood)
    }
}