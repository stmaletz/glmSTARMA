# -----------------------------------------------------------------------------
# File: fitted_glmstarma.R
# Purpose: Implement S3 method 'fitted' for glmSTARMA models
# Author: Steffen Maletz
# Last modified: 2025-12-05
# -----------------------------------------------------------------------------

#' @rdname fitted.glmstarma
#' @name fitted.glmstarma
#' @aliases fitted.dglmstarma
#' @title Fitted values for glmstarma Models
#' @description Compute fitted values for \code{glmstarma} and \code{dglmstarma} models.
#' @param object A \code{glmstarma} or \code{dglmstarma} object.
#' @param drop_init Logical; if \code{TRUE}, initial first \code{max_time_lag} columns of fitted values are dropped.
#' @param return_value Character; return fitted values of the mean model (\code{"mean"}) or the dispersion model (\code{"dispersion"}).
#' @return A matrix of fitted values.
#' @seealso \code{\link{fitted}}, \code{\link{glmstarma}}, \code{\link{dglmstarma}}
#' @examples
#' dat <- load_data("chickenpox")
#' chickenpox <- dat$chickenpox
#' population_hungary <- dat$population_hungary
#' W_hungary <- dat$W_hungary
#'
#' model_autoregressive <- list(past_obs = rep(1, 7))
#' fit <- glmstarma(chickenpox, model_autoregressive, W_hungary, family = vpoisson("log"),
#'                  covariates = list(population = population_hungary))
#' fitted.values(fit)
#'
#' mean_model <- list(past_obs = rep(1, 7))
#' dispersion_model <- list(past_obs = 1)
#' fit2 <- dglmstarma(chickenpox, mean_model, dispersion_model, mean_family = vquasipoisson("log"),
#'                    dispersion_link = "log",
#'                    wlist = W_hungary, 
#'                    mean_covariates = list(population = population_hungary))
#' fitted.values(fit2)
#' fitted.values(fit2, return_value = "dispersion")
#' delete_glmSTARMA_data("chickenpox")  # Clean up cached data
#' @exportS3Method stats::fitted
fitted.glmstarma <- function(object, drop_init = TRUE){
    if(drop_init && object$max_time_lag > 0){
        return(object$fitted.values[, -seq(object$max_time_lag), drop = FALSE])
    }
    return(object$fitted.values)
}

#' @rdname fitted.glmstarma
#' @exportS3Method stats::fitted
fitted.dglmstarma <- function(object, return_value = c("mean", "dispersion"), drop_init = TRUE){
    return_value <- match.arg(return_value)
    # Zugriff auf das passende fitted-Values-Objekt
    fitted_vals <- object[[return_value]]$fitted.values
    max_lag <- object[[return_value]]$max_time_lag

    # Optional: Initialwerte entfernen
    if (drop_init && max_lag > 0) {
        return(fitted_vals[, -seq(object$max_time_lag), drop = FALSE])
    }
    return(fitted_vals)
}
