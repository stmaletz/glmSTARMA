# -----------------------------------------------------------------------------
# File: coef_glmstarma.R
# Purpose: S3 Method 'coef' for glmstarma and dglmstarma objects
# Author: Steffen Maletz
# Last modified: 2025-12-05
# -----------------------------------------------------------------------------

#' @rdname coef.glmstarma
#' @aliases coef.dglmstarma
#' @title Extract Coefficients of glmstarma and dglmstarma Models
#' @description Extracts model coefficients from objects of class \code{glmstarma} and \code{dglmstarma}.
#' @param object An object of class \code{glmstarma} or \code{dglmstarma}
#' @param asList Logical; if \code{TRUE}, returns coefficients as a list, or otherwise as a numeric vector. Default is \code{FALSE}.
#' @return A numeric vector, or a list, of model coefficients.
#' @examples
#' data("chickenpox")
#'
#' model_autoregressive <- list(past_obs = rep(1, 7))
#' fit <- glmstarma(chickenpox, model_autoregressive, W_hungary, family = vpoisson("log"),
#'                  covariates = list(population = population_hungary))
#' coef(fit)
#' coef(fit, asList = TRUE)
#'
#' mean_model <- list(past_obs = rep(1, 7))
#' dispersion_model <- list(past_obs = 1)
#' fit2 <- dglmstarma(chickenpox, mean_model, dispersion_model, mean_family = vquasipoisson("log"),
#'                    dispersion_link = "log",
#'                    wlist = W_hungary, 
#'                    mean_covariates = list(population = population_hungary))
#' coef(fit2)
#' coef(fit2, asList = TRUE)
#' @exportS3Method stats::coef
coef.glmstarma <- function(object, asList = FALSE){
    if(asList){
        return(object$coefficients_list)
    } else {
        return(drop(object$coefficients))
    }
}

#' @rdname coef.glmstarma
#' @exportS3Method stats::coef
coef.dglmstarma <- function(object, asList = FALSE){
    if(asList){
        return(list(mean = object$mean$coefficients_list, 
                    dispersion = object$dispersion$coefficients_list))
    } else {
        return(c(drop(object$mean$coefficients), drop(object$dispersion$coefficients)))
    }
}
