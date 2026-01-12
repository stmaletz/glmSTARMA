# -----------------------------------------------------------------------------
# File: vcov_glmstarma.R
# Purpose: Implement S3 method 'vcov' for glmSTARMA models
# Author: Steffen Maletz
# Last modified: 2025-12-14
# -----------------------------------------------------------------------------


#' @title Variance-Covariance Matrix for glmstarma and dglmstarma objects
#' @description Computes the variance-covariance matrix for \code{glmstarma} and \code{dglmstarma} objects.
#' @param object An object of class \code{glmstarma} or \code{dglmstarma} for which the variance-covariance matrix is to be computed.
#' @param return_value A character string specifying which variance-covariance matrix to return. Options are "mean", "dispersion", or "both". Default is "mean".
#' @return For \code{glmstarma} objects, the function returns the variance-covariance matrix of the mean model coefficients.
#' For \code{dglmstarma} objects, the function return depends on the \code{return_value} argument:
#' \item{mean}{Variance-covariance matrix for the mean model coefficients.}
#' \item{dispersion}{Variance-covariance matrix for the dispersion model coefficients.}
#' \item{both}{A list containing both variance-covariance matrices.}
#' @details
#' The variance-covariance matrix is computed using a sandwich estimator approach, which accounts for potential misspecification of the model.
#' The sandwich variance estimation is defined as \eqn{V = G^{-1} H G^{-1}}, where \eqn{G} is the expected information matrix and \eqn{H} is the empirical covariance matrix of the score functions.
#' In case of \code{dglmstarma} objects, separated variance-covariance matrices are computed for the mean and dispersion models because of the alternating estimation procedure.
#' @seealso \code{\link{vcov}}, \code{\link{glmstarma}}, \code{\link{dglmstarma}}
#' @examples
#' data("chickenpox")
#'
#' model_autoregressive <- list(past_obs = rep(1, 7))
#' fit <- glmstarma(chickenpox, model_autoregressive, W_hungary, family = vpoisson("log"),
#'                  covariates = list(population = population_hungary))
#' vcov(fit)
#'
#' mean_model <- list(past_obs = rep(1, 7))
#' dispersion_model <- list(past_obs = 1)
#' fit2 <- dglmstarma(chickenpox, mean_model, dispersion_model, mean_family = vquasipoisson("log"),
#'                    dispersion_link = "log",
#'                    wlist = W_hungary, 
#'                    mean_covariates = list(population = population_hungary))
#' vcov(fit2)
#' vcov(fit2, return_value = "dispersion")
#' vcov(fit2, return_value = "both")
#' @rdname sandwich_variance
#' @aliases sandwich_variance
#' @exportS3Method stats::vcov
vcov.dglmstarma <- function(object, return_value = c("mean", "dispersion", "both")){
    return_value <- match.arg(return_value)
    if(return_value %in% c("mean", "both")){
        mean_vcov <- object$mean$variance_estimation
    }
    if(return_value %in% c("dispersion", "both")){
        dispersion_vcov <- object$dispersion$variance_estimation
        if(is.null(dispersion_vcov)){
            dispersion_vcov <- 2 * diag(drop(object$dispersion$coefficients)^2) / ifelse(length(object$dispersion$coefficients) == 1, 
            object$mean$n_obs_effective - length(object$mean$coefficients), object$mean$n_obs_effective - length(object$mean$coefficients) / object$target_dim)
        }
    }
    if(return_value == "mean"){
        return(mean_vcov)
    } else if(return_value == "dispersion"){
        return(dispersion_vcov)
    } else {
        return(list(mean = mean_vcov, dispersion = dispersion_vcov))
    }
}


#' @rdname sandwich_variance
#' @exportS3Method stats::vcov
vcov.glmstarma <- function(object){
    return(object$variance_estimation)
}