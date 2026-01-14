# -----------------------------------------------------------------------------
# File: QIC_glmstarma.R
# Purpose: Implement QIC functions for glmSTARMA models
# Author: Steffen Maletz
# Last modified: 2025-12-14
# -----------------------------------------------------------------------------


#' @name QIC
#' @rdname QIC
#' @title Quasi Information Criterion (QIC) for glmstarma and dglmstarma objects
#'
#' @description Generic function to compute the QIC (Pan, 2001), a model selection criterion
#' commonly used for Generalized Estimating Equations (GEE) and related models.
#'
#' @param object Object of class \code{glmstarma} or \code{dglmstarma}.
#' @param adjust Logical; if \code{TRUE} (default), an adjustment for the temporal orders of the model is applied to the likelihood. See Details.
#' @param ... Additional arguments passed to specific methods.
#' @return A numeric value for the QIC.
#' @details
#' The quasi information criterion (QIC) has been proposed by Pan (2001) as alternative to Akaike's information criterion (AIC) which is properly adjusted for regression analysis based on the generalized estimating equations (GEE).
#' It is defined as
#' \deqn{QIC = -2 \cdot \ell + 2 \cdot \left(\mathrm{trace}(G_{\mu}^{-1} H_{\mu}) + \mathrm{trace}(G_{\phi}^{-1} H_{\phi})\right),}
#' where \eqn{\ell} is the (quasi-)log-likelihood of the estimated model, \eqn{G_{\mu}} is the expected information matrix of the regression parameters of the mean model, and \eqn{H_{\mu}} the empirical covariance matrix of the regression parameters of the mean model.
#' Similarly, \eqn{G_{\phi}} and \eqn{H_{\phi}} denote the corresponding matrices for the dispersion model (only for \code{dglmstarma} objects). For \code{glmstarma} objects, the second term reduces to \eqn{\mathrm{trace}(G_{\mu}^{-1} H_{\mu})}.
#' 
#' For more details on the calculation of \eqn{G} and \eqn{H}, see \link{sandwich_variance}
#'
#'
#' During model estimation, the (quasi-)log-likelihood is computed only on the last \code{n_eff} time-points, where \code{n_eff = n - max_time_lag_mean - max_time_lag_dispersion}. 
#' Here \code{n} is the total number of time-points, \code{max_time_lag_mean} the maximum temporal lag in the mean model, and \code{max_time_lag_dispersion} the maximum temporal lag in the dispersion model (for \code{dglmstarma} objects).
#' If no dispersion model is present (class \code{glmstarma}), \code{max_time_lag_dispersion} is zero.
#'
#' To be more specific the (quasi-)log-likelihood calculated during model estimation is given by 
#' \deqn{\ell(\mathbf{\theta}) = \sum_{t = \tau}^n \sum_{i = 1}^p \ell_{i, t}(\mathbf{\theta}),}
#' where \eqn{\ell_{i, t}(\mathbf{\theta})} denotes the (quasi-)log-likelihood of the observation at location \eqn{i} at time \eqn{t},
#' and \eqn{\tau = n - n_{\mathrm{eff}}}.  
#'
#' This calculation of the (quasi-)log-likelihood introduces bias when comparing models of different temporal orders. 
#' If \code{adjust = TRUE}, the (quasi-)log-likelihood is rescaled to \code{n} observations by multiplying with \eqn{n / n_{\mathrm{eff}}}, before calculating the QIC.
#'
#' @seealso \code{\link{AIC}}, \code{\link{BIC}}, \code{\link{logLik}}, \code{\link{information_criteria}} 
#' @references
#' Pan, W. (2001). Akaike's Information Criterion in Generalized Estimating Equations.
#' \emph{Biometrics}, 57(1), 120–125. \doi{10.1111/j.0006-341X.2001.00120.x}
#'
#' @examples
#' dat <- load_data("chickenpox")
#' chickenpox <- dat$chickenpox
#' population_hungary <- dat$population_hungary
#' W_hungary <- dat$W_hungary
#'
#' model_autoregressive <- list(past_obs = rep(1, 7))
#' fit <- glmstarma(chickenpox, model_autoregressive, W_hungary, family = vpoisson("log"),
#'                  covariates = list(population = population_hungary))
#' QIC(fit)
#'
#' mean_model <- list(past_obs = rep(1, 7))
#' dispersion_model <- list(past_obs = 1)
#' fit2 <- dglmstarma(chickenpox, mean_model, dispersion_model, mean_family = vquasipoisson("log"),
#'                    dispersion_link = "log",
#'                    wlist = W_hungary, 
#'                    mean_covariates = list(population = population_hungary))
#' QIC(fit2)
#' delete_glmSTARMA_data("chickenpox")  # Clean up cached data
#' @export
QIC <- function(object, ...) UseMethod("QIC")

#' @rdname QIC
#' @exportS3Method glmSTARMA::QIC
QIC.glmstarma <- function(object, adjust = TRUE, ...) {
    if (adjust) {
        return(object$qic)
    } else {
        ll_weight <- ifelse(adjust, 1, object$n_obs_effective / ncol(object$ts))
        return(object$qic + 2 * object$log_likelihood - 2 * ll_weight * object$log_likelihood)
    }
}

#' @rdname QIC
#' @exportS3Method glmSTARMA::QIC
QIC.dglmstarma <- function(object, adjust = TRUE, ...) {
    if(adjust){
        return(object$qic)
    } else {
        n_obs_effective <- ncol(object$mean$ts) - (object$mean$max_time_lag + object$dispersion$max_time_lag)
        ll_weight <- n_obs_effective / ncol(object$mean$ts)
        return(object$qic + 2 * object$total_log_likelihood - 2 * ll_weight * object$total_log_likelihood)
    }
}