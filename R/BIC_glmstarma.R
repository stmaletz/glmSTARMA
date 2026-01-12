# -----------------------------------------------------------------------------
# File: BIC_glmstarma.R
# Purpose: Functions for computing the Bayesian Information Criterion (BIC)
# Author: Steffen Maletz
# Last modified: 2025-12-14
# -----------------------------------------------------------------------------

# Functions are documented together with AIC, QIC and logLik in AIC_glmstarma.R

#' @rdname information_criteria
#' @exportS3Method stats::BIC
BIC.glmstarma <- function(object, adjust = TRUE){
    if(adjust){
        return(object$bic)
    } else {
        ll_weight <- ifelse(adjust, 1, object$n_obs_effective / ncol(object$ts))
        return(object$bic + 2 * object$log_likelihood - 2 * ll_weight * object$log_likelihood)
    }
}

#' @rdname information_criteria
#' @exportS3Method stats::BIC
BIC.dglmstarma <- function(object, adjust = TRUE) {
    if(adjust){
        return(object$bic)
    } else {
        n_obs_effective <- ncol(object$mean$ts) - (object$mean$max_time_lag + object$dispersion$max_time_lag)
        ll_weight <- n_obs_effective / ncol(object$mean$ts)
        return(object$bic + 2 * object$total_log_likelihood - 2 * ll_weight * object$total_log_likelihood)
    }
}
