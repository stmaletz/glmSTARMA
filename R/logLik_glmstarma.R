# -----------------------------------------------------------------------------
# File: logLik_glmstarma.R
# Purpose: S3 Method 'logLik' for glmstarma and dglmstarma objects
# Author: Steffen Maletz
# Last modified: 2025-12-14
# -----------------------------------------------------------------------------

# Functions are documented together with AIC, QIC and logLik in AIC_glmstarma.R

#' @rdname information_criteria
#' @exportS3Method stats::logLik
logLik.glmstarma <- function(object, adjust = TRUE){
    ll_weight <- ifelse(adjust, 1, object$n_obs_effective / ncol(object$ts))
    return(ll_weight * object$log_likelihood)
}

#' @rdname information_criteria
#' @exportS3Method stats::logLik
logLik.dglmstarma <- function(object, adjust = TRUE){
    if(adjust){
        return(object$total_log_likelihood)
    } else {
        n_obs_effective <- ncol(object$mean$ts) - (object$mean$max_time_lag + object$dispersion$max_time_lag)
        ll_weight <- n_obs_effective / ncol(object$mean$ts)
        return(ll_weight * object$total_log_likelihood)
    }
}
