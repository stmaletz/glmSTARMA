# -----------------------------------------------------------------------------
# File: predict_glmstarma.R
# Purpose: Implement S3 method 'predict' for glmSTARMA models
# Author: Steffen Maletz
# Last modified: 2025-11-16
# -----------------------------------------------------------------------------

#' @exportS3Method stats::predict
predict.glmstarma <- function(object, n.ahead = 1, type = c("response", "link", "sample"), newobs = NULL, newxreg = NULL){
    if(!is.null(object$family$copula)){
        if(object$family$copula %in% c("normal", "t")){
            copula_obj <- copula::ellipCopula(object$family$copula, param = object$family$copula_param, dim = object$target_dim)
        } else {
            copula_obj <- copula::archmCopula(object$family$copula, param = object$family$copula_param, dim = object$target_dim)
        }
    } else {
        copula_obj <- NULL
    }
    type <- match.arg(type)

    glmstarma_predict(n.ahead, type, object$ts, object$covariates, object$model, object$wlist_ar, object$wlist_ma, object$wlist_covariates, 
                      object$family, newxreg, newobs, object$fitted.values, object$control, object$coefficients_list, copula_obj)
}

#' @exportS3Method stats::predict
predict.dglmstarma <- function(object, n.ahead = 1, type = c("response", "link", "sample"), type_dispersion = c("response", "link"),
                                newobs = NULL, newxreg_mean = NULL, newxreg_dispersion = NULL) {
    type <- match.arg(type)
    type_dispersion <- match.arg(type_dispersion)
    if(!is.null(object$family$copula)){
        if(object$family$copula %in% c("normal", "t")){
            copula_obj <- copula::ellipCopula(object$family$copula, param = object$family$copula_param, dim = object$target_dim)
        } else {
            copula_obj <- copula::archmCopula(object$family$copula, param = object$family$copula_param, dim = object$target_dim)
        }
    } else {
        copula_obj <- NULL
    }

    dglmstarma_predict(n.ahead, type, type_dispersion, object$mean$ts, object$mean$model, 
                        object$dispersion$model, object$dispersion$pseudo_type, object$mean$family, object$dispersion$family,
                        object$mean$wlist_ar, object$mean$covariates, object$dispersion$covariates, object$dispersion$ts,
                        object$mean$wlist_past_mean, object$mean$wlist_covariates, object$dispersion$wlist_ar,
                        object$dispersion$wlist_past_mean, object$dispersion$wlist_covariates, object$mean$coefficients_list,
                        object$dispersion$coefficients_list, newxreg_mean, newxreg_dispersion, newobs, object$mean$fitted.values,
                        object$dispersion$fitted.values, object$control, copula_obj)
}

