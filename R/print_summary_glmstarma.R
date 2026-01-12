# -----------------------------------------------------------------------------
# File: print_summary_glmstarma.R
# Purpose: Implement print functions for summary objects of glmSTARMA models
# Author: Steffen Maletz
# Last modified: 2025-11-14
# -----------------------------------------------------------------------------

#' @exportS3Method base::print
print.summary.glmstarma <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), ...){
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    cat("\nCoefficients:\n")
    stats::printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

    cat("\nMarginal Distribution: ", x$distribution, "\n", sep = "")
    cat("Link: ", x$link, "\n", sep = "")
    if(x$estimate_dispersion){
        cat("Dispersion Parameter estimated as ", x$dispersion, "\n", sep = "")
    } else {
        if(length(x$dispersion) > 1)
        {
            time_varying <- is.matrix(x$dispersion) && ncol(x$dispersion) > 1
            time_varying <- ifelse(time_varying, "space-time-varying", "space-varying")
            cat("Dispersion Parameter fixed at ", time_varying, " values\n", sep = "")
        } else {
            cat("Dispersion Parameter fixed at ", x$dispersion, "\n", sep = "")
        }
    }
    cat("\nNumber of coefficients: ", x$df, "\n", sep = "")
    cat("\nQuasi-Log-Likelihood: ", x$log_likelihood, "\n", sep = "")
    cat("AIC: ", x$aic, "\n", sep = "")
    cat("BIC: ", x$bic, "\n", sep = "")
    cat("QIC: ", x$qic, "\n", sep = "")
    cat("\n")
    invisible(x)
}

#' @exportS3Method base::print
print.summary.dglmstarma <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), ...){
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    
    cat_central <- function(title, line_width = 50) {
        padding <- floor((line_width - nchar(title)) / 2)
        cat("\n", strrep("=", line_width), "\n", sep = "")
        cat(strrep(" ", padding), title, "\n", sep = "")
        cat(strrep("=", line_width), "\n", sep = "")
    }
    
    cat_central("Coefficients of Mean Model")
    stats::printCoefmat(x$coefficients_mean, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

    cat("\nMarginal Distribution: ", x$distribution, "\n", sep = "")
    cat("Link: ", x$link, "\n", sep = "")

    cat_central("Coefficients of Dispersion Model")
    stats::printCoefmat(x$coefficients_dispersion, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

    cat("Link: ", x$dispersion_link, "\n", sep = "")
    if(!is.null(x$dispersion_disp_parameter)){
        cat("Dispersion Parameter of dispersion model fixed at ", x$dispersion_disp_parameter, "\n\n", sep = "")
    } else {
        cat("Dispersion Parameters estimated via methods of moments.\n\n", sep = "")
    }
    

    cat("Number of coefficients of mean model: ", x$df_mean, "\n", sep = "")
    cat("Number of coefficients of dispersion model: ", x$df_dispersion, "\n", sep = "")
    cat("\nQuasi-Log-Likelihood: ", x$log_likelihood, "\n", sep = "")
    cat("AIC: ", x$aic, "\n", sep = "")
    cat("BIC: ", x$bic, "\n", sep = "")
    cat("QIC: ", x$qic, "\n", sep = "")
    cat("\n")
    invisible(x)
}