# -----------------------------------------------------------------------------
# File: print_covariates.R
# Purpose: Implement S3 method 'print' for printing covariate objects
# Author: Steffen Maletz
# Last modified: 2025-11-16
# -----------------------------------------------------------------------------


#' @exportS3Method base::print
print.spatial_constant <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
    cat("Spatially constant covariate - Values at time points:\n")
    print.default(x, digits = digits, ...)
}

#' @exportS3Method base::print
print.time_constant <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
    cat("Time constant covariate - Values at locations:\n")
    print.default(x, digits = digits, ...)
}
