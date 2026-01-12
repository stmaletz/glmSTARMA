# -----------------------------------------------------------------------------
# File: covariates.R
# Purpose: Implement covariate functions for glmSTARMA models
# Author: Steffen Maletz
# Last modified: 2025-11-12
# -----------------------------------------------------------------------------

#' @rdname TimeConstant
#' @title Creates a time constant covariate
#'
#' @description This functions assigns a \code{const} attribute set to \code{"time"} to a numeric vector.
#'
#' @param x (numeric vector) covariate values for each location. Values are used for each time-point.
#' @return The input vector with an additional attribute \code{const} set to \code{"time"}.
#' @details A time-constant covariate has the form \deqn{\mathbf{x}_t = \mathbf{x},} for some fixed vector \eqn{\mathbf{x} = (x_1, \ldots, x_p)',} i.e., it has the same value for all time-points \eqn{t}, where \eqn{x_i} is the value of the covariate at location \eqn{i}.
#' @seealso [SpatialConstant]
#' @export
TimeConstant <- function(x){
  stopifnot("x must be a vector" = is.vector(x),
            "x must be numeric" = is.numeric(x),
            "x must have a value for all locations" = length(x) >= 2L,
            "x cannot contain NA values" = !any(is.na(x)),
            "x cannot contain infinite values" = !any(is.infinite(x)))
  class(x) <- "time_constant"
  attr(x, "const") <- "time"
  return(x)
}

#' @rdname SpatialConstant
#' @title Creates a spatial constant covariate
#'
#' @description This function assigns a \code{const} attribute set to \code{"space"} to a numeric vector.
#'
#' @param x (numeric vector) covariate values for each time-point. Values are used for each location.
#' @return The input numeric vector with an additional attribute \code{const} set to \code{"space"}.
#' @details A spatial-constant covariate has the form \deqn{\mathbf{x}_t = x_t \mathbf{1}_p,} i.e., it has the same value \eqn{x_t} for all locations at time-point \eqn{t}.
#' @seealso [TimeConstant]
#' @export
SpatialConstant <- function(x){
  stopifnot("x must be a vector" = is.vector(x),
            "x must be numeric" = is.numeric(x),
            "x must have a value for all time-points" = length(x) >= 2L,
            "x cannot be the same value for all time-points" = length(unique(x)) > 1L,
            "x cannot contain NA values" = !any(is.na(x)),
            "x cannot contain infinite values" = !any(is.infinite(x)))
  class(x) <- "spatial_constant"
  attr(x, "const") <- "space"
  return(x)
}
