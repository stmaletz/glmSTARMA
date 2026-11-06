# -----------------------------------------------------------------------------
# File: simulation_s3.R
# Purpose: S3 methods for the glmstarma_sim and dglmstarma_sim classes
# Author: Steffen Maletz
# Last modified: 2026-08-27
# -----------------------------------------------------------------------------

### Print methods

#' @exportS3Method base::print
print.glmstarma_sim <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  d <- dim(x$observations)
  cat("Simulated spatio-temporal GLM\n")
  cat(sprintf("Time-points: %d | Locations: %d\n", d[2], d[1]))
  if (length(x$parameters) > 0) {
    coeffs <- x$parameters
    
    print_coef_block(coeffs, x$model, digits, "True Coefficients")

  } else {
    cat("No coefficients\n")
  }
  invisible(x)
}


#' @exportS3Method base::print
print.dglmstarma_sim <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  d <- dim(x$observations)
  cat("Simulated spatio-temporal DGLM\n")
  cat(sprintf("Time-points: %d | Locations: %d\n", d[2], d[1]))
  if (length(x$parameters_mean) > 0) {
    coeffs <- x$parameters_mean
    print_coef_block(coeffs, x$mean_model, digits, "True Coefficients - Mean Model")
  } else {
    cat("No mean coefficients\n")
  }
  if (length(x$parameters_dispersion) > 0) {
    coeffs <- x$parameters_dispersion
    print_coef_block(coeffs, x$dispersion_model, digits, "True Coefficients - Dispersion Model")
  } else {
    cat("No dispersion coefficients\n")
  }
  invisible(x)
}





### Accessors


#' @title Accessors for Spatio-Temporal GLM and DGLM Simulation Objects
#' @name accessors_simulation
#' @aliases get_observations get_link_values get_model get_parameters
#' @description
#' These functions provide access to various components of the simulation objects returned by \code{\link{glmstarma.sim}} and \code{\link{dglmstarma.sim}}.
#'
#' @param object Object of class \code{"glmstarma_sim"} or \code{"dglmstarma_sim"}.
#' @param type (character) Specifies which model (\code{"mean"} or \code{"dispersion"}) the return value corresponds to. 
#'
#' @return The return value depends on the specific accessor function called:
#' - \code{get_observations}: Returns the simulated observations as a numeric matrix.
#' - \code{get_link_values}: Returns the linear predictor values as a numeric matrix.
#' - \code{get_parameters}: Returns the model parameters as a list.
#' - \code{get_model_orders}: Returns the model orders as a list.
#' - \code{get_pseudo_observations}: Returns the pseudo-observations for the dispersion model as a numeric matrix (only for \code{"dglmstarma_sim"} objects).
#' - \code{get_dispersion_values}: Returns the dispersion values as a numeric matrix.
#' - \code{get_family}: Returns the \code{stfamily} object used in the simulation.
#'
#' @examples
#' set.seed(42)
#' n_obs <- 200L
#' W <- generateW("rectangle", 100, 2, 10)
#' model_orders_mean <- list(intercept = "homogeneous", 
#'                           past_obs = 2, past_mean = 1, 
#'                           covariates = c(0, 0))
#' model_orders_dispersion <- list(intercept = "homogeneous", 
#'                                 past_obs = 1, 
#'                                 covariates = c(0, 0))
#'
#' covariates_mean <- list(season = SpatialConstant(sin(2 * pi / 12 * seq(n_obs))),
#'                   location = TimeConstant(rnorm(100, sd = 0.81)))
#'
#' covariates_dispersion <- list(season = SpatialConstant(sin(2 * pi / 24 * seq(n_obs))),
#'                         location = TimeConstant(runif(100)))
#'
#' params_mean <- list(intercept = 0.6, 
#'                     past_mean = matrix(c(0.2, 0.1), nrow = 2), 
#'                     past_obs = matrix(c(0.2, 0.1, 0.05), nrow = 3), 
#'                     covariates = matrix(c(0.9, 0.2), ncol = 2))
#' params_dispersion <- list(intercept = 0.5, 
#'                     past_obs = matrix(c(0.5, 0.2), nrow = 2), 
#'                     covariates = matrix(c(0.1, 0.75), ncol = 2))
#' family <- vnormal(copula = "frank", copula_param = 2)
#' sim <- dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
#'               model_orders_dispersion, mean_family = family, 
#'               wlist = W, pseudo_observations = "deviance", 
#'               mean_covariates = covariates_mean, 
#'               dispersion_covariates = covariates_dispersion)
#'
#'
#' get_observations(sim)
#' get_link_values(sim)
#' get_parameters(sim)
#' get_model_orders(sim)
#' get_pseudo_observations(sim)
#' get_dispersion_values(sim)
#' get_family(sim)
#'
#' @seealso [glmstarma.sim], [dglmstarma.sim]
#'
#' @export
get_observations <- function(object, ...) UseMethod("get_observations")

#' @exportS3Method
get_observations.glmstarma_sim <- function(object, ...) {
  object$observations
}

#' @exportS3Method
get_observations.dglmstarma_sim <- function(object, ...) {
  object$observations
}

#' @rdname accessors_simulation
#' @export
get_link_values <- function(object, ...) UseMethod("get_link_values")

#' @exportS3Method
get_link_values.glmstarma_sim <- function(object, ...) {
  object$link_values
}

#' @exportS3Method
get_link_values.dglmstarma_sim <- function(object, ...) {
  object$link_values
}



#' @rdname accessors_simulation
#' @export
get_parameters <- function(object, ...) UseMethod("get_parameters")

#' @exportS3Method 
get_parameters.glmstarma_sim <- function(object, ...) {
  object$parameters
}

#' @exportS3Method
get_parameters.dglmstarma_sim <- function(object, type = c("mean", "dispersion"), ...) {
    type <- match.arg(type)
    if (type == "mean") return(object$parameters_mean)
    if (type == "dispersion") return(object$parameters_dispersion)
}



#' @rdname accessors_simulation
#' @export
get_model_orders <- function(object, ...) UseMethod("get_model_orders")

#' @exportS3Method
get_model_orders.glmstarma_sim <- function(object, ...) {
  object$model
}

#' @exportS3Method
get_model_orders.dglmstarma_sim <- function(object, type = c("mean", "dispersion"), ...) {
    type <- match.arg(type)
    if (type == "mean") return(object$mean_model)
    if (type == "dispersion") return(object$dispersion_model)
}

#' @rdname accessors_simulation
#' @export
get_pseudo_observations <- function(object, ...) UseMethod("get_pseudo_observations")

#' @exportS3Method
get_pseudo_observations.dglmstarma_sim <- function(object, ...) {
  object$pseudo_observations
}


#' @rdname accessors_simulation
#' @export
get_dispersion_values <- function(object, ...) UseMethod("get_dispersion_values")

#' @exportS3Method
get_dispersion_values.glmstarma_sim <- function(object, ...) {
  if(length(object$family$dispersion) == 1){
    return(matrix(object$family$dispersion, nrow = nrow(object$observations), ncol = ncol(object$observations)))
  }
  if(length(object$family$dispersion) == nrow(object$observations)){
    return(matrix(object$family$dispersion, nrow = nrow(object$observations), ncol = ncol(object$observations)))
  }
  return(object$family$dispersion)
}

#' @exportS3Method
get_dispersion_values.dglmstarma_sim <- function(object, ...) {
  object$dispersion_values
}

#' @rdname accessors_simulation
#' @export
get_family <- function(object, ...) UseMethod("get_family")

#' @exportS3Method
get_family.glmstarma_sim <- function(object, ...) {
  object$family
}

#' @exportS3Method
get_family.dglmstarma_sim <- function(object, type = c("mean", "dispersion"), ...) {
  type <- match.arg(type)
  if (type == "mean") return(object$mean_family)
  if (type == "dispersion") return(object$dispersion_family)
}


### plot

#' @name plot.glmstarma_sim
#' @title Plot Simulation Objects
#' @description
#' Plot simulated time series returned by \code{\link{glmstarma.sim}} and \code{\link{dglmstarma.sim}}.
#' The methods support optional location filtering and aggregation across locations before plotting.
#'
#' @param x An object of class \code{glmstarma_sim} or \code{dglmstarma_sim}.
#' @param locations Optional numeric vector of row indices selecting the locations to include in the plot.
#' @param aggregate_fun Aggregation function applied across locations for each time point. The default is \code{mean}.
#' @param which_plot For \code{dglmstarma_sim} objects, selects whether to plot the observed values, pseudo-observations, or dispersion values.
#' @param ... Additional arguments passed to \code{aggregate_fun} and the underlying \code{plot()} call.
#'
#' @details
#' For \code{glmstarma_sim} objects, the method plots the aggregated observed series.
#' For \code{dglmstarma_sim} objects, \code{which_plot} controls whether the aggregated observed values,
#' pseudo-observations, or dispersion values are shown.
#'
#' @seealso \code{\link{glmstarma.sim}}, \code{\link{dglmstarma.sim}}
#' @exportS3Method base::plot
plot.glmstarma_sim <- function(x, locations = NULL, aggregate_fun = mean, ...) {
  observations <- get_observations(x)
  if(!is.null(locations)){
    stopifnot("locations must be a numeric vector" = is.numeric(locations),
          "locations must be finite" = all(is.finite(locations)),
          "locations must be row indices" = all(locations > 0 & locations == floor(locations)),
          "locations must not exceed the number of locations in the data" = all(locations <= nrow(observations)))
    observations <- observations[locations, , drop = FALSE]
  }
  aggregated_obs <- apply(observations, 2, aggregate_fun, ...)
  plot(aggregated_obs, type = "l", xlab = "Time", ylab = "Aggregated Observations", ...)
}

#' @rdname plot.glmstarma_sim
#' @exportS3Method base::plot
plot.dglmstarma_sim <- function(x, locations = NULL, which_plot = c("observations", "pseudo", "dispersion"), aggregate_fun = mean, ...) {
  which_plot <- match.arg(which_plot)
  observations <- switch(which_plot,
             observations = get_observations(x),
             pseudo = get_pseudo_observations(x),
             dispersion = get_dispersion_values(x))
  if(!is.null(locations)){
    stopifnot("locations must be a numeric vector" = is.numeric(locations),
          "locations must be finite" = all(is.finite(locations)),
          "locations must be row indices" = all(locations > 0 & locations == floor(locations)),
          "locations must not exceed the number of locations in the data" = all(locations <= nrow(observations)))
    observations <- observations[locations, , drop = FALSE]
  }
  aggregated_obs <- apply(observations, 2, aggregate_fun, ...)
  ylabel <- switch(which_plot,
           observations = "Aggregated Observations",
           pseudo = "Aggregated Pseudo-Observations",
           dispersion = "Aggregated Dispersion Values")
  plot(aggregated_obs, type = "l", xlab = "Time", ylab = ylabel, ...)
}




### summary

#' @exportS3Method base::summary
summary.glmstarma_sim <- function(object, ...) {
  obs <- object$observations
  obs_vector <- as.vector(obs)
  d <- dim(obs)
  nloc <- d[1]
  ntime <- d[2]

  family_label <- NA_character_
  if (!is.null(object$family)) {
    if(!is.null(object$family$distribution)) {
      family_label <- as.character(object$family$distribution)
    }
  }
  average_dispersion <- mean(object$family$dispersion)
  link <- object$family$link

  summary_list <- list(
    timepoints = ntime,
    locations = nloc,
    family = family_label,
    link = link,
    mean = mean(obs_vector),
    var = var(obs_vector),
    dispersion = average_dispersion,
    quantiles = stats::quantile(obs_vector, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  )

  class(summary_list) <- "summary.glmstarma_sim"
  return(summary_list)
}

#' @exportS3Method base::print
print.summary.glmstarma_sim <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Summary of simulated spatio-temporal GLM\n")
  cat(sprintf("Time-points: %d | Locations: %d\n", x$timepoints, x$locations))
  cat(sprintf("Family: %s\n", x$family))
  cat(sprintf("Link function: %s\n", x$link))
  cat(sprintf("Average dispersion: %.*f\n", digits, x$dispersion))

  cat(sprintf("Overall mean: %.*f | Var: %.*f\n", digits, x$overall$mean, digits, x$overall$var))
  cat("Quantiles (min, 25%, 50%, 75%, max):\n")
  print(format(x$overall$quantiles, digits = digits))

  invisible(x)
}


#' @exportS3Method base::summary
summary.dglmstarma_sim <- function(object, ...) {
  summarize_vector <- function(values) {
    values <- as.vector(values)
    list(
      mean = mean(values),
      var = var(values),
      quantiles = stats::quantile(values, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    )
  }

  family_label <- object$mean_family$distribution

  observations <- object$observations
  d <- dim(observations)
  nloc <- d[1]
  ntime <- d[2]

  summary_list <- list(
    timepoints = ntime,
    locations = nloc,
    mean_family = get_family_label(object$mean_family),
    mean_link = object$mean_family$link,
    dispersion_link =  object$dispersion_family$link,
    observations = summarize_vector(observations),
    pseudo_observations = summarize_vector(object$pseudo_observations),
    dispersion_values = summarize_vector(object$dispersion_values)
  )

  class(summary_list) <- "summary.dglmstarma_sim"
  return(summary_list)
}

#' @exportS3Method base::print
print.summary.dglmstarma_sim <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print_stats_block <- function(title, stats_object) {
    if (is.null(stats_object)) {
      return(invisible(NULL))
    }

    cat(sprintf("\n%s:\n", title))
    cat(sprintf("  Mean: %.*f | Var: %.*f\n", digits, stats_object$mean, digits, stats_object$var))
    cat("  Quantiles (min, 25%, 50%, 75%, max):\n")
    print(format(stats_object$quantiles, digits = digits))
    invisible(NULL)
  }

  cat("Summary of simulated spatio-temporal DGLM\n")
  cat(sprintf("Time-points: %d | Locations: %d\n", x$timepoints, x$locations))
  cat(sprintf("Mean family: %s\n", x$mean_family))
  cat(sprintf("Mean link: %s\n", x$mean_link))
  cat(sprintf("Dispersion link: %s\n", x$dispersion_link))

  print_stats_block("Observed values", x$observations)
  print_stats_block("Pseudo-observations", x$pseudo_observations)
  print_stats_block("Dispersion values", x$dispersion_values)

  invisible(x)
}
