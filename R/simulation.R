# -----------------------------------------------------------------------------
# File: simulation.R
# Purpose: Implement simulation functions for glmSTARMA models
# Author: Steffen Maletz
# Last modified: 2025-11-13
# -----------------------------------------------------------------------------

#' @rdname glmstarma.sim
#' @title Simulate spatial time-series based on generalized linear models
#' @description Generates a simulated multivariate time series based on a GLM-like model (see \code{\link{glmstarma}} for details)
#'
#' @param ntime Number of observation times to be simulated
#' @param parameters a named list specifying the parameters of the model to be simulated, which has the following elements:
#'   - `intercept` (numeric): Intercept parameter. If an inhomogeneous model is simulated, a value must be specified for each component of the time series.
#'   - `past_obs` (numeric matrix): Parameter values for the past observations. 
#'   - `past_mean` (numeric matrix): Parameter values for the past means.
#'   - `covariates` (numeric matrix): Parameter values for the covariates. 
#' @param model a named list specifying the model for the linear predictor, which can be of the following elements:
#'   - `intercept` (character): 'homogenous' (default) for a homogenous model, i.e. the same intercept for all components, 'inhomogenous' for inhomogenous models, i.e. an individual intercept for each component.
#'   - `past_obs` (integer vector/binary matrix): Maximal spatial orders for the time lags in `past_obs_time_lags`. A binary matrix can be passed as an alternative, with the entry in row \eqn{i} and column \eqn{j} indicating whether the \eqn{(i - 1)}-spatial lag for the \eqn{j}-th time lag is included in the model. If not specified, no regression on past observations is performed.
#'   - `past_obs_time_lags` (optional integer vector) indicates the time lags for regression on past observations. Defaults are `seq(length(past_obs))` (for vectors) and `seq(ncol(past_obs))` (for a matrix)
#'   - `past_mean` (integer vector/binary matrix): Spatial orders for the regression on the (latent) feedback process. Values can be entered in the same format as in `past_obs`. If not specified, no regression to the feedback process is performed.
#'   - `past_mean_time_lags` (optional integer vector) indicates the time lags for regression on past values of the feedback process. Defaults are `seq(length(past_mean))` (for vectors) and `seq(ncol(past_mean))` (for a matrix)
#'   - `covariates` (integer vector/binary matrix) spatial orders for the covariate processes passed in the argument `covariates`. The values can be passed as in `past_obs` and `past_means`, where the jth entry or column represents the jth covariable. If no values are specified but covariates are included, the spatial order 0 is used by default, which corresponds to the first matrix in argument `wlist_covariates`.
#' @param wlist A list of quadratic matrices, with the same dimension as the time series, which describe the spatial dependencies. Row-normalized matrices are recommended.
#' @param covariates List of covariates, containing matrices or returns of the covariate functions of this package (see also \code{\link{TimeConstant}}, \code{\link{SpatialConstant}}). The matrices must have the same dimension as \code{ts}
#' @param wlist_past_mean (Optional) List of matrices, which describes spatial dependencies for the past mean. If this is \code{NULL}, the matrices from \code{wlist} are used.
#' @param wlist_covariates (Optional) List of matrices, which describes spatial dependencies for the covariates. If this is \code{NULL}, the matrices from \code{wlist} are used.
#' @param family An object of class \code{stfamily} that specifies the marginal distributions and the type of model fitted.
#' @param n_start Number of observations to be used for the burn-in period
#' @param control A list of parameters for controlling the fitting process. This list is passed to \code{\link{glmstarma.control}}.
#' @return a named list with the following elements:
#'   - `observations` (numeric matrix): The simulated time series
#'   - `link_values` (numeric matrix): The underlying linear predictor resulting from the model and simulation
#'   - `model` (list): The model used for the simulation
#'   - `parameters` (list): The true parameters used for the simulation
#' @examples
#' set.seed(42)
#' n_obs <- 200L
#' W <- generateW("rectangle", 100, 2, 10)
#' model_orders <- list(intercept = "homogeneous", past_obs = 2, past_mean = 1)
#' parameter <- list(intercept = 0.5, past_obs = matrix(c(0.3, 0.2, 0.05), nrow = 3), 
#'                   past_mean = matrix(c(0.1, 0.05), nrow = 2),
#'                   covariates = c(0.75, 0.5))
#' covariates <- list(season = SpatialConstant(sin(2* pi / 12 * seq(n_obs))),
#'                    location = TimeConstant(rnorm(100, sd = 0.81)))
#' # Simulation using marginal poisson distribution
#' glmstarma.sim(n_obs, parameter, model_orders, W, covariates, family = vpoisson("log"))
#' # Simulation using negative binomial marginals
#' glmstarma.sim(n_obs, parameter, model_orders, W, covariates, 
#'               family = vnegative.binomial(dispersion = 3))
#' @export
glmstarma.sim <- function(ntime, parameters, model, family = NULL, wlist, covariates = list(), wlist_past_mean = NULL, wlist_covariates = NULL, n_start = 100L, control = list()){
  stopifnot("ntime must be a numeric value" = is.numeric(ntime),
            "ntime must be a finite numeric value" = is.finite(ntime),
            "ntime must be a positive integer value" = (ntime > 0 & ntime == floor(ntime)),
            "family must be specified" = !is.null(family), 
            "family must be of class stfamily" = inherits(family, "stfamily"),
            "covariates must be submitted in a list" = is.list(covariates),
            "wlist must be a list of numeric matrices" = is.list(wlist),
            "wlist_covariates must be a list of matrices" = is.null(wlist_covariates) | is.list(wlist_covariates),
            "wlist_past_mean must be a list of matrices" = is.null(wlist_past_mean) | is.list(wlist_past_mean),
            "n_start must be numeric" = is.numeric(n_start),
            "n_start must be a finite numeric value" = is.finite(n_start),
            "n_start must be a positive integer value" = (n_start > 0 & n_start == floor(n_start)),
            "parameters must be submitted in a list" = is.list(parameters),
            "control must be a list" = is.list(control),
            "wlist must not be an empty list" = length(wlist) > 0)

    dim <- c(wlist_check(wlist),
        wlist_check(wlist_past_mean),
        wlist_check(wlist_covariates))
    dim <- dim[!is.na(dim)]
    dim <- unique(dim)
    if(length(dim) != 1){
        stop("All wlist matrices must have the same dimension.")
    }

    if(is.null(model$covariates) && length(covariates) > 0){
        model$covariates <- rep(0, length(covariates))
    }

    temp <- model_and_parameter_check(model, parameters, dim)
    model <- temp$model
    parameters <- temp$parameters
    stopifnot("n_start must be larger than the maximum time lag in the model" = n_start >= max(c(0, model$past_obs_time_lags, model$past_mean_time_lags)))

    if(ifelse(is.null(parameters$past_obs), 0, sum(abs(parameters$past_obs))) + ifelse(is.null(parameters$past_mean), 0, sum(abs(parameters$past_mean))) >= 1){
        warning("Parameters of the model may not ensure stability of the process. Use with caution.")
    }
    if(family$non_negative_parameters && any(unlist(parameters) < 0)){
        stop("All parameters in the mean model must be non-negative for the selected mean family.")
    }

    wlist_length <- length(wlist)
    check_length_of_wlist(nrow(model$past_obs), wlist_length, 0, "wlist")
    check_length_of_wlist(nrow(model$past_mean), length(wlist_past_mean), wlist_length, "wlist_past_mean")
    if(length(covariates) > 0)
    {
        stopifnot("Model orders for covariates do not match the number of covariates" = ncol(model$covariates) == length(covariates),
                      "covariates must be submitted in a list" = covariate_check(covariates, ntime, dim, family))
        check_length_of_wlist(nrow(model$covariates), length(wlist_covariates), wlist_length, "wlist_covariates")
        if(is.null(names(covariates))){
            colnames(model$covariates) <- paste0("X", seq(length(covariates)))
            colnames(parameters$covariates) <- paste0("X", seq(length(covariates)))
            names(covariates) <- paste0("X", seq(length(covariates)))
        } else {
            names_temp <- names(covariates)
            names_temp[names_temp == ""] <- paste0("X", which(names_temp == ""))
            colnames(model$covariates) <- names_temp
            colnames(parameters$covariates) <- names_temp
            names(covariates) <- names_temp
        }

    }

    control <- do.call("glmstarma_sim.control", control)
    control$parameter_init <- parameters

    if(!is.null(family$copula)){
        if(family$copula %in% c("normal", "t")){
            copula_obj <- copula::ellipCopula(family$copula, param = family$copula_param, dim = nrow(wlist[[1]]))
        } else {
            copula_obj <- copula::archmCopula(family$copula, param = family$copula_param, dim = nrow(wlist[[1]]))
        }
    } else {
        copula_obj <- NULL
    }
    result <- glmstarma_sim_cpp(ntime, parameters, model, wlist, wlist_past_mean, family, covariates, wlist_covariates, copula_obj, n_start, control)
    return(result)
}


#' @rdname dglmstarma.sim
#' @title Simulate spatial time-series based on double generalized linear models
#' @description Generates a simulated multivariate time series based on a GLM-like model (see \code{\link{dglmstarma}} for details)
#'
#' @param ntime Number of observation times to be simulated
#' @param parameters_mean a named list specifying the parameters of the model to be simulated:
#'   - `intercept` (numeric): Intercept parameter. If an inhomogeneous model is simulated, a value must be specified for each component of the time series.
#'   - `past_obs` (numeric matrix): Parameter values for the past observations. 
#'   - `past_mean` (numeric matrix): Parameter values for the past means.
#'   - `covariates` (numeric matrix): Parameter values for the covariates.
#' @param parameters_dispersion a named list specifying the parameters of the dispersion model to be simulated, with the same possible elements as in \code{parameters_mean}.
#' @param model_mean a named list specifying the model for the linear predictor, which can be of the following elements:
#'   - `intercept` (character): 'homogenous' (default) for a homogenous model, i.e. the same intercept for all components, 'inhomogenous' for inhomogenous models, i.e. an individual intercept for each component.
#'   - `past_obs` (integer vector/binary matrix): Maximal spatial orders for the time lags in `past_obs_time_lags`. A binary matrix can be passed as an alternative, with the entry in row \eqn{i} and column \eqn{j} indicating whether the \eqn{(i - 1)}-spatial lag for the \eqn{j}-th time lag is included in the model. If not specified, no regression on past observations is performed.
#'   - `past_obs_time_lags` (optional integer vector) indicates the time lags for regression on past observations. Defaults are `seq(length(past_obs))` (for vectors) and `seq(ncol(past_obs))` (for a matrix)
#'   - `past_mean` (integer vector/binary matrix): Spatial orders for the regression on the (latent) feedback process. Values can be entered in the same format as in `past_obs`. If not specified, no regression to the feedback process is performed.
#'   - `past_mean_time_lags` (optional integer vector) indicates the time lags for regression on past values of the feedback process. Defaults are `seq(length(past_mean))` (for vectors) and `seq(ncol(past_mean))` (for a matrix)
#'   - `covariates` (integer vector/binary matrix) spatial orders for the covariate processes passed in the argument `covariates`. The values can be passed as in `past_obs` and `past_means`, where the \eqn{j}-th entry or column represents the \eqn{j}-th covariable. If no values are specified but covariates are included, the spatial order 0 is used by default, which corresponds to the first matrix in argument `wlist_covariates`.
#' @param model_dispersion a named list specifying the model for the dispersion linear predictor, with the same possible elements as in \code{model_mean}. Orders supplied in \code{past_obs} are applied to the pseudo-observations.
#' @param mean_family An object of class \code{stfamily} that specifies the marginal distributions of the observations and the link-function for the mean model.
#' @param dispersion_link Link function for the dispersion model. Possible values are \code{"log"} (default), \code{"identity"}, and \code{"inverse"}.
#' @param wlist A list of quadratic matrices, with the same dimension as the time series, which describe the spatial dependencies. Row-normalized matrices are recommended.
#' @param mean_covariates List of covariates included in the mean model, containing matrices or returns of the covariate functions of this package (see also \code{\link{TimeConstant}}, \code{\link{SpatialConstant}}).
#' @param dispersion_covariates List of covariates included in the dispersion model.
#' @param pseudo_observations Method to generate the  pseudo-observations for the dispersion model. Possible values are \code{"deviance"} (default) and \code{"pearson"}.
#' @param wlist_past_mean (Optional) List of matrices, which describes spatial dependencies for the past mean. If this is \code{NULL}, the matrices from \code{wlist} are used. 
#' @param wlist_covariates (Optional) List of matrices, which describes spatial dependencies for the covariates. If this is \code{NULL}, the matrices from \code{wlist} are used.
#' @param wlist_pseudo_obs (Optional) List of matrices, which describes spatial dependencies for the pseudo-observations in the dispersion model. If this is \code{NULL}, the matrices from \code{wlist} are used.
#' @param wlist_past_dispersion (Optional) List of matrices, which describes spatial dependencies for the past dispersion in the dispersion model. If this is \code{NULL}, the matrices from \code{wlist} are used.
#' @param wlist_covariates_dispersion (Optional) List of matrices, which describes spatial dependencies for the covariates in the dispersion model. If this is \code{NULL}, the matrices from \code{wlist} are used.
#' @param n_start Number of observations to be used for the burn-in period
#' @param control A list of parameters for controlling the fitting process. This list is passed to \code{\link{dglmstarma.control}}.
#' @return a named list with the following elements:
#'   - `observations` (numeric matrix): The simulated time series
#'   - `link_values` (numeric matrix): The underlying linear predictor resulting from the model and simulation
#'   - `pseudo_observations` (numeric matrix): The pseudo-observations generated for the dispersion model
#'   - `dispersion_values` (numeric matrix): The dispersion values resulting from the dispersion model
#'   - `mean_model` (list): The mean model used for the simulation
#'   - `dispersion_model` (list): The dispersion model used for the simulation
#'   - `parameters_mean` (list): The true parameters used for the mean model
#'   - `parameters_dispersion` (list): The true parameters used for the dispersion model
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
#' dglmstarma.sim(n_obs, params_mean, params_dispersion, model_orders_mean, 
#'               model_orders_dispersion, mean_family = family, 
#'               wlist = W, pseudo_observations = "deviance", 
#'               mean_covariates = covariates_mean, 
#'               dispersion_covariates = covariates_dispersion)
#' @export
dglmstarma.sim <- function(ntime, parameters_mean, parameters_dispersion, model_mean, model_dispersion, 
                           mean_family = NULL, dispersion_link = c("log", "identity", "inverse"),
                           wlist = NULL, mean_covariates = list(), dispersion_covariates = list(), pseudo_observations = c("deviance", "pearson"),
                           wlist_past_mean = NULL, wlist_covariates = NULL, wlist_pseudo_obs = NULL, 
                           wlist_past_dispersion = NULL, wlist_covariates_dispersion = NULL, n_start = 100L, control = list()){
  stopifnot("ntime must be a numeric value" = is.numeric(ntime),
            "ntime must be a finite numeric value" = is.finite(ntime),
            "ntime must be a positive integer value" = (ntime > 0 & ntime == floor(ntime)),
            "mean_family must be specified" = !is.null(mean_family),
            "mean_family must be of class stfamily" = inherits(mean_family, "stfamily"),
            "The vpoisson family is not supported for dglmstarma models. Please use vquasipoisson instead." = mean_family$distribution != "poisson",
            "The vbinomial family is not supported for dglmstarma models. Please use vquasibinomial instead." = mean_family$distribution != "binomial",
            "mean_covariates must be submitted in a list" = is.list(mean_covariates),
            "dispersion_covariates must be submitted in a list" = is.list(dispersion_covariates),
            "wlist must be a list of numeric matrices" = is.list(wlist),
            "wlist_covariates must be a list of matrices" = is.null(wlist_covariates) | is.list(wlist_covariates),
            "n_start must be numeric" = is.numeric(n_start),
            "n_start must be a finite numeric value" = is.finite(n_start),
            "n_start must be a positive integer value" = (n_start > 0 & n_start == floor(n_start)),
            "parameters must be submitted in a list" = is.list(parameters_mean) & is.list(parameters_dispersion),
            "control must be a list" = is.list(control),
            "wlist must not be an empty list" = length(wlist) > 0)

        
    ## Check all wlist-Arguments
    dim <- c(wlist_check(wlist),
        wlist_check(wlist_past_mean),
        wlist_check(wlist_covariates),
        wlist_check(wlist_pseudo_obs),
        wlist_check(wlist_past_dispersion),
        wlist_check(wlist_covariates_dispersion))
    dim <- dim[!is.na(dim)]
    dim <- unique(dim)
    if(length(dim) != 1){
        stop("All wlist matrices must have the same dimension.")
    }

    dispersion_link <- match.arg(dispersion_link)
    dispersion_family <- vgamma(dispersion_link, dispersion = 2)

    pseudo_observations <- match.arg(pseudo_observations)
    if(mean_family$distribution == "negative_binomial" && pseudo_observations == "deviance"){
        pseudo_observations <- "pearson"
    }

    # default covariate order if not specified
    if(is.null(model_mean$covariates) && length(mean_covariates) > 0){
        model_mean$covariates <- rep(0, length(mean_covariates))
    }
    if(is.null(model_dispersion$covariates) && length(dispersion_covariates) > 0){
        model_dispersion$covariates <- rep(0, length(dispersion_covariates))
    }

    temp <- model_and_parameter_check(model_mean, parameters_mean, dim)
    model_mean <- temp$model
    parameters_mean <- temp$parameters
    if(ifelse(is.null(parameters_mean$past_obs), 0, sum(abs(parameters_mean$past_obs))) + ifelse(is.null(parameters_mean$past_mean), 0, sum(abs(parameters_mean$past_mean))) >= 1){
        warning("Parameters of the mean model may not ensure stability of the process. Use with caution.")
    }
    if(mean_family$non_negative_parameters && any(unlist(parameters_mean) < 0)){
        stop("All parameters in the mean model must be non-negative for the selected mean family.")
    }
    param_unlist <- unlist(parameters_mean)
    stopifnot("All parameters must be finite values" = all(is.finite(param_unlist)),
              "All parameters must be numeric values" = is.numeric(param_unlist))

    temp <- model_and_parameter_check(model_dispersion, parameters_dispersion, dim)
    model_dispersion <- temp$model
    parameters_dispersion <- temp$parameters
    if(ifelse(is.null(parameters_dispersion$past_obs), 0, sum(abs(parameters_dispersion$past_obs))) + ifelse(is.null(parameters_dispersion$past_mean), 0, sum(abs(parameters_dispersion$past_mean))) >= 1){
        warning("Parameters of the dispersion model may not ensure stability of the process. Use with caution.")
    }
    if(dispersion_family$non_negative_parameters && any(unlist(parameters_dispersion) < 0)){
        stop("All parameters in the dispersion model must be non-negative for the selected dispersion link.")
    }
    param_unlist <- unlist(parameters_dispersion)
    stopifnot("All parameters must be finite values" = all(is.finite(param_unlist)),
              "All parameters must be numeric values" = is.numeric(param_unlist))

    ## Check if enough matrices in wlist-Arguments
    wlist_length <- length(wlist)
    check_length_of_wlist(nrow(model_mean$past_obs), wlist_length, 0, "wlist")
    check_length_of_wlist(nrow(model_mean$past_mean), length(wlist_past_mean), wlist_length, "wlist_past_mean")
    if(length(mean_covariates) > 0)
    {
        stopifnot("Model orders for mean_covariates do not match the number of mean_covariates" = ncol(model_mean$covariates) == length(mean_covariates),
                  "mean_covariates must be submitted in a list" = covariate_check(mean_covariates, ntime, dim, mean_family))
        check_length_of_wlist(nrow(model_mean$covariates), length(wlist_covariates), wlist_length, "wlist_covariates")
        if(is.null(names(mean_covariates))){
            colnames(model_mean$covariates) <- paste0("X", seq(length(mean_covariates)))
            colnames(parameters_mean$covariates) <- paste0("X", seq(length(mean_covariates)))
            names(mean_covariates) <- paste0("X", seq(length(mean_covariates)))
        } else {
            names_temp <- names(mean_covariates)
            names_temp[names_temp == ""] <- paste0("X", which(names_temp == ""))
            colnames(model_mean$covariates) <- names_temp
            colnames(parameters_mean$covariates) <- names_temp
            names(mean_covariates) <- names_temp
        }
    }
    check_length_of_wlist(nrow(model_dispersion$past_obs), length(wlist_pseudo_obs), wlist_length, "wlist_pseudo_obs")
    check_length_of_wlist(nrow(model_dispersion$past_mean), length(wlist_past_dispersion), wlist_length, "wlist_past_dispersion")
    if(length(dispersion_covariates) > 0)
    {
        stopifnot("Model orders for dispersion_covariates do not match the number of dispersion_covariates" = ncol(model_dispersion$covariates) == length(dispersion_covariates),
                  "dispersion_covariates must be submitted in a list" = covariate_check(dispersion_covariates, ntime, dim, dispersion_family))
        check_length_of_wlist(nrow(model_dispersion$covariates), length(wlist_covariates_dispersion), wlist_length, "wlist_covariates_dispersion")
        if(is.null(names(dispersion_covariates))){
            colnames(model_dispersion$covariates) <- paste0("X", seq(length(dispersion_covariates)))
            colnames(parameters_dispersion$covariates) <- paste0("X", seq(length(dispersion_covariates)))
            names(dispersion_covariates) <- paste0("X", seq(length(dispersion_covariates)))
        } else {
            names_temp <- names(dispersion_covariates)
            names_temp[names_temp == ""] <- paste0("X", which(names_temp == ""))
            colnames(model_dispersion$covariates) <- names_temp
            colnames(parameters_dispersion$covariates) <- names_temp
            names(dispersion_covariates) <- names_temp
        }
    }
    
    control <- do.call("glmstarma_sim.control", control)
    control$parameter_init <- parameters_mean

    if(!is.null(mean_family$copula)){
        if(mean_family$copula %in% c("normal", "t")){
            copula_obj <- copula::ellipCopula(mean_family$copula, param = mean_family$copula_param, dim = nrow(wlist[[1]]))
        } else {
            copula_obj <- copula::archmCopula(mean_family$copula, param = mean_family$copula_param, dim = nrow(wlist[[1]]))
        }
    } else {
        copula_obj <- NULL
    }

    result <- dglmstarma_sim_cpp(ntime, parameters_mean, parameters_dispersion,
                                 mean_family, dispersion_family,
                                 model_mean, model_dispersion, pseudo_observations, wlist, 
                                 mean_covariates, dispersion_covariates, wlist_past_mean, wlist_covariates,
                                 wlist_pseudo_obs, wlist_past_dispersion, wlist_covariates_dispersion,
                                 control, copula_obj, n_start)
    return(result)
}