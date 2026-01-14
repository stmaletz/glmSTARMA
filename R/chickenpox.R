#' @title Chickenpox Infections in Hungary
#'
#' @description Multivariate count time series consisting of weekly chickenpox infections in the districts of Hungary.
#'
#' @usage
#' load_data("chickenpox")
#'
#' @format
#' \describe{
#'   \item{chickenpox}{
#'     A matrix with counts of chickenpox infections (rows = districts, columns = time points).
#'   }
#'   \item{W_hungary}{
#'     A list of matrices containing spatial weight matrices:
#'     \enumerate{
#'      \item Identity matrix.  
#'      \item Row-normalized adjacency matrix of the districts.
#'   }}
#'  \item{population_hungary}{
#'    A numeric matrix containing the population per 10000 inhabitants of each district over time.}
#' }
#' @details
#' This dataset contains chickenpox counts in the 20 districts (NUTS 3) of Hungary over a time period of 522 weeks (from 2005 to 2014).
#'
#' The row-normalized adjacency matrix indicates which districts share a common border.
#'
#' The population data is only availabyle on a yearly basis and has been linearly interpolated by us to obtain weekly estimates.
#'
#'
#' @source
#' The data originate from the UCI Machine Learning Repository and the 
#' Hungarian Central Statistical Office and are licensed under the 
#' Creative Commons Attribution 4.0 International (CC BY 4.0) license:
#'   * \url{https://archive.ics.uci.edu/dataset/580/hungarian+chickenpox+cases}
#'   * \url{https://www.ksh.hu/stadat_files/nep/en/nep0034.html}
#'
#' @examples
#' dat <- load_data("chickenpox")
#' chickenpox <- dat$chickenpox
#' population_hungary <- dat$population_hungary
#' W_hungary <- dat$W_hungary
#'
#' covariates <- list(population = population_hungary, 
#'                  season_cos = SpatialConstant(cos(2 * pi / 52 * 1:522)),
#'                  season_sin = SpatialConstant(sin(2 * pi / 52 * 1:522)))
#' glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
#'           covariates = covariates, family = vpoisson("log"))
#' glmstarma(chickenpox, list(past_obs = 1), wlist = W_hungary, 
#'           covariates = covariates, family = vnegative.binomial("log"))
#'
#' mean_model <- list(past_obs = rep(1, 7))
#' dispersion_model <- list(past_obs = 1)
#' dglmstarma(chickenpox, mean_model, dispersion_model, mean_family = vquasipoisson("log"), 
#'            dispersion_link = "log", wlist = W_hungary, mean_covariates = covariates)
#' @docType data
#' @name chickenpox
#' @keywords datasets
"chickenpox"


