#' @title Rota Virus Infections in Germany
#'
#' @description Multivariate count time series with rota virus infections in the counties of Germany.
#'
#' @usage
#' data("rota")
#'
#' @format
#' \describe{
#'   \item{rota}{
#'     A matrix with counts of rota virus infections (rows = counties, columns = time points).
#'   }
#'   \item{gdr_feature}{
#'     A numeric vector indicating whether a county was part of the former German Democratic Republic (1) or not (0). Counties in the Capital Berlin have the value 0.5.
#'   }
#'  \item{population_germany}{
#'    A numeric matrix containing the population in 10000 inhabitants of each district over time.}
#'   \item{W_germany}{
#'     A list of matrices containing spatial weight matrices:
#'     \enumerate{
#'      \item Identity matrix.  
#'      \item Row-normalized adjacency matrix of the districts. See details.
#'      \item Row-normalized adjacency matrix of order 2. See details.
#'   }}
#' }
#'
#' @details
#' This dataset contains weekly rota virus counts in 411 counties of Germany over a time period of 1252 weeks (from 2001 to 2024).
#' Estimates of the population size of each county are only available on a yearly basis and have been linearly interpolated to obtain weekly estimates.
#'
#' Throughout the observation period, there were several territorial reforms in which certain regions were merged or reorganized. The data source (RKI) for infection counts takes into account the regional divisions that were in existence at the time the data was queried from the website.
#' The data on population size contains the data for the divisions valid at the time of data collection. For the purpose of standardisation, the data was aggregated in order to reflect the territorial reforms as accurately as possible and to align it with the infection figures.
#' For Berlin, annual population estimates are only available for Berlin as a whole. We use data from the 2022 census to divide the population proportionally among the 12 districts. 
#'
#' The row-normalized adjacency matrix of first order indicates which districts share a common border. The weight matrix of order 2 indicates which districts can be reached in two steps (i.e., via a common neighbor).
#'
#' @source
#' The data originate from the Robert Koch Institute (RKI) and the Federal Statistical Office of Germany (Destatis).
#'
#' **Infection data:** Robert Koch Institute: SurvStat@RKI 2.0, 
#' \url{https://survstat.rki.de}, retrieved on 2025-02-03.
#' 
#' **Population data:** Federal Statistical Office of Germany (Destatis), 
#' Data licence Germany – attribution – version 2.0, 
#' \url{https://www-genesis.destatis.de/datenbank/online/statistic/12411/table/12411-0015/table-toolbar}, 
#' retrieved on 2025-12-08.
#' 
#' **Census data (Berlin):** Federal Statistical Office of Germany (Destatis), 
#' Data licence Germany – attribution – version 2.0, 
#' \url{https://ergebnisse.zensus2022.de/datenbank/online/statistic/1000A/table/1000A-0000}, 
#' retrieved on 2025-12-08.
#' @examples
#' \dontrun{
#' data("rota")
#'
#' covariates <- list(population = population_germany,
#'                   gdr = TimeConstant(gdr_feature),
#'                   season_cos = SpatialConstant(cos(2 * pi / 52 * 1:1252)),
#'                   season_sin = SpatialConstant(sin(2 * pi / 52 * 1:1252)),
#'                   vaccine_west = (gdr_feature == 0) %*% t(seq(ncol(rota)) >= 654),
#'                   vaccine_east = (gdr_feature > 0) %*% t(seq(ncol(rota)) >= 654))
#' fit <- glmstarma(rota, list(past_obs = rep(2, 4)), wlist = W_germany,
#'                covariates = covariates, family = vpoisson("log"))
#'
#'
#' mean_model <- list(past_obs = rep(2, 4))
#' dispersion_model <- list(past_obs = 2)
#' fit2 <- dglmstarma(rota, mean_model, dispersion_model, mean_covariates = covariates,
#'                    dispersion_covariates = covariates,
#'                    mean_family = vquasipoisson("log"), 
#'                     dispersion_link = "log", W_germany)
#' }  
#' @docType data
#' @name rota
#' @keywords datasets
"rota"


#' Internal auxiliary data for rota
#'
#' @name population_germany
#' @seealso \code{\link{rota}}
#' @docType data
#' @keywords internal
"population_germany"

#' Internal auxiliary data for rota
#'
#' @name gdr_feature
#' @seealso \code{\link{rota}}
#' @docType data
#' @keywords internal
"gdr_feature"
#' Internal auxiliary data for rota
#'
#' @name W_germany
#' @seealso \code{\link{rota}}
#' @docType data
#' @keywords internal
"W_germany"