#' @docType package
#' @name glmSTARMA-package
#' @details
#' The implemented models are based on spatio-temporal autoregressive moving average (STARMA) models. They incorporate spatial and temporal dependencies by spatial lagging, via spatial weight matrices, and temporal lagging via past observations and past values of the linear predictor.
#'
#' The main functions for fitting such models are \code{\link{glmstarma}} and \code{\link{dglmstarma}}.
#' The main difference between the two functions is that \code{glmstarma} fits a model for the (conditional) mean of the spatio-temporal process and \code{dglmstarma} fits two models, one for the (conditional) mean and another one for the (conditional) dispersion.
#' The mean model in both functions generalizes the structure of spatio-temporal Poisson autoregressions, and allows for various distributions from the exponential dispersion family.
#' The dispersion model can be seen as an generalization of an spatio-temporal GARCH or log-GARCH model.
#' Data can be simulated with \code{\link{glmstarma.sim}} and \code{\link{dglmstarma.sim}}.
#'
#' For more details on the models see the documentation of the fitting functions \code{\link{glmstarma}} and \code{\link{dglmstarma}}.
#'
#' @references
#' - Armillotta, M., Tsagris, M., & Fokianos, K. (2024). Inference for Network Count Time Series with the R Package PNAR. The R Journal, 15(4), 255–269. \doi{10.32614/RJ-2023-094}
#' - Barreto‐Souza, W., Piancastelli, L. S., Fokianos, K., & Ombao, H. (2025). Time‐Varying Dispersion Integer‐Valued GARCH Models. Journal of Time Series Analysis. \doi{10.1111/jtsa.12838}
#' - Cliff, A. D., & Ord, J. K. (1975). Space-Time Modelling with an Application to Regional Forecasting. Transactions of the Institute of British Geographers, 64, 119–128. \doi{10.2307/621469}
#' - Jahn, M., Weiß, C.H., Kim, H.Y. (2023), Approximately linear INGARCH models for spatio-temporal counts, Journal of the Royal Statistical Society Series C: Applied Statistics, 72(2), 476-497, \doi{10.1093/jrsssc/qlad018}
#' - Jørgensen, B. (1987), Exponential Dispersion Models. Journal of the Royal Statistical Society: Series B (Methodological), 49(2), 127-145. \doi{10.1111/j.2517-6161.1987.tb01685.x}
#' - Knight, M., Leeming, K., Nason, G., & Nunes, M. (2020). Generalized Network Autoregressive Processes and the GNAR Package. Journal of Statistical Software, 96(5), 1–36. \doi{10.18637/jss.v096.i05}
#' - Maletz, S., Fokianos, K., & Fried, R. (2024). Spatio-Temporal Count Autoregression. Data Science in Science, 3(1). \doi{10.1080/26941899.2024.2425171}
#' - Meyer, S., Held, L., & Höhle, M. (2017). Spatio-Temporal Analysis of Epidemic Phenomena Using the R Package surveillance. Journal of Statistical Software, 77(11), 1–55. \doi{10.18637/jss.v077.i11}
#' - Otto, P. (2024). A multivariate spatial and spatiotemporal ARCH Model. Spatial Statistics, 60. \doi{10.1016/j.spasta.2024.100823}
#' - Pfeifer, P. E., & Deutsch, S. J. (1980). A Three-Stage Iterative Procedure for Space-Time Modeling Phillip. Technometrics, 22(1), 35–47. \doi{10.2307/1268381}
#' - Smyth, G.K. (1989), Generalized Linear Models with Varying Dispersion. Journal of the Royal Statistical Society: Series B (Methodological), 51(1), 47-60. \doi{10.1111/j.2517-6161.1989.tb01747.x}
#' @import nloptr
#' @importFrom Rcpp evalCpp
#' @useDynLib glmSTARMA, .registration = TRUE
"_PACKAGE"