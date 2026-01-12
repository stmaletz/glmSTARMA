# -----------------------------------------------------------------------------
# File: family.R
# Purpose: Implements family functions for glmSTARMA models
# Author: Steffen Maletz
# Last modified: 2026-01-10
# -----------------------------------------------------------------------------


#' @title Families for spatio-temporal GLMs
#' @name stfamily
#' @aliases vpoisson vquasipoisson vnegative.binomial vbinomial vquasibinomial vgamma vinverse.gaussian vnormal
#' @description
#' These functions create family objects for various distributions used in
#' spatio-temporal generalized linear models (STGLMs). 
#' Each function returns an object of class \code{"stfamily"} describing a
#' (conditional) marginal distribution, a link function, and optional dispersion values.
#' The output is intended only for use within this package and is processed internally in the functions of this package.
#'
#' @param link Character string specifying the link function.
#'   Options depend on the distribution (see Details).
#' @param dispersion Optional dispersion parameter(s). Can be either a numerical scalar describing a global, 
#' time-invariant dispersion parameter, a vector with (temporal) constant dispersion parameters for each location, 
#' or a matrix with dispersion parameters for each location and time. (Rows are locations, columns are time points)
#' If \code{NULL}, dispersion will be estimated as a scalar where applicable. 
#' @param const Optional numeric constant used in some link functions.
#' @param size Number of trials for binomial-type families.
#' @param copula Optional copula family to model dependence between responses. Has no effect on parameter estimation, but only on situations in which data is generated with this family.
#' @param copula_param Parameter for the copula. (Numeric scalar of length 1)
#' @param sampling_method Sampling algorithm for Poisson and quasi-Poisson.
#'
#' @details
#' The \code{link} argument specifies the link function to be used for the family. The available link functions depend on the distribution:
#'  \itemize{
#'   \item Poisson, Quasi-Poisson, Negative-Binomial: "log", "identity", "sqrt", "softplus"
#'   \item Binomial, Quasi-Binomial: "softclipping", "identity", "logit", "probit"
#'   \item Gamma: "inverse", "log", "identity"
#'   \item Inverse Gaussian: "1/mu^2", "inverse", "identity", "log"
#'   \item Gaussian/Normal: "identity", "log", "inverse"
#' }
#'
#'
#' The following families are available:
#' \itemize{
#'   \item \code{vpoisson()} – Poisson distribution
#'   \item \code{vquasipoisson()} – Quasi-Poisson, i.e. Poisson like, but dispersion can differ from 1
#'   \item \code{vnegative.binomial()} – Negative binomial distribution
#'   \item \code{vbinomial()} – Binomial distribution
#'   \item \code{vquasibinomial()} – Quasi-Binomial, i.e. Binomial like, but dispersion can differ from 1
#'   \item \code{vgamma()} – Gamma distribution
#'   \item \code{vinverse.gaussian()} – Inverse Gaussian distribution
#'   \item \code{vnormal()} – Gaussian distribution
#'   \item \code{vgarch()} – GARCH distribution
#' }
#'
#' The following copulas are available:
#' \itemize{
#'   \item "normal" – Gaussian copula
#'   \item "t" – t copula
#'   \item "clayton" – Clayton copula
#'   \item "frank" – Frank copula
#'   \item "gumbel" – Gumbel copula
#'   \item "joe" – Joe copula
#' }
#' The data generating processes of each distribution rely on (sequences of) uniform marginals, which are transformed to obtain the observed data.
#' The copula specifies the dependence structure between these uniform marginals. If no copula (\code{copula = NULL}) is specified, the uniform marginals are generated independent.
#'
#' For most distributions, the data generation is based on inversion of the cumulative distribution function (CDF).
#' For the Poisson distribution, two different sampling methods are implemented:
#' \code{sampling_method = "inversion"} results in the inversion method using \code{qpois} and \code{sampling_method = "poisson_process"} implements the Poisson process method described in Fokianos et al. (2020).
#' For a quasi-Poisson model, four different sampling methods are implemented: \code{sampling_method = "build_up"}, \code{sampling_method = "chop_down"}, \code{sampling_method = "branching"} and \code{sampling_method = "negbin"}.
#' The first three methods generate data from a generalized Poisson distribution, see Consul and Jain (1973) and are described in Demirtas (2017). We translated the code from the R package \code{RNGRNGforGPD} to C++ for these three methods. The \code{sampling_method = "negbin"} uses the Inversion method on properly parameterized negative binomial distribution to generate the data.
#' If the "branching" or "negbin" method is used, only overdispersion can be generated. Dispersion values resulting in underdispersion will be set to 1, i.e. a standard Poisson case.
#' For Quasi-Binomial models, in case of overdispersion, the data is generated from a sequence of positive correlated Bernoulli trials, see Ahn and Chen (1995) for a discussion. In case of underdispersion, data is generated using a normal approximation.
#' In case of the inverse Gaussian distribution, data is generated using the Michael-Schucany-Haas method, see Michael et al. (1976).
#'
#' Note that for the negative binomial family, the dispersion parameter corresponds to the shape parameter of the negative binomial distribution.
#'
#' @return An object of class \code{"stfamily"} containing elements
#'   such as \code{link}, \code{distribution}, \code{variance}, \code{dev.resids}.
#'
#' @references
#' - Ahn, H., & Chen, J. J. (1995). Generation of Over-Dispersed and Under-Dispersed Binomial Variates. Journal of Computational and Graphical Statistics, 4(1), 55–64. \doi{10.1080/10618600.1995.10474665}
#' - Consul, P. C., & Jain, G. C. (1973). A Generalization of the Poisson Distribution. Technometrics, 15(4), 791–799. \doi{10.1080/00401706.1973.10489112}
#' - Demirtas, H. (2017). On accurate and precise generation of generalized Poisson variates. Communications in Statistics - Simulation and Computation, 46(1), 489–499. \doi{10.1080/03610918.2014.968725}
#' - Fokianos, K., Støve, B., Tjøstheim, D., & Doukhan, P. (2020). Multivariate count autoregression. Bernoulli, 26(1), 471–499. \doi{10.3150/19-BEJ1132}
#' - Michael, J. R., Schucany, W. R., & Haas, R. W. (1976). Generating Random Variates Using Transformations with Multiple Roots. The American Statistician, 30(2), 88–90. \doi{10.1080/00031305.1976.10479147}
#' @examples
#' fam <- vpoisson(link = "log")
#' print(fam)
#'
#' fam2 <- vbinomial(link = "logit", size = 10)
#' print(fam2)
#'
#' @seealso \code{\link[stats]{family}}
#'
#' @export
vpoisson <- function(link = c("log", "identity", "sqrt", "softplus"), const = 1, copula = NULL, copula_param = NULL, sampling_method = c("inversion", "poisson_process")){
  fam <- list()
  fam$link <- match.arg(link)
  stopifnot("const must be a positive scalar" = (fam$link %in% c("identity", "sqrt")) || (is.numeric(const) && length(const) == 1 && const > 0),
            "const must not be infinite" = (fam$link %in% c("identity", "sqrt")) || (!is.infinite(const)),
            "copula_param must be a numeric scalar" = is.null(copula_param) || (is.numeric(copula_param) && length(copula_param) == 1),
            "copula_param must be positive" = is.null(copula_param) || copula_param > 0,
            "copula_param must be specified when copula is specified" = is.null(copula) || !is.null(copula_param))
  fam$distribution <- "poisson"
  fam$non_negative_parameters <- FALSE
  if(fam$link %in% c("identity", "sqrt")){
    fam$non_negative_parameters <- TRUE
  }
  fam$const <- ifelse(is.null(const), 1, const)
  fam$sampling_method <- match.arg(sampling_method)
  fast_sampling <- fam$sampling_method == "inversion"
  if(!is.null(copula)){
    fam$copula <- match.arg(copula, c("normal", "t", "clayton", "frank", "gumbel", "joe"))
    fam$copula_param <- copula_param
    fam$fast_sampling <- fast_sampling
  } else {
    fam$copula <- NULL
    fam$copula_param <- NULL
    fam$fast_sampling <- NULL
  }
  fam$dispersion <- 1
  fam$estimate_dispersion <- FALSE
  fam$variance <- function(mu, dispersion, ignore_dispersion = FALSE){
    mu
  }
  fam$dev.resids <- function (y, mu, dispersion, ignore_dispersion = FALSE) 
  {
    r <- mu
    p <- which(y > 0)
    r[p] <- ((y * log(y/mu) - (y - mu)))[p]
    2 * r
  }
  class(fam) <- "stfamily"
  return(fam)
}

#' @rdname stfamily
#' @export
vquasipoisson <- function(link = c("log", "identity", "sqrt", "softplus"), dispersion = NULL, const = 1, copula = NULL, copula_param = NULL, sampling_method = c("build_up", "chop_down", "branching", "negbin")){
  fam <- list()
  fam$link <- match.arg(link)
  stopifnot("const must be a positive scalar" = (fam$link %in% c("identity", "sqrt")) || (is.numeric(const) && length(const) == 1 && const > 0),
            "const must not be infinite" = (fam$link %in% c("identity", "sqrt")) || (!is.infinite(const)),
            "copula_param must be a numeric scalar" = is.null(copula_param) || (is.numeric(copula_param) && length(copula_param) == 1),
            "copula_param must be positive" = is.null(copula_param) || copula_param > 0,
            "copula_param must be specified when copula is specified" = is.null(copula) || !is.null(copula_param),
            "dispersion must be positive" = is.null(dispersion) || (is.numeric(dispersion) && all(dispersion > 0)),
            "dispersion must not contain infinite values" = is.null(dispersion) || !any(is.infinite(dispersion)))
  fam$distribution <- "quasipoisson"
  fam$non_negative_parameters <- FALSE
  fam$sampling_method <- match.arg(sampling_method)
  if(fam$link %in% c("identity", "sqrt")){
    fam$non_negative_parameters <- TRUE
  }
  fam$const <- ifelse(is.null(const), 1, const)
  if(!is.null(copula)){
    fam$copula <- match.arg(copula, c("normal", "t", "clayton", "frank", "gumbel", "joe"))
    fam$copula_param <- copula_param
  } else {
    fam$copula <- NULL
    fam$copula_param <- NULL
  }
  if(is.null(dispersion)){
    fam$dispersion <- 1
    fam$estimate_dispersion <- TRUE
  } else {
    fam$dispersion <- dispersion
    fam$estimate_dispersion <- FALSE
    if(any(fam$dispersion < 1) && fam$sampling_method %in% c("branching", "negbin")){
      warning("The 'branching' and 'negbin' sampling methods can only generate overdispersed data. In the (theoretical) case of underdispersion a Poisson distribution is used instead (dispersion = 1).")
    }
  }

  fam$variance <- function(mu, dispersion, ignore_dispersion = FALSE){
    mu * ifelse(ignore_dispersion, 1, dispersion)
  }
  fam$dev.resids <- function (y, mu, dispersion, ignore_dispersion = FALSE) 
  {
    r <- mu
    p <- which(y > 0)
    r[p] <- ((y * log(y/mu) - (y - mu)))[p]
    if(ignore_dispersion) {
      return(2 * r)
    } else {
      return(2 * r / dispersion)
    }
  }
  class(fam) <- "stfamily"
  return(fam)
}

#' @rdname stfamily
#' @export
vnegative.binomial <- function(link = c("log", "identity", "sqrt", "softplus"), dispersion = NULL, const = 1, copula = NULL, copula_param = NULL){
  fam <- list()
  fam$link <- match.arg(link)
  stopifnot("const must be a positive scalar" = (fam$link %in% c("identity", "sqrt")) || (is.numeric(const) && length(const) == 1 && const > 0),
            "const must not be infinite" = (fam$link %in% c("identity", "sqrt")) || (!is.infinite(const)),
            "copula_param must be a numeric scalar" = is.null(copula_param) || (is.numeric(copula_param) && length(copula_param) == 1),
            "copula_param must be positive" = is.null(copula_param) || copula_param > 0,
            "copula_param must be specified when copula is specified" = is.null(copula) || !is.null(copula_param),
            "dispersion must be positive" = is.null(dispersion) || (is.numeric(dispersion) && all(dispersion > 0)),
            "dispersion must not contain infinite values" = is.null(dispersion) || !any(is.infinite(dispersion)))
  fam$distribution <- "negative_binomial"
  if(is.null(dispersion)){
    fam$dispersion <- 0
    fam$estimate_dispersion <- TRUE
  } else {
    fam$dispersion <- dispersion
    fam$estimate_dispersion <- FALSE
  }
  fam$non_negative_parameters <- (fam$link %in% c("identity", "sqrt"))
  fam$const <- const
  if(!is.null(copula)){
    fam$copula <- match.arg(copula, c("normal", "t", "clayton", "frank", "gumbel", "joe"))
    fam$copula_param <- copula_param
  } else {
    fam$copula <- NULL
    fam$copula_param <- NULL
  }
  fam$variance <- function(mu, dispersion, ignore_dispersion = FALSE){
    mu + dispersion * mu^2 * (1 - ignore_dispersion)
  }
  fam$dev.resids <- function (y, mu, dispersion, ignore_dispersion = FALSE) 
  {
    theta <- dispersion * 0^(ignore_dispersion)
    2 * (y * log(pmax(1, y)/mu) - (y + theta) * log((y + theta)/(mu + theta)))
  }
  class(fam) <- "stfamily"
  return(fam)
}

#' @rdname stfamily
#' @export
vbinomial <- function(link = c("softclipping", "identity", "logit", "probit"), size = 1, const = 1, copula = NULL, copula_param = NULL){
  fam <- list()
  fam$link <- match.arg(link)
  stopifnot("const must be a positive scalar" = (fam$link != "softclipping") || (is.numeric(const) && length(const) == 1 && const > 0),
            "const must not be infinite" = (fam$link != "softclipping") || (!is.infinite(const)),
            "copula_param must be a numeric scalar" = is.null(copula_param) || (is.numeric(copula_param) && length(copula_param) == 1),
            "copula_param must be positive" = is.null(copula_param) || copula_param > 0,
            "copula_param must be specified when copula is specified" = is.null(copula) || !is.null(copula_param),
            "size must be a positive integer" = (is.numeric(size) && is.vector(size) && all(size > 0) && all(size == floor(size))),
            "size must not contain infinite values" = is.null(size) || !any(is.infinite(size)))
  
  fam$distribution <- "binomial"
  fam$non_negative_parameters <- FALSE
  fam$dispersion <- 1
  fam$estimate_dispersion <- FALSE
  fam$const <- ifelse(is.null(const), 1, const)
  fam$non_negative_parameters <- (fam$link == "identity")
  fam$size <- size
  if(!is.null(copula)){
    fam$copula <- match.arg(copula, c("normal", "t", "clayton", "frank", "gumbel", "joe"))
    fam$copula_param <- copula_param
  } else {
    fam$copula <- NULL
    fam$copula_param <- NULL
  }

  fam$variance <- function(mu, dispersion, ignore_dispersion = FALSE){
      p <- mu / size
      size * p * (1 - p)
  }
  fam$dev.resids <- function (y, mu, dispersion, ignore_dispersion = FALSE) 
  {
    r <- 2 * (y * log(y / mu) + (size - y) * log((size - y) / (size - mu)))
    y_0 <- which(y == 0)
    y_n <- which(y == size)
    r[y_0] <- 2 * size * log((size - mu) / (size - mu))
    r[y_n] <- 2 * size * log(size / mu)
    r
  }

  class(fam) <- "stfamily"
  return(fam)
}

#' @rdname stfamily
#' @export
vquasibinomial <- function(link = c("softclipping", "identity", "logit", "probit"), size = 1, dispersion = NULL, const = 1, copula = NULL, copula_param = NULL){
  fam <- list()
  fam$link <- match.arg(link)
  stopifnot("const must be a positive scalar" = (fam$link != "softclipping") || (is.numeric(const) && length(const) == 1 && const > 0),
            "const must not be infinite" = (fam$link != "softclipping") || (!is.infinite(const)),
            "copula_param must be a numeric scalar" = is.null(copula_param) || (is.numeric(copula_param) && length(copula_param) == 1),
            "copula_param must be positive" = is.null(copula_param) || copula_param > 0,
            "copula_param must be specified when copula is specified" = is.null(copula) || !is.null(copula_param),
            "size must be a positive integer" = (is.numeric(size) && is.vector(size) && all(size > 0) && all(size == floor(size))),
            "size must not contain infinite values" = (is.null(size) || !any(is.infinite(size))),
            "dispersion must be positive" = is.null(dispersion) || (is.numeric(dispersion) && all(dispersion > 0)),
            "dispersion must not contain infinite values" = is.null(dispersion) || !any(is.infinite(dispersion)))

  fam$distribution <- "quasibinomial"
  if(is.null(dispersion)){
    fam$dispersion <- 0
    fam$estimate_dispersion <- TRUE
  } else {
    fam$dispersion <- dispersion
    fam$estimate_dispersion <- FALSE
    if(is.matrix(dispersion)){
      stopifnot("dispersion must have the same number of rows as the length of size" = (length(size) == 1) || nrow(dispersion) == length(size))
    } else if(is.vector(dispersion)){
      stopifnot("dispersion must have the same length as size" = (length(size) == 1) || length(dispersion) == length(size))
    }
    if(any(dispersion > size)) {
        warning("The dispersion parameter for quasibinomial models should not exceed the number of trials (size).")
      } 
  }
  fam$const <- ifelse(is.null(const), 1, const)
  fam$non_negative_parameters <- (fam$link == "identity")
  fam$size <- size
  if(!is.null(copula)){
    fam$copula <- match.arg(copula, c("normal", "t", "clayton", "frank", "gumbel", "joe"))
    fam$copula_param <- copula_param
  } else {
    fam$copula <- NULL
    fam$copula_param <- NULL
  }
  fam$variance <- function(mu, dispersion, ignore_dispersion = FALSE){
      p <- mu / size
      size * p * (1 - p) * dispersion^(1 - ignore_dispersion)
  }
  fam$dev.resids <- function (y, mu, dispersion, ignore_dispersion = FALSE) 
  {
    r <- 2 * (y * log(y / mu) + (size - y) * log((size - y) / (size - mu)))
    y_0 <- which(y == 0)
    y_n <- which(y == size)
    r[y_0] <- 2 * size * log((size - mu) / (size - mu))
    r[y_n] <- 2 * size * log(size / mu)
    if(ignore_dispersion){
      return(r)
    } else {
      return(r/dispersion)
    }
  }

  class(fam) <- "stfamily"
  return(fam)
}

#' @rdname stfamily
#' @export
vgamma <- function(link = c("inverse", "log", "identity"), dispersion = NULL, copula = NULL, copula_param = NULL){
  stopifnot("copula_param must be a numeric scalar" = is.null(copula_param) || (is.numeric(copula_param) && length(copula_param) == 1),
            "copula_param must be positive" = is.null(copula_param) || copula_param > 0,
            "copula_param must be specified when copula is specified" = is.null(copula) || !is.null(copula_param),
            "dispersion must be positive" = is.null(dispersion) || (is.numeric(dispersion) && all(dispersion > 0)),
            "dispersion must not contain infinite values" = is.null(dispersion) || !any(is.infinite(dispersion)))
  fam <- list()
  fam$link <- match.arg(link)
  fam$distribution <- "gamma"
  fam$non_negative_parameters <- (fam$link %in% c("identity", "inverse"))

  if(is.null(dispersion)){
    fam$dispersion <- 1.0
    fam$estimate_dispersion <- TRUE
  } else {
    fam$dispersion <- dispersion
    fam$estimate_dispersion <- FALSE
  }
  if(!is.null(copula)){
    fam$copula <- match.arg(copula, c("normal", "t", "clayton", "frank", "gumbel", "joe"))
    fam$copula_param <- copula_param
  } else {
    fam$copula <- NULL
    fam$copula_param <- NULL
  }
  fam$variance <- function(mu, dispersion, ignore_dispersion = FALSE){
      mu^2 * dispersion^(1 - ignore_dispersion)
  }
  fam$dev.resids <- function(y, mu, dispersion, ignore_dispersion = FALSE) 
  {
    r <- -2  * (log(ifelse(y == 0, 1, y/mu)) - (y - mu)/mu)
    if(ignore_dispersion){
      return(r)
    } else {
      return(r/dispersion)
    }
  }
  class(fam) <- "stfamily"
  return(fam)
}

#' @rdname stfamily
#' @export
vinverse.gaussian <- function(link = c("1/mu^2", "inverse", "identity", "log"), dispersion = NULL, copula = NULL, copula_param = NULL){
  stopifnot("copula_param must be a numeric scalar" = is.null(copula_param) || (is.numeric(copula_param) && length(copula_param) == 1),
            "copula_param must be positive" = is.null(copula_param) || copula_param > 0,
            "copula_param must be specified when copula is specified" = is.null(copula) || !is.null(copula_param),
            "dispersion must be positive" = is.null(dispersion) || (is.numeric(dispersion) && all(dispersion > 0)),
            "dispersion must not contain infinite values" = is.null(dispersion) || !any(is.infinite(dispersion)))
  fam <- list()
  fam$link <- match.arg(link)
  fam$distribution <- "inverse_gaussian"
  fam$non_negative_parameters <- (fam$link %in% c("1/mu^2", "identity", "inverse"))

  if(is.null(dispersion)){
    fam$dispersion <- 1.0
    fam$estimate_dispersion <- TRUE
  } else {
    fam$dispersion <- dispersion
    fam$estimate_dispersion <- FALSE
  }
  if(!is.null(copula)){
    fam$copula <- match.arg(copula, c("normal", "t", "clayton", "frank", "gumbel", "joe"))
    fam$copula_param <- copula_param
  } else {
    fam$copula <- NULL
    fam$copula_param <- NULL
  }

  fam$variance <- function(mu, dispersion, ignore_dispersion = FALSE){
      mu^3 * dispersion^(1 - ignore_dispersion)
  }
  fam$dev.resids <- function (y, mu, dispersion, ignore_dispersion = FALSE) 
  {
    r <- ((y - mu)^2)/(y * mu^2)
    if(ignore_dispersion){
      return(r)
    } else {
      return(r/dispersion)
    }
  }

  class(fam) <- "stfamily"
  return(fam)
}

#' @rdname stfamily
#' @export
vnormal <- function(link = c("identity", "log", "inverse"), dispersion = NULL, copula = NULL, copula_param = NULL){
  stopifnot("copula_param must be a numeric scalar" = is.null(copula_param) || (is.numeric(copula_param) && length(copula_param) == 1),
            "copula_param must be positive" = is.null(copula_param) || copula_param > 0,
            "copula_param must be specified when copula is specified" = is.null(copula) || !is.null(copula_param),
            "dispersion must be positive" = is.null(dispersion) || (is.numeric(dispersion) && all(dispersion > 0)),
            "dispersion must not contain infinite values" = is.null(dispersion) || !any(is.infinite(dispersion)))
  fam <- list()
  fam$link <- match.arg(link)
  fam$distribution <- "gaussian"
  fam$non_negative_parameters <- FALSE
  if(is.null(dispersion)){
    fam$dispersion <- 1.0
    fam$estimate_dispersion <- TRUE
  } else {
    fam$dispersion <- dispersion
    fam$estimate_dispersion <- FALSE
  }
  if(!is.null(copula)){
    fam$copula <- match.arg(copula, c("normal", "t", "clayton", "frank", "gumbel", "joe"))
    fam$copula_param <- copula_param
  } else {
    fam$copula <- NULL
    fam$copula_param <- NULL
  }

  fam$variance <- function(mu, dispersion, ignore_dispersion = FALSE){
      if(ignore_dispersion){
        return(matrix(1, nrow = nrow(mu), ncol = ncol(mu)))
      } else {
        return(dispersion)
      }
  }
  fam$dev.resids <- function (y, mu, dispersion, ignore_dispersion = FALSE) 
  {
    r <- (y - mu)^2
    if(ignore_dispersion){
      return(r)
    } else {
      return(r/dispersion)
    }
  }

  class(fam) <- "stfamily"
  return(fam)
}

#' @method print stfamily
#' @exportS3Method base::print
#' @keywords internal
print.stfamily <- function(x, ...){
  cat("\nMarginal distribution:", x$distribution, "\n")
  cat("Link function:", x$link, "\n\n")
  invisible(x)
}

