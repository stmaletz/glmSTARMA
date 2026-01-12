# -----------------------------------------------------------------------------
# File: generateW.R
# Purpose: Implement function to generate spatial weight matrices for simulation
# Author: Steffen Maletz
# Last modified: 2025-11-16
# -----------------------------------------------------------------------------

#' @rdname generateW
#' @title Generate spatial weight matrices for simulation
#'
#' @description This function generates row-normalized spatial weight matrices for different types of neighborhood structures.
#'
#' @param method (character scalar) Defines type of neighborhood structure. Options are "rectangle", "line", "circle", "full", and "independent". Default is "rectangle".
#' @param dim (integer scalar) Number of locations, i.e. dimension of the time series.
#' @param maxOrder (integer scalar) Maximum spatial order up to which the spatial weight matrices are generated. Ignored if `method`` is "full" or "independent".
#' @param width (integer scalar) Width of the rectangular grid. Must be a divisor of `dim`. Ignored if `method` is not "rectangle".
#' @param ... Additional arguments passed to specific methods.
#' @return A list of (row normalized) spatial weight matrices.
#' @details 
#' The function generates spatial weight matrices for different types of neighborhood structures. The options are:
#' * "rectangle" - A regular rectangular grid (2 dimensional) with `width` columns and `dim / width` rows. The spatial order is defined by the Euclidean distances between locations.
#' * "line" - Locations are placed on a line (1 dimensional). The spatial order is defined by the Euclidean distances between locations.
#' * "circle" - Locations are placed on a circle. The spatial order is defined by the Euclidean distances between locations. In contrast to the "line" neighborhood, there are no boundary locations.
#' * "full" - Generates a list with `dim^2` matrices. Allows simulation/fitting of a full time series model without any restrictions in dependencies between the locations. Not recommended if `dim` is large.
#' * "independent" - Generates a list with `dim` matrices. Each matrix is a spatial weight matrix with a single 1 in the diagonal. Allows simultaneously simulation/fitting of `dim` univariate time series models without spatial dependencies.
#'
#' @examples
#' generateW(method = "rectangle", dim = 100, maxOrder = 2, width = 5)
#' generateW(method = "full", dim = 4)
#' @references
#' For more advanced spatial weight matrices, consider using the \code{spdep} package.
#' - Bivand R, Pebesma E, Gómez-Rubio V (2013). *Applied spatial data analysis with R*, Second edition. Springer, NY. \url{https://asdar-book.org/}.
#' - Pebesma E, Bivand R (2023). *Spatial Data Science With Applications in R*. Chapman & Hall. \url{https://r-spatial.org/book/}.
#' @export
generateW <- function(method = c("rectangle", "line", "circle", "full", "independent"),
                      dim,
                      maxOrder = NULL,
                      width = NULL,
                      ...) {
  method <- match.arg(method)
  toCall <- get(paste0("generateW.", method))
  
  # Build argument list
  args <- list(dim = dim, ...)
  if (!is.null(maxOrder)) args$maxOrder <- maxOrder
  if (!is.null(width)) args$width <- width
  
  do.call(toCall, args)
}



generateW.full <- function(dim, ...) {
  W <- vector(mode = "list", length = dim^2)
  for (i in seq_len(dim^2)) {
    vec <- numeric(dim^2)
    vec[i] <- 1
    W[[i]] <- matrix(vec, dim, dim)
  }
  W
}



generateW.independent <- function(dim, ...) {
  W <- vector("list", dim)
  idx <- 1 - (dim + 1)
  for (i in seq_len(dim)) {
    vec <- numeric(dim^2)
    idx <- idx + (dim + 1)
    vec[idx] <- 1
    W[[i]] <- matrix(vec, dim, dim)
  }
  W
}


generateW.line <- function(dim, maxOrder, ...) {
  stopifnot(!missing(maxOrder))
  if (maxOrder < 0){
    warning("maxOrder is negative. Returning empty list.")
    return(list())
  }
  if (maxOrder == 0){
    return(list(diag(dim)))
  }
  if (maxOrder > floor(dim / 2))
  {
    stop("maxOrder must not be larger than floor(dim/2), as no higher order neighbors exist.")
  }

  distances <- as.matrix(stats::dist(seq_len(dim)))
  result <- vector("list", maxOrder + 1)

  for (l in seq(from = 0, to = maxOrder, by = 1L)) {
    M <- 1 * (distances == l)
    result[[l + 1]] <- M / rowSums(M)
  }
  result
}


generateW.circle <- function(dim, maxOrder, ...) {
  stopifnot(!missing(maxOrder))
  if (maxOrder < 0){
    warning("maxOrder is negative. Returning empty list.")
    return(list())
  }
  if (maxOrder == 0){
    return(list(diag(dim)))
  }
  if (maxOrder > floor(dim / 2)){
    warning("maxOrder should not exceed floor(dim/2). As these results in duplicate matrices, consider using maxOrder = floor(dim/2).")
  }
  result <- vector("list", maxOrder + 1)
  result[[1]] <- diag(dim)

  for (l in seq(from = 1, to = maxOrder, by = 1L)) {
    M <- matrix(0, dim, dim)
    for (i in seq_len(dim)) {
      idx <- ((i + c(-l, l) - 1) %% dim) + 1
      M[i, idx] <- ifelse(diff(idx) == 0, 1, 0.5)
      # indexes <- c(-l, l) + i
      # indexes <- ifelse(indexes < 1, indexes + dim, indexes)
      # indexes <- ifelse(indexes > dim, indexes - dim, indexes)
      # result[[l + 1]][i, indexes] <- ifelse(diff(indexes) == 0, 1, 0.5)
    }
    result[[l + 1]] <- M
  }
  result
}

generateW.rectangle <- function(dim, width, maxOrder, ...) {
  stopifnot(!missing(width), !missing(maxOrder))
  if (maxOrder < 0){
    warning("maxOrder is negative. Returning empty list.")
    return(list())
  }
  if (maxOrder == 0){
    return(list(diag(dim)))
  }

  if (dim %% width != 0)
    stop("'dim' must be width * k")

  height <- dim / width
  
  coords <- expand.grid(x = seq_len(width), y = seq_len(height))
  distances <- as.matrix(stats::dist(coords)^2)
  uni_dist <- sort(unique(as.vector(distances)))

  result <- vector("list", maxOrder + 1)
  for (l in seq(from = 0, to = maxOrder, by = 1L)) {
    M <- 1 * (distances == uni_dist[l + 1])
    result[[l + 1]] <- M / rowSums(M)
  }
  result
}


