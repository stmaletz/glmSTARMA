# -----------------------------------------------------------------------------
# File: generateW.R
# Purpose: Implement function to generate spatial weight matrices for simulation
# Author: Steffen Maletz
# Last modified: 2026-08-12
# -----------------------------------------------------------------------------

#' @rdname generateW
#' @title Generate Spatial Weight Matrices for Simulation
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
  
  result_list <- do.call(toCall, args) 
  
  structure(result_list, 
            class = c("generateW", "list"), 
            method = method,
            n_locations = dim,
            width = width)
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

#' @exportS3Method base::print
print.generateW <- function(x, ...) {
  cat("Spatial Weight Matrices\n")
  cat("Spatial domain:", attr(x, "method"), "\n")
  cat("Number of locations:", attr(x, "n_locations"), "\n")
  if(attr(x, "method") == "rectangle"){
    cat("Width of grid:", attr(x, "n_locations") / attr(x, "width") ,"x" , attr(x, "width"), "\n")
  }
  cat("Number of matrices: ", length(x), ifelse(length(x) == 1, "matrix\n", "matrices\n"))
  invisible(x)
}

#' @exportS3Method base::summary
summary.generateW <- function(object, ...) {
  densities <- sapply(object, function(m) sum(m != 0) / length(m))

  res <- list(
    method = attr(object, "method"),
    n_locations = attr(object, "n_locations"),
    length = length(object),
    densities = densities
  )
  
  class(res) <- "summary.generateW"
  return(res)
}

#' @exportS3Method base::print
print.summary.generateW <- function(x, ...) {
  cat("Summary of Spatial Weight Matrices:\n")
  cat("- Spatial domain:", x$method, "\n")
  cat("- Number of locations:", x$n_locations, "\n")
  cat("- Number of matrices:", x$length, "\n")
  cat("- Non-zero proportions:\n")
  print(round(x$densities, 4))
  invisible(x)
}



#' @rdname plot.generateW
#' @title Plot Spatial Weight Matrices for Simulation
#'
#' @description This function plots spatial weight matrices and neighborhoods generated by the \code{generateW} function. 
#'
#' @param x An object of class \code{generateW}, returned by the \code{generateW} function.
#' @param order (integer scalar) Specifies which spatial order of the neighborhood to plot. Default is \code{0}.
#' @param type (character scalar) Specifies the type of plot. Options are "matrix" for a heatmap representation of the weight matrix, and "network" for a network graph representation of the neighborhood structure. Default is "matrix".
#' @param col (character vector) Color palette for the heatmap. Default is a blue gradient. Ignored if \code{type} is "network".
#' @param node_col (character scalar) Color of the nodes in the network graph. Default is "black". Ignored if \code{type} is "matrix".
#' @param edge_col (character scalar) Color of the edges in the network graph. Default is "gray". Ignored if \code{type} is "matrix".
#' @param ... Additional arguments passed to the underlying plotting functions.
#' @details 
#' Plotting neighborhood structures of type "independent" or "full" is not meaningful, as these structures do not represent spatial dependencies in a conventional sense. 
#' For "independent", each location is treated as an isolated node, while "full" represents a complete graph where every location is connected to every other location, but each connection is treated independently.
#' @examples
#' W_rect <- generateW(method = "rectangle", dim = 100, maxOrder = 2, width = 5)
#' plot(W_rect, order = 1, type = "matrix")
#' plot(W_rect, order = 1, type = "network")
#'
#' @exportS3Method base::plot
plot.generateW <- function(x, order = 0, type = c("matrix", "network"), 
                          col = NULL, node_col = "black", edge_col = "gray", ...) {
  type <- match.arg(type)
  elem <- order + 1
  if (elem > length(x) || elem < 1) {
    stop("Requested spatial order exceeds the number of available matrices.")
  }
  
  mat <- x[[elem]]
  method <- attr(x, "method")
  dim <- attr(x, "n_locations")
  
  if (type == "matrix") {
    if (is.null(col)) {
      col <- colorRampPalette(c("transparent", "darkblue"))(100)
    } else if(is.character(col) && length(col) == 1) {
      col <- colorRampPalette(c("transparent", col))(100)
    }
    mat_rot <- t(mat[nrow(mat):1, ])

    image(1:ncol(mat), 1:nrow(mat), mat_rot,
          main = paste("Weight Matrix - Method:", method, "- Order:", order),
          xlab = "Columns", ylab = "Rows",
          col = col, axes = FALSE, ...)
    box()
    
  } else if (type == "network") {
    if (method == "line") {
      x_coords <- 1:dim
      y_coords <- rep(0, dim)
      
    } else if (method == "rectangle") {
      width <- attr(x, "width")
      height <- dim / width
      
      coords <- expand.grid(x = seq_len(width), y = seq_len(height))
      
      x_coords <- coords$x
      y_coords <- max(coords$y) - coords$y + 1
      
    } else {
      angles <- seq(0, 2 * pi, length.out = dim + 1)[-(dim + 1)]
      x_coords <- cos(angles)
      y_coords <- sin(angles)
    }
    
    plot(x_coords, y_coords, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = paste("Neighborhood Graph - Type:", method, "- Order:", order), ...)
    edges <- which(mat > 0, arr.ind = TRUE)
    
    if (nrow(edges) > 0) {
      segments(x0 = x_coords[edges[, "row"]], 
               y0 = y_coords[edges[, "row"]],
               x1 = x_coords[edges[, "col"]], 
               y1 = y_coords[edges[, "col"]],
               col = edge_col, lwd = 1.5)
    }
    points(x_coords, y_coords, pch = 19, col = node_col, cex = 2)
    
    if (dim <= 50) {
      text(x_coords, y_coords, labels = 1:dim, col = "white", cex = 0.7)
    }
  }
  
  invisible(x)
}

