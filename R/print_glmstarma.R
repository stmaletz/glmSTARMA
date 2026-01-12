# -----------------------------------------------------------------------------
# File: print_glmstarma.R
# Purpose: Implement print functions for glmSTARMA models
# Author: Steffen Maletz
# Last modified: 2025-11-14
# -----------------------------------------------------------------------------

#' @exportS3Method base::print
print.glmstarma <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  # Hilfsfunktion für den wiederholten Blockdruck
  print_coef_block <- function(name, coef_matrix, model_matrix) {
    if (!is.null(coef_matrix) && nrow(coef_matrix) > 0) {
      cat("\n", name, " Coefficients:\n", sep = "")
      max_width <- max(nchar(formatC(coef_matrix, digits = digits, format = "f")))
      mat <- formatC(coef_matrix, digits = digits, format = "f", flag = "", width = max_width)
      mat[model_matrix == 0] <- "."
      colnames(mat) <- colnames(model_matrix)
      rownames(mat) <- rownames(model_matrix)
      print.default(mat, print.gap = 2L, quote = FALSE)
    }
  }

  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")

  if (length(stats::coef(x)) > 0) {
    coeffs <- stats::coef(x, TRUE)
    
    cat("Estimated Coefficients:\n")
    cat("Intercept:\n")
    print.default(format(drop(coeffs$intercept), digits = digits), print.gap = 2L, quote = FALSE)

    print_coef_block("Autoregressive", coeffs$past_obs, x$model$past_obs)
    print_coef_block("Moving-Average", coeffs$past_mean, x$model$past_mean)
    print_coef_block("Covariate", coeffs$covariates, x$model$covariates)

  } else {
    cat("No coefficients\n")
  }

  invisible(x)
}

#' @exportS3Method base::print
print.dglmstarma <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  print_coef_block <- function(coef_list, model_list, digits, title){
    title_length <- nchar(title)
    left_padding <- floor((50 - title_length) / 2)
    right_padding <- ceiling((50 - title_length) / 2)
    cat(strrep("=", 50), "\n")
    cat(strrep(" ", left_padding), title, strrep(" ", right_padding), "\n", sep = "")
    cat(strrep("=", 50), "\n\n")
    cat("Intercept:\n")
    print(format(drop(coef_list$intercept), digits = digits), quote = FALSE)

    print_subblock <- function(name, coef_matrix, model_matrix){
      if(!is.null(coef_matrix) && nrow(coef_matrix) > 0){
        cat("\n", name, ":\n", sep = "")
        max_width <- max(nchar(formatC(coef_matrix, digits = digits, format = "f")))
        mat <- formatC(coef_matrix, digits = digits, format = "f", flag = "", width = max_width)
        mat[model_matrix == 0] <- "."
        colnames(mat) <- colnames(model_matrix)
        rownames(mat) <- rownames(model_matrix)
        print.default(mat, print.gap = 2L, quote = FALSE)
      }
    }

    print_subblock("Autoregressive Coefficients", coef_list$past_obs, model_list$past_obs)
    print_subblock("Moving-Average Coefficients", coef_list$past_mean, model_list$past_mean)
    print_subblock("Covariate Coefficients", coef_list$covariates, model_list$covariates)

    cat("\n")
  }

  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n")

  if(length(x$mean$coefficients) > 0){
    print_coef_block(x$mean$coefficients_list, x$mean$model, digits, "Mean Model:")
  } else {
    cat("No mean coefficients\n\n")
  }

  if(length(x$dispersion$coefficients) > 0){
    print_coef_block(x$dispersion$coefficients_list, x$dispersion$model, digits, "Dispersion Model:")
  } else {
    cat("No dispersion coefficients\n")
  }

  invisible(x)
}







