# -----------------------------------------------------------------------------
# File: load_data.R
# Purpose: Methods to load example datasets from the glmSTARMA GitHub repository
# Author: Steffen Maletz
# Last modified: 2026-01-14
# -----------------------------------------------------------------------------


#' @title Load example datasets
#'
#' @description Download and return datasets from the glmSTARMA GitHub repository 
#'
#' @param name Name of the dataset to load. One of \code{"rota"}, \code{"chickenpox"}, or \code{"sst"}.
#' @param refresh Logical; re-download the dataset if it already exists locally.
#' @param directory Directory where the dataset should be cached. Defaults to a user-specific data directory of the \code{glmSTARMA} package.
#'
#' @return A named list of objects
#' @details
#' This function downloads example datasets from the glmSTARMA GitHub repository and caches them in a directory specified by the user.
#' The default directory is a user-specific data directory of the \code{glmSTARMA} package.
#' If the dataset has already been downloaded to the specified directory, it is loaded from the local cache unless \code{refresh = TRUE} is specified.
#'
#' @seealso \code{\link{delete_glmSTARMA_data}}, \code{\link{rota}}, \code{\link{chickenpox}}, \code{\link{sst}} 
#' @examples
#' \donttest{
#' # Load the 'chickenpox' dataset
#' chickenpox_data <- load_data("chickenpox", directory = tempdir())
#' str(chickenpox_data)
#' }
#' @export
load_data <- function(name = NULL, refresh = FALSE, directory = tools::R_user_dir("glmSTARMA", which = "data")) {
  stopifnot("Only one dataset can be loaded at a time" = length(name) == 1,
            "Parameter 'refresh' must be of type logical" = is.logical(refresh),
            "Parameter 'refresh' must be of length 1" = length(refresh) == 1,
            "name must be in 'rota', 'chickenpox', or 'sst')" = name %in% c("rota", "chickenpox", "sst"))

  # GitHub raw URL
  base_url <- "https://raw.githubusercontent.com/stmaletz/glmSTARMA/main/data-raw"
  file_name <- paste0(name, ".rda")
  url <- paste0(base_url, "/", file_name)

  # user data directory
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  data_file <- file.path(directory, file_name)

  # Download only if file is missing or refresh is requested
  if (!file.exists(data_file) || isTRUE(refresh)) {
    tryCatch(
      utils::download.file(url, data_file, mode = "wb", quiet = TRUE),
      error = function(e) {
        if (!file.exists(data_file)) {
          stop(
            "Data could not be downloaded.\n",
            "No internet connection or file does not exist.\n\n",
            "URL: ", url,
            call. = FALSE
          )
        }
      }
    )
  }

  # Load the .rda file into an isolated environment
  env <- new.env(parent = emptyenv())
  objs <- load(data_file, envir = env)

  mget(objs, envir = env)
}



#' @title Delete cached example datasets
#' @description Delete one or more cached example datasets downloaded via \code{load_data()}.
#' @param name Name(s) of the dataset(s) to delete. One or more of \code{"rota"}, \code{"chickenpox"}, or \code{"sst"}. If \code{NULL} (default), no action is taken.
#' @return Invisibly returns \code{TRUE} if all specified datasets were deleted, \code{FALSE} otherwise.
#' @details
#' This function deletes datasets that were previously downloaded and cached using \code{load_data()}
#' from the user-specific data directory.
#' If no datasets are found in the cache, a message is printed and no action is taken.
#' @seealso \code{\link{load_data}}, \code{\link{rota}}, \code{\link{chickenpox}}, \code{\link{sst}} 
#' @examples
#' delete_glmSTARMA_data("chickenpox")  # Only gives a message if dataset is not cached
#' @export
delete_glmSTARMA_data <- function(name = NULL) {
    data_dir <- tools::R_user_dir("glmSTARMA", which = "data")

    if(!dir.exists(data_dir)){
        message("There is no dataset to delete.")
        return(invisible(FALSE))
    }
    files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
    if(length(files) == 0){
        message("There is no dataset to delete.")
        return(invisible(FALSE))
    }

    if(!is.null(name)){
        stopifnot("You can only delete the datasets downloaded with functions of the glmSTARMA package" = all(name %in% c("rota", "chickenpox", "sst")))
        rtrn_val <- TRUE
        for(nam in name){
            target <- file.path(data_dir, paste0(nam, ".rda"))

            if(!file.exists(target)){
                message("The dataset '", nam, "' is not in the cache.", call. = FALSE)
                rtrn_val <- FALSE
            }
            unlink(target)
            msg <- paste0("Deleted: ", nam)
            message(msg)
        }
        return(invisible(rtrn_val))
    }
}
