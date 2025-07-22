#' Get the kernel matrix for a dataset
#' @name get_kernel_matrix
#'
#' @description Given a dataset, kernel function, and tuning parameter, will
#' return the n x n kernel matrix
#'
#' @param data data vector
#' @param kern the kernel to use, one of ('gaussian', 'poly', 'sigmoid', 'laplacian'), can use first
#' letter
#' @param param1 first parameter to pass to kernel function.
#' @param param2 second parameter to pass to kernel function.
#' @export
get_kernel_matrix <- function(data, kern = "g", param = 1, param2 = 1) {

  valid_kerns = c("gaussian", "poly", "sigmoid", "laplacian")
  valid_prefs = c("g", "p", "s", "l")

  # some light error checking
  if ( !(kern %in% valid_kerns) & !(kern %in% valid_prefs) )
    stop(paste0("`kern` must be one of ", paste(valid_kerns, collapse = ", "), "."))
  if ( !is.matrix(data) ) {
    # warning("Converting data to matrix.")
    data <- as.matrix(data)
  }

  kernel_matrix <- .Call('get_k_matrix', data, kern, param, param2)

  return(kernel_matrix)
}
