#' Function to get jump statistic for varying values of k
#' @name jump_stat
#'
#' @description Obtains the jump statistic for a particular kernel for the
#' specified number of clusters
#'
#' @param data Numeric data to cluster. This will be converted to a matrix using `as.matrix`.
#' @param kern The kernel to use.
#' @param param The parameter value to pass to the kernel
#' @param k_max The maximum number of clusters to consider
#' @param eta Power for the jump statistic
#' @param iter_max Maximum number of iterations to use in `kkmeans` call
#' @return Sum of squares and value of jump statistic for 1, ..., K chosen
#' clusters
#' @export
jump_stat <- function(data, kern = "g", param = 1, k_max, eta, iter_max = 1000L) {
  valid_kerns = c("gaussian", "poly")
  valid_prefs = c("g", "p")

  # some light error checking
  if ( !(kern %in% valid_kerns) & !(kern %in% valid_prefs) )
    stop(paste0("`kern` must be one of ", paste(valid_kerns, collapse = ", "), "."))
  if ( !is.matrix(data) ) {
    warning("Converting data to matrix.")
    data <- as.matrix(data)
  }

  n <- nrow(data)
  p <- ncol(data)
  k_vals <- seq(k_max)
  min_wss <- Inf
  wss <- numeric(k_max)
  jump_vals <- numeric(k_max)

  retlist <- NULL
  for (curr_k in k_vals) {
    curr_out <- .Call('kkmeans', data, curr_k, kern, param, iter_max)
    wss[curr_k] <- sum(curr_out[[3]])
  }
  for (i in k_vals) {
    if (i == 1)
      jump_vals[i] <- (wss[i] / (n * p)) ^ (-eta)
    else
      jump_vals[i] <- (wss[i] / (n * p)) ^ (-eta) - (wss[i - 1] / (n * p)) ^ (-eta)
  }
  retlist$sse <- wss
  retlist$jump_stat <- jump_vals
  return(retlist)
}
