#' Function to get jump statistic for varying values of k
#' @name jump_stat
#'
#' @description
#'
#' @param data
#' @param kern
#' @param param
#' @param k_max
#' @param eta
#' @param iter_max
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
