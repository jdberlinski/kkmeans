#' An Efficient Kernel K-Means Algorithm
#' @name kkmeans
#'
#' @description Performs kernel k-means with the specified kernel using an
#' algorithm similar to Hartigan and Wong's k-means algorithm.
#'
#' @param data data to cluster
#' @param k the number of clusters. if jump == TRUE, then this is ignored
#' @param kern the kernel to use, one of ('gaussian', 'poly'), can use first
#' letter
#' @param param value of parameter to pass to kernel function.(eg sigma in
#' gaussian kernel)
#' @export
kkmeans <- function(data, k = FALSE, kern = "g", param = 1, jump = FALSE,
                    k_max = FALSE, eta = FALSE, iter_max = 1000L) {

  valid_kerns = c("gaussian", "poly")
  valid_prefs = c("g", "p")

  # some light error checking
  if ( !(kern %in% valid_kerns) & !(kern %in% valid_prefs) )
    stop(paste0("`kern` must be one of ", paste(valid_kerns, collapse = ", "), "."))
  if ( !is.integer(k) )
    k <- as.integer(k)
  if ( !is.matrix(data) ) {
    warning("Converting data to matrix.")
    data <- as.matrix(data)
  }
  if ( !is.integer(iter_max) )
    iter_max <- as.integer(iter.max)

  # TODO: I think that the jump statistic should be its own function, since
  # these aren't returning the same type of list.
  if (!jump) {
    retlist <- .Call('kkmeans', data, k, kern, param, iter_max)
    names(retlist) <- c("cluster", "centers", "sse")
  }
  else {
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
        jump_vals[i] <- (wss[i] / (n * p)) ^ (-eta) - 0^(-eta)
      else
        jump_vals[i] <- (wss[i] / (n * p)) ^ (-eta) - (wss[i - 1] / (n * p)) ^ (-eta)
    }
    retlist$jump_stat <- jump_vals
  }

  return(retlist)
}
