#' An Efficient Kernel K-Means Algorithm
#' @name kkmeans
#'
#' @description Performs kernel k-means with the specified kernel using an
#' algorithm similar to Hartigan and Wong's k-means algorithm.
#'
#' @param data data to cluster
#' @param k the number of clusters
#' @param kern the kernel to use, one of ('gaussian', 'poly'), can use first
#' letter
#' @param params parameters to pass to kernel function.
#' @export
kkmeans <- function(data, k, kern = "g", param = 1, iter_max = 1000L) {

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

  retlist <- .Call('kkmeans', data, k, kern, param, iter_max)

  names(retlist) <- c("cluster", "centers", "sse")
  return(retlist)
}
