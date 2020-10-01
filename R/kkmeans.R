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
#' @param nstart number of times to run the algorithm. the run with the lowest
#' total within cluster SSE (in feature space) will be returned
#' @param iter_max the maximum number of iterations to allow
#' @export
kkmeans <- function(data, k, kern = "g", param = 1, nstart = 10, iter_max = 1000L) {

  valid_kerns = c("gaussian", "poly")
  valid_prefs = c("g", "p")

  # some light error checking
  if ( !(kern %in% valid_kerns) & !(kern %in% valid_prefs) )
    stop(paste0("`kern` must be one of ", paste(valid_kerns, collapse = ", "), "."))
  if ( !is.integer(k) )
    k <- as.integer(k)
  if ( !is.integer(nstart) )
    nstart <- as.integer(nstart)
  if ( !is.matrix(data) ) {
    warning("Converting data to matrix.")
    data <- as.matrix(data)
  }
  if ( !is.integer(iter_max) )
    iter_max <- as.integer(iter.max)

  lowest_wss <- Inf
  lowest_res <- NULL
  for (i in 1:nstart) {
    retlist <- .Call('kkmeans', data, k, kern, param, iter_max)
    if (sum(retlist[[3]]) < lowest_wss) {
      lowest_wss <- sum(retlist[[3]])
      lowest_res <- retlist
    }
  }
  names(lowest_res) <- c("cluster", "centers", "wss")

  return(lowest_res)
}
