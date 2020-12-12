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
#' @param depth either 0 (no parameter estimation) or the maximum iteration
#' depth for parameter estimation
#' @param iter_max the maximum number of iterations to allow
#' @param estimate which form of estimation to do -- "otg" or "mknn"
#' @param nn which neighbor to consider for mknn
#' @export
kkmeans <- function(data, k, kern = "g", param = 1, nstart = 10, iter_max = 1000L, estimate = F,
                    nn = 0, init_centers = sample(1:k, size = nrow(data), replace = TRUE),
                    method = c("otqt", "macqueen")) {

  valid_kerns = c("gaussian", "poly")
  valid_prefs = c("g", "p")
  valid_methods = c("otqt", "macqueen")

  # some light error checking
  if ( !(kern %in% valid_kerns) & !(kern %in% valid_prefs) )
    stop(paste0("`kern` must be one of ", paste(valid_kerns, collapse = ", "), "."))
  if ( !is.integer(k) )
    k <- as.integer(k)
  if ( !is.integer(nstart) )
    nstart <- as.integer(nstart)
  # if ( !is.integer(depth) )
  #   depth <- as.integer(depth)
  if ( !is.matrix(data) ) {
    warning("Converting data to matrix.")
    data <- as.matrix(data)
  }
  if ( !is.integer(iter_max) )
    iter_max <- as.integer(iter_max)

  # if (depth > 0 && kern != "g" && kern != "gaussian")
  #   stop("If depth > 0, `kern` must be gaussian.")

  if (!is.integer(init_centers))
    init_centers <- as.integer(init_centers - 1)
  else
    init_centers <- init_centers - 1L
  if (max(init_centers) > k - 1 || min(init_centers) < 0)
    stop("Initial centers must be between 1 and k")

  if (!(method[[1]] %in% valid_methods))
    stop(paste("`method` must be one of", paste(valid_methods, collapse = ", "), "."))


  if (tolower(estimate) == "mknn") {
    if (!nn) nn <- round(log2(nrow(data)) + 1)
    param <- get_mknn_dist(data, nn)
  }

  # if (depth > 0) {
  #   lowest_res <- .Call('kkmeans_est', data, k, depth, param, iter_max)
  #   names(lowest_res) <- c("cluster", "centers", "wss", "param")
  # }
  # else {
  method <- which(valid_methods == method[[1]])
  lowest_wss <- Inf
  lowest_res <- NULL
  for (i in 1:nstart) {
    retlist <- .Call('kkmeans', data, k, kern, param, iter_max, init_centers, method)
    if (sum(retlist[[3]]) < lowest_wss) {
      lowest_wss <- sum(retlist[[3]])
      lowest_res <- retlist
    }
  }
  names(lowest_res) <- c("cluster", "centers", "wss")
  lowest_res$param <- param
  # }


  return(lowest_res)
}
