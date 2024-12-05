#' An Efficient Kernel K-Means Algorithm
#' @name kkmeans
#'
#' @description Performs kernel k-means with the specified kernel using an
#' optimal-transfer quick-transfer algorithm.
#'
#' @param data Numeric data to cluster. This will be converted to a matrix using `as.matrix`.
#' @param k Number of clusters.
#' @param kern Kernel to use, one of ('gaussian', 'poly').
#' @param param value of parameter to pass to kernel function.(eg sigma in
#' gaussian kernel). The Gaussian kernel is K(x, y) = exp(- ||x - y||^2 / (2*`param`))),
#' and the polynomial kernel is K(x, y) = (x'y + 1) ^ `param`
#' @param nstart Number of times to run the algorithm. The run with the lowest
#' total within cluster SSE (in feature space) will be returned
#' @param iter_max The maximum number of iterations to allow.
#' @param estimate If using the Gaussian kernel, specifying `estimate = "mknn"` will use
#' an `nn`-nearest neighbor method for estimating `param`.
#' @param nn How many neighbors to consider for mknn estimation.
#' @param init_centers The initial values for cluster membership. If `nstart` is greater
#' than 1, any start beyond the first iteration will use randomized centers.
#' @param method Which method to use for kernel k-means iteration. One of ("otqt", "macqueen", "lloyd").
#' "otqt" is a method using optimal-transfer and quick-transfer heuristics similar to the Hartigan and
#' Wong algorithm for k-means clustering.
#' @param trueest Whether or not the within-cluster sum of squares should be
#' recomputed in R after clustering is finished
#' @param kmat kernel matrix, if using a custom kernel
#' @return A list containing the following useful information
#' \describe{
#'   \item{cluster}{The final cluster membership.}
#'   \item{centers}{A k x p matrix, the rows of which contain the centers of the clusters in R^n (not to be confused
#' with the clusters in feature space)}
#'   \item{wss}{The within-cluster sum of squares for each cluster in feature space.}
#'   \item{param}{The parameter value used.}
#' }
#' @export
#' @examples
#' data <- as.matrix(iris[, 1:4])
#'
#' # cluster using linear kernel (normal k-means)
#' result <- kkmeans(data, k = 3, kern = "poly", param = 1)
#'
#' # cluster using gaussian kernel
#' # estimating the parameter with 3-nearest neighbors
#' result <- kkmeans(data, k = 3, kern = "g", estimate = "mknn", nn = 3)
kkmeans <- function(data, k, kern = "g", param = 1, nstart = 10, iter_max = 1000L, estimate = F,
                    nn = 0, init_centers = sample(1:k, size = nrow(data), replace = TRUE),
                    method = c("otqt", "macqueen", "lloyd", "ot"), trueest = F, kmat = NULL, random_centers = TRUE) {

  valid_kerns = c("gaussian", "poly")
  valid_prefs = c("g", "p")
  valid_methods = c("otqt", "macqueen", "lloyd", "ot")

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
    # warning("Converting data to matrix.")
    data <- as.matrix(data)
  }
  if ( !is.integer(iter_max) )
    iter_max <- as.integer(iter_max)

  # if (depth > 0 && kern != "g" && kern != "gaussian")
  #   stop("If depth > 0, `kern` must be gaussian.")
  use_centers <- FALSE

  if (length(init_centers) == k) {
    candidate_centers <- init_centers
    init_centers <- rep(0, nrow(data))
    nstart <- 1
    random_centers <- FALSE
    use_centers <- TRUE
  } else if (!is.integer(init_centers) || min(init_centers) > 0) {
    init_centers <- as.integer(init_centers - 1)
  } else {
    init_centers <- init_centers - 1L
  }
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
  if (is.null(kmat))
    K <- get_kernel_matrix(data, kern, param)
  else
    K <- kmat

  for (i in 1:nstart) {
    # if random centers, then choose k points to be the initial centers, and
    # assign other points to their closest one
    if (random_centers | use_centers) {
      if (random_centers)
        candidate_centers <- sample(nrow(data), k)
      init_centers[candidate_centers] <- 0:(k - 1)
      for (j in seq_len(nrow(data))) {
        if (j %in% candidate_centers) next
        dvec <- numeric(k)
        for (l in seq_len(k))
          dvec[l] <- K[j, j] + K[candidate_centers[l], candidate_centers[l]] - 2*K[j, candidate_centers[l]]
        init_centers[j] <- which.min(dvec) - 1
      }
      init_centers <- as.integer(init_centers)
    }

    retlist <- .Call('kkmeans', data, k, kern, param, iter_max, init_centers, method, K)
    # fix WSS?
    if (trueest) {
      wss_vals <- Map(function(k) {
                        inclust <- which(retlist[[1]] == k)
                        Ks <- K[inclust, inclust]
                        nk <- length(inclust)
                        J <- matrix(1, nrow = nk, ncol = nk)
                        C <- Ks - J %*% Ks / nk - Ks %*% J / nk + J %*% Ks %*% J / nk^2
                        sum(diag(C))
                      },
                      1:k)
      retlist[[3]] <- Reduce(c, wss_vals)
    }
    # re-randomize clusters if running more than once
    if (!random_centers)
      init_centers <- as.integer(sample(1:k, size = nrow(data), replace = TRUE) - 1)
    if (sum(retlist[[3]]) < lowest_wss) {
      lowest_wss <- sum(retlist[[3]])
      lowest_res <- retlist
    }
  }
  names(lowest_res) <- c("cluster", "centers", "wss", "niter")
  lowest_res$param <- param
  # }


  return(lowest_res)
}
