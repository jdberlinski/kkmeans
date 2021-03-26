#' Estimate the bandwidth parameter for a gaussian kernel using MATr
#' @name matr
#'
#' @description given data and number of clusters K, choose a bandwidth from a
#' grid according to the MATR algorithm
#'
#' @param dat data vector
#' @param k number of clusters
#' @param grid parameter grid to search on
#' @param tol tolerence for choosing the "best" sigma values. Reducing the
#' tolerance will give larger values of the bandwidth parameter
#' @export
matr <- function(dat, k, grid, tol = 0.01) {
  l <- numeric(length(grid))
  rough1 <- quantile(dist(dat), .1)
  rough2 <- get_mknn_dist(dat)
  rough <- max(rough1, rough2)
  S <- get_kernel_matrix(dat, 'g', rough)
  for (i in seq_along(grid)) {
    res <- kkmeans(dat, k, 'g', grid[i])
    Z <- model.matrix(~ 0 + factor(res$cluster))
    X <- Z %*% solve(crossprod(Z)) %*% t(Z)
    l[i] <- sum(diag(crossprod(S, X)))
  }
  best_inds <- which(max(l)*(1 - tol) <=  l)
  best <- grid[best_inds[length(best_inds)]]
  return(list(l = l, best = best))
}
