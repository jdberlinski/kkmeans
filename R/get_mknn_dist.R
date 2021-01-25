#' Get the average distance to each points k-nearest neighbor
#' @name get_mknn_dist
#'
#' @description Given a dataset and a value k, will return the value of the
#' average distance from each point to it's k-nearest neighbor
#'
#' @param data the data vector
#' @param k which neighbor to average over
#' @export
get_mknn_dist <- function(data, k = F) {
  if (!k) k <- round(log2(nrow(data)) + 1)
  dm <- as.matrix(dist(data, upper = T)) ^ 2
  knn_dist <- apply(dm, 1, function(x) x[order(x)[k + 1]])
  return(mean(knn_dist))
}
