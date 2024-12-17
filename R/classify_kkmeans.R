#' Classify new data based on kkmeans result
#' @name cluster_new
#'
#' @description Given training data, test data, and kkmeans-results,
#'  get which partitions a new set of data belong to
#'
#' @param data data vector
#' @param kern the kernel to use, one of ('gaussian', 'poly'), can use first
#' letter
#' @param param parameter to pass to kernel function.
#' @export
cluster_new <- function(test_data, train_data, train_result) {
  if (!is.matrix(test_data))
    test_data <- as.matrix(test_data)
  if (!is.matrix(train_data))
    train_data <- as.matrix(train_data)
  if (!inherits(train_result, "kkmeans_result"))
    stop("`train_result` must be a result from calling `kkmeans::kkmeans()`")
  if (nrow(train_data) != length(train_result$cluster))
    stop("Number of rows in `train_data` must match length of result in `train_result`")

  test_labels <- .Call(
    "classify_kkmeans",
    test_data, train_data, train_result$cluster - 1L,
    train_result$kern, train_result$params[1],
    train_result$params[2], train_result$k
  )
  return(test_labels)
}
