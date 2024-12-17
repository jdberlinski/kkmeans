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
cluster_new <- function(test_data, train_data, train_labels, kern = "g", param = 1, param2 = 1) {
  test_labels <- .Call(
    "classify_kkmeans",
    test_data, train_data, train_labels - 1L,
    kern, param, param2, max(train_labels)
  )
  return(test_labels)
}
