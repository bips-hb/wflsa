#' Apply soft thresholding to a numeric vector
#'
#' This function applies soft thresholding to each element of a numeric vector
#' using a specified threshold value (\code{lambda1}). Soft thresholding shrinks
#' values toward zero by applying a penalty when the absolute value of an element
#' exceeds the threshold.
#'
#' @param x A numeric vector.
#' @param lambda1 Threshold value for soft thresholding.
#' @return A numeric vector with each element softened using the specified threshold.
#' @examples
#' x <- c(3, -2, 5, -4, 1)
#' lambda1 <- 2
#' soft_thresholded <- soft_threshold(x, lambda1)
#' print(soft_thresholded)
#' @export
soft_threshold <- function(x, lambda1) {
  return(sign(x) * pmax(abs(x) - lambda1, 0))
}
