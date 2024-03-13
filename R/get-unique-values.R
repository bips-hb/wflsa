#' Get unique absolute values of a vector after rounding
#'
#' This function takes a numeric vector, calculates its absolute values,
#' rounds them to a specified number of digits, finds the unique values,
#' sorts them, and returns the sorted unique absolute values.
#'
#' @param vector A numeric vector.
#' @param digits Number of digits to round to (default is 10).
#' @return A numeric vector containing the sorted unique absolute values.
#' @examples
#' vector <- c(1.23456789, -1.23456788, 1.23456787, -1.23456786, 1.23456785, -1.23456784)
#' sorted_unique_abs_vals <- get_sorted_unique_abs_values(vector)
#' print(sorted_unique_abs_vals)
#' @export
get_unique_values <- function(vector, digits = 10) {
  
  abs_vector <- abs(vector)  # Take the absolute values of the vector
  
  rounded_abs_vector <- round(abs_vector, digits = digits)  # Round the absolute values
  
  sorted_unique_abs_values <- sort(unique(rounded_abs_vector))  # Sort and find unique values
  
  return(sorted_unique_abs_values)
}
