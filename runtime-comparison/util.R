#' Util. Function for Runtime Analysis of the FLSA and wFLSA Package
#'
#' This script contains various utility functions used to check the runtime 
#' of the FLSA and WFLSA packages. It includes functions for creating a 
#' connListObj for the FLSA package, copying the upper triangle of a square 
#' matrix to its lower triangle, and generating a band matrix with a band around 
#' the diagonal.
#' 

#' Create a connListObj for FLSA package
#' 
#' This function creates a connListObj for the FLSA package given a weight matrix.
#' 
#' @param W A weight matrix.
#' @return A connListObj for FLSA package and a logical value indicating whether 
#'         all elements of connList are NULL.
create_connListObj <- function(W) {
  m <- nrow(W)
  
  connList <- vector("list", m)
  class(connList) <- "connListObj"
  
  lapply(1:m, function(i) {
    indices <- which(W[i, ] != 0)
    if (length(indices) != 0) {
      connList[[i]] <<- as.integer(which(W[i, ] != 0) - 1)
    }
  })
  
  names(connList) <- as.character(0:(m-1))
  
  connList_is_null <- all(sapply(connList, function(l) is.null(l)))
  
  return(list(connList = connList, connList_is_null = connList_is_null))
}


#' Copy Upper Triangle to Lower Triangle
#'
#' Copies the upper triangle of a square matrix to its lower triangle.
#'
#' @param mat A square matrix.
#' @return A matrix with the upper triangle copied to the lower triangle.
#' @export
#'
#' @examples
#' mat <- matrix(1:9, nrow = 3)
#' copy_upper_to_lower(mat)
copy_upper_to_lower <- function(mat) {
  p <- nrow(mat)
  for (i in 1:(p - 1)) {
    mat[(i + 1):p, i] <- mat[i, (i + 1):p]
  }
  return(mat)
}

#' Band Matrix Generator
#'
#' \code{band_matrix} function creates a square matrix with a band around the diagonal.
#'
#' @param p Size of the square matrix.
#' @return A square matrix with a band around the diagonal.
#' @usage band_matrix(p)
#' Default band width is set to 1, but it can be adjusted as needed.
#'
#' @examples
#' \dontrun{
#' # Generate a band matrix of size 5
#' band_matrix(5)
#' }
band_matrix <- function(p) {
  
  # Create a matrix with a band around the diagonal
  my_matrix <- matrix(0, nrow = p, ncol = p)
  
  for (i in 1:p) {
    lower <- max(1, i - 1)
    upper <- min(p, i + 1)
    my_matrix[i, lower:upper] <- 1
  }
  
  diag(my_matrix) <- 0  # Set diagonal elements to 0 
  my_matrix <- copy_upper_to_lower(my_matrix)
  return(my_matrix)
}

#' Generate Random Binary Weight Matrix
#'
#' \code{random_binary_weight_matrix} function creates a random binary weight matrix 
#' of size \code{p} with a specified density.
#'
#' @param p Size of the square weight matrix.
#' @param density Density of the binary values (probability of being 1).
#' @return A random binary weight matrix of size \code{p} with the specified density.
#'
#' @examples
#' \dontrun{
#' # Generate a random binary weight matrix of size 5 with density 0.3
#' random_binary_weight_matrix(5, 0.3)
#' }
random_binary_weight_matrix <- function(p, density) {
  W <- matrix(rbinom(p*p, 1, density), nrow = p)
  W[upper.tri(W)] <- t(W)[upper.tri(W)]
  diag(W) <- 0
  return(W)
}
