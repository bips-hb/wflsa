#' Using a random graph with the flsa package

library(flsa)

set.seed(1)

# The parameters
p <- 200
density <- .5
lambda1 <- 1
lambda2 <- 1
  
#' Create a connListObj for FLSA package given Adjacency Matrix
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

# Create the random binary adjacency matrix
W <- random_binary_weight_matrix(p, density)
  
# Create the corresponding connListObj needed by the flsa package
connListObj <- create_connListObj(W)
  
# generate the raw data 
y <- rnorm(p)
  
flsa::flsa(y, lambda1 = lambda1, lambda2 = lambda2,
           connListObj = connListObj$connList)

