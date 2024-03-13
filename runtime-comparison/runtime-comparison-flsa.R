#' Runtime Assessment for FLSA and wFLSA Packages
#'
#' This script is used to assess the runtime of the FLSA and wFLSA packages for 
#' both an arbitrary binary matrix and the 1D case. It performs a time comparison 
#' between the two packages using the 'microbenchmark' package. Additionally, it 
#' utilizes utility functions defined in the 'util.R' script to facilitate the 
#' assessment process.

# Load necessary libraries
library(flsa)
library(wflsa)
library(dplyr)
library(microbenchmark)
library(progress)
library(ggplot2)

# Source utility functions
source("runtime-comparison/util.R")

# Constants
lambda1 <- 1
lambda2 <- 1

set.seed(1)

################################################################################
# Random binary graph 
################################################################################

#' This tibble contains all the parameter settings. It includes combinations of 
#' 'p' (size of the square weight matrix) and 'density' (density of the binary values).
parameter_settings <- tibble(expand.grid(
  p = seq(5,100,by=1), 
  density = c(.5)#seq(.1, 1, by = .1)
))

# number of parameter settings
n_parameter_settings <- nrow(parameter_settings)

# Create a progress bar
pb <- progress::progress_bar$new(total = n_parameter_settings)

#' Go over each of the parameter settings and check the run time
res <- lapply(1:n_parameter_settings, function(i) {
  
  # Get the parameters
  p <- as.integer(parameter_settings[i, 'p'])
  density <- as.double(parameter_settings[i,'density'])
  
  # Create the random binary matrix weight matrix
  W <- random_binary_weight_matrix(p, density)
  
  # Create the corresponding connListObj needed by the flsa package
  connListObj <- create_connListObj(W)
  
  # generate the raw data 
  y <- rnorm(p)
  
  # Applies the flsa function to the data
  FLSA <- function() {
    
    # if there are no connections what so ever
    if (connListObj$connList_is_null) {
      # lambda2 is zero, since there is no smoothness penalty
      flsa::flsa(y, lambda1 = lambda1, lambda2 = 0) 
    } else {
      flsa::flsa(y, lambda1 = lambda1, lambda2 = lambda2,
               connListObj = connListObj$connList)
    }
  }

  # Applies the wflsa package to the data
  wFLSA <- function() {
    wflsa::wflsa(y, W, lambda1 = lambda1, lambda2 = lambda2)
  }
  
  # Assesses the runtime 
  time <- microbenchmark::microbenchmark(wFLSA(), 
                                         FLSA(), times = 100)
 
  pb$tick()
  
  return(time)
})

# save the raw runtime 
saveRDS(res, "runtime-comparison/runtime-random-weight-matrix.rds")

################################################################################
# Classic 1-D FLSA
################################################################################

#' This tibble contains all the parameter settings
parameter_settings_1D <- tibble(expand.grid(
  p = seq(5,100, by = 1)
))

# number of parameter settings
n_parameter_settings <- nrow(parameter_settings_1D)

# Create a progress bar
pb <- progress::progress_bar$new(total = n_parameter_settings)

#' Go over each of the parameter settings and check the run time
res <- lapply(1:n_parameter_settings, function(i) {
  
  # Get the parameters
  p <- as.integer(parameter_settings_1D[i, 'p'])
  
  # Create the random binary matrix weight matrix
  W <- band_matrix(p)
  
  # generate the raw data 
  y <- rnorm(p)
  
  # Applies the flsa function to the data
  FLSA <- function() {
      flsa::flsa(y, lambda1 = lambda1, lambda2 = lambda2)
  }
  
  # Applies the wflsa package to the data
  wFLSA <- function() {
    wflsa::wflsa(y, W, lambda1 = lambda1, lambda2 = lambda2)
  }
  
  # Assesses the runtime 
  time <- microbenchmark::microbenchmark(wFLSA(), 
                                         FLSA(), times = 100)
  
  pb$tick()
  
  return(time)
})

# save the raw runtime 
saveRDS(res, "runtime-comparison/runtime-1D.rds")

################################################################################
# Process results
################################################################################

#' Create Plot for Runtime Assessment Results
#'
#' \code{create_plot} function reads in the runtime assessment results from a file, 
#' processes the data, and creates a plot to visualize the performance of different 
#' algorithms based on the parameter settings.
#'
#' @param filename The filename of the file containing the runtime assessment results.
#' @param parameter_settings A tibble containing the parameter settings used for the assessment.
#' @param title A title for the plot.
#' @return A plot visualizing the performance of different algorithms based on the parameter settings.
create_plot <- function(filename, 
                        parameter_settings = parameter_settings, 
                        title = "") {
  
  # read in the data
  res <- readRDS(filename)
  
  # go over all the results for the different parameter settings
  res <- lapply(1:length(res), function(i) {
    
    df <- as.data.frame(res[[i]])
    df$expr <- as.character(df$expr)
    df$expr <- substring(df$expr, 1, nchar(df$expr) - 2)
    colnames(df) <- c("algorithm", "time")
    
    df %>% mutate(p = as.integer(parameter_settings[i, 1]))
                  #density = as.double(parameter_settings[i,2]))
  })
  
  # combine all the results into a single tibble
  res <- do.call(rbind, res)
  
  # turn the time data into seconds
  res <- res %>% mutate(time = time / 1e9)
  
  # create the plot
  ggplot2::ggplot(res, aes(x = p, y = time, color = algorithm)) + 
    geom_smooth() +
    ylab("time (s)") + 
    ggtitle(title) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() + 
    guides(color=guide_legend(override.aes=list(fill=NA)))
}

p_random_graph <- create_plot("runtime-comparison/runtime-random-weight-matrix.rds", 
            parameter_settings = parameter_settings, 
            title = "Random Graph")

p_1d <- create_plot(filename = "runtime-comparison/runtime-1D.rds", 
                              parameter_settings = parameter_settings_1D, 
                              title = "1-D FLSA")

saveRDS(list(
    p_random_graph = p_random_graph, 
    p_1d = p_1d
  ), 'runtime-comparison/figures.rds')
  
ggsave("runtime-comparison/runtime-random-weight-matrix.pdf", 
       plot = p_random_graph, 
       width = 5, height = 3) 

ggsave("runtime-comparison/runtime-1D.pdf", 
       plot = p_1d, 
       width = 5, height = 3) 
