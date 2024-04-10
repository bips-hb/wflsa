pkgname <- "wflsa"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "wflsa-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('wflsa')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("get_estimate_lambda1")
### * get_estimate_lambda1

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_estimate_lambda1
### Title: Get Estimates for Beta Coefficients for a Given Lambda1 Value
### Aliases: get_estimate_lambda1

### ** Examples

# Example usage:
# fit <- wflsa(data, lambda1_values, lambda2_values)
# estimate <- get_estimate_lambda1(fit, lambda1_value)
# Now you can access estimated beta coefficients for different lambda2 values



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_estimate_lambda1", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_unique_values")
### * get_unique_values

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_unique_values
### Title: Get unique absolute values of a vector after rounding
### Aliases: get_unique_values

### ** Examples

vector <- c(1.23456789, -1.23456788, 1.23456787, -1.23456786, 1.23456785, -1.23456784)
sorted_unique_abs_vals <- get_sorted_unique_abs_values(vector)
print(sorted_unique_abs_vals)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_unique_values", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("soft_threshold")
### * soft_threshold

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: soft_threshold
### Title: Apply soft thresholding to a numeric vector
### Aliases: soft_threshold

### ** Examples

x <- c(3, -2, 5, -4, 1)
lambda1 <- 2
soft_thresholded <- soft_threshold(x, lambda1)
print(soft_thresholded)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("soft_threshold", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("wflsa")
### * wflsa

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: wflsa
### Title: The Weighted Fused Lasso Signal Approximator (wFLSA) Algorithm
### Aliases: wflsa

### ** Examples

# Example usage of the wflsa function
y <- c(1, 2, 3)
W <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
lambda2 <- c(0.1, 0.2)
result <- wflsa::wflsa(y, W, lambda2)
wflsa::get_estimate_lambda1(result, lambda1 = .5)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("wflsa", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
