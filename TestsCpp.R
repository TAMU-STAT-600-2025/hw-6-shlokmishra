
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

# Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Test 1: Positive value
test_soft_1 <- function() {
  a <- 5.0
  lambda <- 2.0
  expected <- soft(a, lambda)
  actual <- soft_c(a, lambda)
  cat("Test 1 (positive): Expected =", expected, "Actual =", actual, "\n")
  all.equal(expected, actual, tolerance = 1e-10)
}

# Test 2: Negative value
test_soft_2 <- function() {
  a <- -3.0
  lambda <- 1.0
  expected <- soft(a, lambda)
  actual <- soft_c(a, lambda)
  cat("Test 2 (negative): Expected =", expected, "Actual =", actual, "\n")
  all.equal(expected, actual, tolerance = 1e-10)
}

# Run tests
test_soft_1()
test_soft_2()


# Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
#################################################


# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################

# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################

# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Make sure riboflavin$x is treated as matrix later in the code for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Standardize the data
out <- standardizeXY(riboflavin$x, riboflavin$y)

# This is just to create lambda_seq, can be done faster, but this is simpler
outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)

# The code below should assess your speed improvement on riboflavin data
microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
  fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
  times = 10
)