
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

# Generate test data
set.seed(123)
n <- 100
p <- 10
X <- matrix(rnorm(n * p), nrow = n)
Y <- rnorm(n)
std_result <- standardizeXY(X, Y)
Xtilde <- std_result$Xtilde
Ytilde <- std_result$Ytilde

# Test 1: Random beta
test_lasso_1 <- function() {
  beta <- rnorm(p)
  lambda <- 0.5
  expected <- lasso(Xtilde, Ytilde, beta, lambda)
  actual <- lasso_c(Xtilde, Ytilde, beta, lambda)
  cat("Test 1 (random beta): Expected =", expected, "Actual =", actual, "\n")
  all.equal(expected, actual, tolerance = 1e-10)
}

# Test 2: Zero beta
test_lasso_2 <- function() {
  beta <- rep(0, p)
  lambda <- 1.0
  expected <- lasso(Xtilde, Ytilde, beta, lambda)
  actual <- lasso_c(Xtilde, Ytilde, beta, lambda)
  cat("Test 2 (zero beta): Expected =", expected, "Actual =", actual, "\n")
  all.equal(expected, actual, tolerance = 1e-10)
}

# Run tests
test_lasso_1()
test_lasso_2()


# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################

# Test 1: Small lambda
test_fitLASSO_1 <- function() {
  lambda <- 0.1
  result_r <- fitLASSOstandardized(Xtilde, Ytilde, lambda)
  beta_r <- as.numeric(result_r$beta)
  
  beta_cpp <- as.numeric(fitLASSOstandardized_c(Xtilde, Ytilde, lambda, numeric(0), eps = 0.001))
  
  cat("Test 1 (lambda = 0.1): Max difference =", max(abs(beta_r - beta_cpp)), "\n")
  all.equal(beta_r, beta_cpp, tolerance = 1e-6, check.attributes = FALSE)
}

# Test 2: Large lambda
test_fitLASSO_2 <- function() {
  lambda <- 2.0
  result_r <- fitLASSOstandardized(Xtilde, Ytilde, lambda)
  beta_r <- as.numeric(result_r$beta)
  
  beta_cpp <- as.numeric(fitLASSOstandardized_c(Xtilde, Ytilde, lambda, numeric(0), eps = 0.001))
  
  cat("Test 2 (lambda = 2.0): Max difference =", max(abs(beta_r - beta_cpp)), "\n")
  all.equal(beta_r, beta_cpp, tolerance = 1e-6, check.attributes = FALSE)
}

# Run tests
test_fitLASSO_1()
test_fitLASSO_2()


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