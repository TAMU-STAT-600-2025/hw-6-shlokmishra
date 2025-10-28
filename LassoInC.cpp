#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Soft-thresholding function, returns scalar
// [[Rcpp::export]]
double soft_c(double a, double lambda){
  return std::copysign(std::max(std::abs(a) - lambda, 0.0), a);
}

// Lasso objective function, returns scalar
// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  arma::colvec r = Ytilde - Xtilde * beta;
  double n = Xtilde.n_rows;
  double data_term = arma::accu(r % r) / (2.0 * n);
  double pen_term = lambda * arma::accu(arma::abs(beta));
  return data_term + pen_term;
}

// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.001){
  int n = Xtilde.n_rows;
  int p = Xtilde.n_cols;
  
  // Initialize beta
  arma::colvec beta;
  if (beta_start.n_elem == p) {
    beta = beta_start;
  } else {
    beta = arma::zeros(p);
  }
  
  // Initialize residual
  arma::colvec r = Ytilde - Xtilde * beta;
  
  double f_prev = std::numeric_limits<double>::infinity();
  double f_curr;
  
  repeat {
    // One full cyclic sweep over coordinates
    for (int j = 0; j < p; j++) {
      // Add back current contribution of feature j to the residual
      r += Xtilde.col(j) * beta(j);
      
      // With standardized columns, the update is soft(mean(xj * r), lambda)
      double rho = arma::accu(Xtilde.col(j) % r) / n;
      double beta_new_j = soft_c(rho, lambda);
      
      // Update residual and coefficient
      r -= Xtilde.col(j) * beta_new_j;
      beta(j) = beta_new_j;
    }
    
    // Check convergence
    f_curr = lasso_c(Xtilde, Ytilde, beta, lambda);
    if ((f_prev - f_curr) < eps) {
      break;
    }
    f_prev = f_curr;
  }
  
  return beta;
}  

// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.001){
  int p = Xtilde.n_cols;
  int n_lambda = lambda_seq.n_elem;
  
  // Initialize beta matrix
  arma::mat beta_mat(p, n_lambda);
  
  // Warm start strategy: use previous solution as starting point for next lambda
  arma::colvec beta_start;  // empty initially
  
  for (int i = 0; i < n_lambda; i++) {
    arma::colvec beta = fitLASSOstandardized_c(Xtilde, Ytilde, lambda_seq(i), beta_start, eps);
    beta_mat.col(i) = beta;
    beta_start = beta;  // Warm start for next iteration
  }
  
  return beta_mat;
}