# [ToDo] Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y){
  # [ToDo] Center Y
  Ymean <- mean(Y)
  Ytilde <- as.numeric(Y - Ymean)
  
  # [ToDo] Center and scale X
  Xmeans <- colMeans(X)
  Xcentered <- sweep(X, 2, Xmeans, "-")
  n <- nrow(Xcentered)
  # weights are sqrt((X_j^T X_j)/n) after centering
  weights <- sqrt(colSums(Xcentered * Xcentered) / n)
  # avoid division by zero for zero-variance columns
  safe_weights <- ifelse(weights == 0, 1, weights)
  Xtilde <- sweep(Xcentered, 2, safe_weights, "/")
  if (any(weights == 0)) {
    Xtilde[, weights == 0] <- 0
  }
  
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# [ToDo] Soft-thresholding of a scalar a at level lambda 
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda){
  return(sign(a) * pmax(abs(a) - lambda, 0))
}

# [ToDo] Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda){
  Xmat <- as.matrix(Xtilde)
  n <- nrow(Xmat)
  r <- as.numeric(Ytilde - Xmat %*% beta)
  data_term <- sum(r * r) / (2 * n)
  pen_term <- lambda * sum(abs(beta))
  return(data_term + pen_term)
}

# [ToDo] Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.001){
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  Xmat <- as.matrix(Xtilde)
  n <- nrow(Xmat)
  p <- ncol(Xmat)
  if (length(Ytilde) != n) stop("Xtilde and Ytilde must have matching n")
  
  #[ToDo]  Check that lambda is non-negative
  if (lambda < 0) stop("lambda must be non-negative")
  
  #[ToDo]  Check for starting point beta_start. 
  # If none supplied, initialize with a vector of zeros.
  # If supplied, check for compatibility with Xtilde in terms of p
  if (is.null(beta_start)) {
    beta <- rep(0, p)
  } else {
    if (length(beta_start) != p) stop("beta_start length must match number of columns in Xtilde")
    beta <- as.numeric(beta_start)
  }
  
  #[ToDo]  Coordinate-descent implementation. 
  # Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999,
  # your should return fmin = 0.99999, and not have another iteration
  
  # Maintain residual for efficiency: r = Y - X %*% beta
  r <- as.numeric(Ytilde - Xmat %*% beta)
  f_prev <- Inf
  repeat {
    # One full cyclic sweep over coordinates
    for (j in seq_len(p)) {
      xj <- Xmat[, j]
      # Add back current contribution of feature j to the residual
      r <- r + xj * beta[j]
      # With standardized columns (n^{-1} xj^T xj = 1), the update is soft(mean(xj * r), lambda)
      rho <- sum(xj * r) / n
      beta_new_j <- soft(rho, lambda)
      # Update residual and coefficient
      r <- r - xj * beta_new_j
      beta[j] <- beta_new_j
    }
    f_curr <- lasso(Xmat, Ytilde, beta, lambda)
    if ((f_prev - f_curr) < eps) {
      fmin <- f_curr
      break
    }
    f_prev <- f_curr
  }
  
  # Return 
  # beta - the solution (a vector)
  # fmin - optimal function value (value of objective at beta, scalar)
  return(list(beta = beta, fmin = fmin))
}

# [ToDo] Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence,
#             is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Check that n is the same between Xtilde and Ytilde
  Xmat <- as.matrix(Xtilde)
  n <- nrow(Xmat)
  p <- ncol(Xmat)
  if (length(Ytilde) != n) stop("Xtilde and Ytilde must have matching n")
 
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >= 0,
  # and make sure the values are sorted from largest to smallest.
  # If none of the supplied values satisfy the requirement,
  # print the warning message and proceed as if the values were not supplied.
  if (!is.null(lambda_seq)) {
    lambda_seq <- lambda_seq[lambda_seq >= 0]
    if (length(lambda_seq) == 0) {
      warning("No valid lambda values supplied, generating default sequence")
      lambda_seq <- NULL
    } else {
      lambda_seq <- sort(lambda_seq, decreasing = TRUE)
    }
  }
  
  # If lambda_seq is not supplied, calculate lambda_max 
  # (the minimal value of lambda that gives zero solution),
  # and create a sequence of length n_lambda as
  if (is.null(lambda_seq)) {
    # lambda_max = max_j |X_j^T Y| / n for standardized X
    lambda_max <- max(abs(t(Xmat) %*% Ytilde)) / n
    lambda_seq <- exp(seq(log(lambda_max), log(0.01), length = n_lambda))
  }
  
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda 
  # (make sure supplied eps is carried over). 
  # Use warm starts strategy discussed in class for setting the starting values.
  n_lambda_actual <- length(lambda_seq)
  beta_mat <- matrix(0, nrow = p, ncol = n_lambda_actual)
  fmin_vec <- numeric(n_lambda_actual)
  
  beta_start <- NULL
  for (i in seq_len(n_lambda_actual)) {
    result <- fitLASSOstandardized(Xmat, Ytilde, lambda_seq[i], beta_start, eps)
    beta_mat[, i] <- result$beta
    fmin_vec[i] <- result$fmin
    # Warm start: use solution from current lambda as starting point for next
    beta_start <- result$beta
  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

# [ToDo] Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  std_result <- standardizeXY(X, Y)
  Xtilde <- std_result$Xtilde
  Ytilde <- std_result$Ytilde
  Ymean <- std_result$Ymean
  Xmeans <- std_result$Xmeans
  weights <- std_result$weights
 
  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq
  # (make sure the parameters carry over)
  fit_result <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq, n_lambda, eps)
  lambda_seq <- fit_result$lambda_seq
  beta_tilde_mat <- fit_result$beta_mat
  fmin_vec <- fit_result$fmin_vec
 
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector
  # for each lambda
  p <- ncol(X)
  n_lambda_actual <- length(lambda_seq)
  beta_mat <- matrix(0, nrow = p, ncol = n_lambda_actual)
  beta0_vec <- numeric(n_lambda_actual)
  
  for (i in seq_len(n_lambda_actual)) {
    # Back-scale coefficients: beta_original = beta_tilde / weights
    beta_mat[, i] <- beta_tilde_mat[, i] / weights
    # Back-center intercept: beta0 = Ymean - Xmeans^T * beta_original
    beta0_vec[i] <- Ymean - sum(Xmeans * beta_mat[, i])
  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


# [ToDo] Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored 
# eps - precision level for convergence assessment, default 0.001
cvLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, k = 5, fold_ids = NULL, eps = 0.001){
  # [ToDo] Fit Lasso on original data using fitLASSO
  full_fit <- fitLASSO(X, Y, lambda_seq, n_lambda, eps)
  lambda_seq <- full_fit$lambda_seq
  beta_mat <- full_fit$beta_mat
  beta0_vec <- full_fit$beta0_vec
  n_lambda_actual <- length(lambda_seq)
 
  # [ToDo] If fold_ids is NULL, split the data randomly into k folds.
  # If fold_ids is not NULL, split the data according to supplied fold_ids.
  n <- nrow(X)
  if (is.null(fold_ids)) {
    # Random assignment to k folds
    fold_ids <- sample(rep(seq_len(k), length.out = n))
  } else {
    # Use supplied fold_ids and update k
    k <- max(fold_ids)
  }
  
  # [ToDo] Calculate LASSO on each fold using fitLASSO,
  # and perform any additional calculations needed for CV(lambda) and SE_CV(lambda)
  cv_errors <- matrix(0, nrow = k, ncol = n_lambda_actual)
  
  for (fold in seq_len(k)) {
    # Split data into training and validation sets
    train_idx <- which(fold_ids != fold)
    val_idx <- which(fold_ids == fold)
    
    X_train <- X[train_idx, , drop = FALSE]
    Y_train <- Y[train_idx]
    X_val <- X[val_idx, , drop = FALSE]
    Y_val <- Y[val_idx]
    
    # Fit LASSO on training data
    fold_fit <- fitLASSO(X_train, Y_train, lambda_seq, n_lambda_actual, eps)
    
    # Calculate prediction errors on validation set for each lambda
    for (i in seq_len(n_lambda_actual)) {
      y_pred <- fold_fit$beta0_vec[i] + X_val %*% fold_fit$beta_mat[, i]
      cv_errors[fold, i] <- mean((Y_val - y_pred)^2)
    }
  }
  
  # Calculate CV mean and standard error for each lambda
  cvm <- colMeans(cv_errors)
  cvse <- apply(cv_errors, 2, sd) / sqrt(k)
  
  # [ToDo] Find lambda_min
  min_idx <- which.min(cvm)
  lambda_min <- lambda_seq[min_idx]

  # [ToDo] Find lambda_1SE
  # Find the largest lambda within 1 SE of the minimum
  min_cv <- cvm[min_idx]
  min_cv_se <- cvse[min_idx]
  threshold <- min_cv + min_cv_se
  
  # Find all lambdas with CV error <= threshold
  valid_idx <- which(cvm <= threshold)
  if (length(valid_idx) > 0) {
    # Take the largest lambda (smallest index since lambda_seq is decreasing)
    lambda_1se_idx <- min(valid_idx)
    lambda_1se <- lambda_seq[lambda_1se_idx]
  } else {
    # Fallback to lambda_min if no lambda satisfies 1SE rule
    lambda_1se <- lambda_min
  }
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec, fold_ids = fold_ids, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse))
}

