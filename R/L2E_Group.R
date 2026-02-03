
#------------------------------------------------------------------------------
# 1. Beta update in L2E group-lasso regression
#------------------------------------------------------------------------------
#' Beta update in L2E group-lasso regression
#'
#' \code{update_beta_group} updates beta for robust L2E regression using a
#' group lasso penalty via the grpreg package.
#'
#' @param y Response vector
#' @param X Design matrix
#' @param beta Initial vector of regression coefficients
#' @param tau Initial precision estimate
#' @param lambda Tuning parameter
#' @param penalty The penalty to be applied to the model. For group selection, one of grLasso, grMCP, or grSCAD
#' @param group Group indicator vector (assigning each predictor to a group)
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @return Returns a list object containing the updated beta (vector) and the number of iterations used.
#' @importFrom Matrix Diagonal
#' @importFrom grpreg grpreg
#' 
#' @keywords internal
update_beta_group <- function(y, X, beta, tau, lambda, penalty, group, max_iter = 1e2, tol = 1e-4) {
  n <- nrow(X)
  
  for (i in 1:max_iter) {
    beta_last <- beta
    Xbeta <- X %*% beta
    r <- y - Xbeta
    # Compute weights (same as in the sparse update)
    w <- as.vector(exp(-0.5 * (tau * r)^2))
    
    # Form the weighted design matrix and response
    W <- Diagonal(n = n, x = sqrt(w))
    Xtilde <- as.matrix(W %*% X)
    ytilde <- as.vector(W %*% y)
    
    # grpreg
    fit <- grpreg(X = Xtilde, y = ytilde, group = group, penalty = penalty, family = "gaussian",
                  lambda = lambda)
    beta <- as.vector(fit$beta)[-1]
    
    # Check convergence on beta
    if (norm(as.matrix(beta_last - beta), 'f') < tol * (1 + norm(as.matrix(beta_last), 'f'))) break
  }
  
  return(list(beta = beta, iter = i))
}




#------------------------------------------------------------------------------
# 2. L2E group-lasso regression: update both beta and tau with Nesterov's Acceleration
#------------------------------------------------------------------------------
#' L2E group-lasso regression with
#'
#' \code{l2e_regression_group} performs robust L2E regression using a group lasso penalty.
#'
#' @param y Response vector
#' @param X Design matrix
#' @param beta Initial vector of regression coefficients
#' @param tau Initial precision estimate
#' @param lambda Tuning parameter
#' @param group Group indicator vector
#' @param penalty The penalty to be applied to the model. For group selection, one of grLasso, grMCP, or grSCAD
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Logical; if TRUE, prints the computing time
#' @return Returns a list containing the updated beta, tau, and number of iterations.
#' 
#' @keywords internal
l2e_regression_group_accel <- function(y, X, beta, tau, lambda, group, penalty,
                                       max_iter = 1e2, tol = 1e-4, Show.Time = TRUE) {
  if (tau <= 0) stop("Entered non-positive initial tau")
  
  # initialize FISTA momentum parameters
  t_curr <- 1
  beta_prev <- beta
  beta_curr <- beta
  
  start_time <- proc.time()
  for (i in 1:max_iter) {
    # compute new momentum weight 
    t_next <- (1 + sqrt(1 + 4 * t_curr^2)) / 2
    momentum <- (t_curr - 1) / t_next
    
    # extrapolate beta
    beta_extrap <- beta_curr + momentum * (beta_curr - beta_prev)
    
    # save for next iteration
    beta_prev <- beta_curr
    t_curr <- t_next
    
    # --- MM step 1: update beta (starting from extrapolated guess) ---
    sol_beta <- update_beta_group(y, X, beta_extrap, tau, lambda, penalty, group,
                                  max_iter = 1, tol = tol)
    beta_new <- sol_beta$beta
    
    
    # Update tau using update_eta_bktk function
    r <- y - X %*% beta_new
    eta_last <- log(tau)
    res_eta <- update_eta_bktk(r, eta_last, tol = tol)
    eta <- res_eta$eta
    tau <- max(exp(eta), 1e-10)
    
    # check convergence on both beta and tau
    conv_beta <- norm(as.matrix(beta_curr - beta_new), 'f') < tol * (1 + norm(as.matrix(beta_curr), 'f'))
    conv_tau  <- abs(eta_last - eta)              < tol * (1 + abs(eta_last))
    beta_curr <- beta_new
    if (conv_beta && conv_tau) break
  }
  
  if (Show.Time) print(proc.time() - start_time)
  return(list(beta = beta_curr, tau = tau, iter = i))
}




#------------------------------------------------------------------------------
# 3. Compute the solution path over a sequence of lambda values
#------------------------------------------------------------------------------
#' L2E group-lasso solution path with
#'
#' \code{L2E_group} computes the solution path for robust L2E regression with a group lasso penalty.
#'
#' @param y Response vector
#' @param X Design matrix
#' @param group Group indicator vector
#' @param b Initial vector of regression coefficients (optional)
#' @param tau Initial precision estimate (optional)
#' @param lambdaSeq A decreasing sequence of tuning parameter lambda values (optional)
#' @param penalty The penalty to be applied to the model. For group selection, one of grLasso, grMCP, or grSCAD
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Logical; if TRUE, prints the computing time
#' @return Returns a list with estimated Beta (matrix), Tau (vector), run time, and lambdaSeq.
#' 
#' @importFrom stats mad
#' @export
L2E_group <- function(y, X, group, b, tau, lambdaSeq, penalty = 'grLasso',
                          max_iter = 1e2, tol = 1e-4, Show.Time = TRUE) {
  if (missing(b)) {
    b <- double(ncol(X))
  }
  if (missing(tau)) {
    tau <- 1 / mad(y)
  }
  if (tau <= 0) stop("Entered non-positive initial tau")
  if (missing(lambdaSeq)) {
    lambdaSeq <- 10^seq(1, -4, length.out = 20)
  }
  
  Nlambda <- length(lambdaSeq)
  Beta <- matrix(0, nrow = ncol(X), ncol = Nlambda)
  Tau <- double(Nlambda)
  
  start_time <- proc.time()
  for (j in 1:Nlambda) {
    res <- l2e_regression_group_accel(y, X, b, tau, lambda = lambdaSeq[j], group = group, penalty = penalty,
                                    max_iter = max_iter, tol = tol, Show.Time = FALSE)
    Beta[, j] <- b <- res$beta
    Tau[j] <- tau <- res$tau
  }
  runtime <- proc.time() - start_time
  if (Show.Time) print(runtime)
  
  return(list(Beta = Beta, Tau = Tau, runtime = runtime, lambdaSeq = lambdaSeq))
}


#------------------------------------------------------------------------------
# 4. One CV fold for group-lasso L2E regression
#------------------------------------------------------------------------------
#' CV fold for L2E group-lasso regression
#'
#' \code{cv_fold_l2e_group} computes the loss for one cross-validation fold.
#'
#' @param i Fold index
#' @param y Response vector
#' @param X Design matrix
#' @param fold Vector indicating fold membership for each observation
#' @param cv.args List of arguments for L2E_group
#' @param method Character string ("median" or "mean") to compute the objective loss
#' @return A numeric vector of loss values for each lambda in lambdaSeq.
#' 
#' @keywords internal
cv_fold_l2e_group <- function(i, y, X, fold, cv.args, method = "median") {
  # Training data (all folds except fold i)
  cv.args$y <- y[fold != i]
  cv.args$X <- X[fold != i, , drop = FALSE]
  fit_i <- do.call("L2E_group", cv.args)
  
  # Hold-out data (fold i)
  y_out <- y[fold == i]
  X_out <- X[fold == i, , drop = FALSE]
  
  L <- length(fit_i$lambdaSeq)
  loss <- double(L)
  
  for (l in 1:L) {
    bhat <- fit_i$Beta[, l]
    Xbeta <- X_out %*% bhat
    r <- y_out - Xbeta
    tauhat <- fit_i$Tau[l]
    
    # Compute the objective loss using objective_tau
    loss[l] <- objective_tau(tau = tauhat, r = r, method = method)
  }
  return(loss)
}

#------------------------------------------------------------------------------
# 5. Cross-validation wrapper for group-lasso L2E regression
#------------------------------------------------------------------------------
#' Cross-validation for L2E group-lasso regression
#'
#' \code{CV_L2E_group} performs k-fold cross-validation for robust L2E regression
#' with a group lasso penalty.
#'
#' @param y Response vector
#' @param X Design matrix
#' @param group Group indicator vector
#' @param b Initial vector of regression coefficients (optional)
#' @param tau Initial precision estimate (optional)
#' @param lambdaSeq A decreasing sequence of tuning parameter lambda values (optional)
#' @param nfolds The number of cross-validation folds (default is 5)
#' @param seed Seed for reproducibility (optional)
#' @param penalty The penalty to be applied to the model. For group selection, one of grLasso, grMCP, or grSCAD
#' @param method Character string ("median" or "mean") to compute the objective loss
#' @param max_iter Maximum number of iterations for the regression update
#' @param tol Relative tolerance
#' @param trace Logical; if TRUE, prints progress for each fold
#' @return Returns a list containing the CV error (cve), standard error (cvse),
#' the index and value of the lambda with minimum CV error, the 1-standard-error lambda,
#' the sequence of lambda values, and the fold assignments.
#' 
#' @importFrom stats mad
#' @importFrom stats sd
#' @export
#' 
CV_L2E_group <- function(y, X, group, b, tau, lambdaSeq, nfolds = 5, seed = 1234, penalty = "grLasso",
                             method = "median", max_iter = 1e2, tol = 1e-4, trace = TRUE) {
  if (missing(lambdaSeq)) {
    lambdaSeq <- 10^seq(1, -4, length.out = 20)
  }
  if (missing(b)) {
    b <- double(ncol(X))
  }
  if (missing(tau)) {
    tau <- 1 / mad(y)
  }
  if (tau <= 0) stop("Entered non-positive tau")
  
  set.seed(seed)
  n <- length(y)
  fold <- sample(1:n %% nfolds)
  fold[fold == 0] <- nfolds
  
  # Set common arguments for L2E_group
  cv.args <- list()
  cv.args$b <- b
  cv.args$tau <- tau
  cv.args$lambdaSeq <- lambdaSeq
  cv.args$max_iter <- max_iter
  cv.args$tol <- tol
  cv.args$Show.Time <- FALSE
  cv.args$group <- group
  cv.args$penalty <- penalty
  
  Loss <- matrix(0, nrow = nfolds, ncol = length(lambdaSeq))
  for (i in 1:nfolds) {
    if (trace) cat("Starting CV fold #", i, "\n")
    res <- cv_fold_l2e_group(i, y, X, fold, cv.args, method = method)
    Loss[i, ] <- res
  }
  
  # Compute the mean and standard error of the loss across folds
  cve <- apply(Loss, 2, mean)
  cvse <- apply(Loss, 2, sd) / sqrt(nfolds)
  min_index <- which.min(round(cve, 8))
  
  # Determine the lambda using the one-standard-error rule
  for (i in min_index:1) {
    if (cve[i] > cve[min_index] + cvse[min_index])
      break
  }
  if (min_index == 1) {
    lambda.1se <- lambdaSeq[1]
    min_1se <- 1
  } else {
    lambda.1se <- lambdaSeq[i + 1]
    min_1se <- i + 1
  }
  
  return(list(cve = cve, cvse = cvse, min = min_index, lambda.min = lambdaSeq[min_index],
              min_1se = min_1se, lambda.1se = lambda.1se, lambdaSeq = lambdaSeq, fold = fold))
}





