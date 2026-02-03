#' Beta update in L2E sparse regression - NCV
#'
#' \code{update_beta_sparse_ncv} updates beta for L2E sparse regression using existing penalization methods
#'
#' @param y Response vector
#' @param X Design matrix
#' @param beta Initial vector of regression coefficients
#' @param tau Initial precision estimate
#' @param lambda Tuning parameter
#' @param penalty Available penalties include lasso, MCP and SCAD.
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @return Returns a list object containing the new estimate for beta (vector) and the number of iterations (scalar).
#' @importFrom Matrix Diagonal
#' @importFrom ncvreg ncvfit
#'
#' @keywords internal
#'
update_beta_sparse_ncv <- function(y,X,beta,tau,lambda, penalty, max_iter=1e2,tol=1e-4) {
  
  n <- nrow(X)
  
  for (i in 1:max_iter) {
    
    beta_last <- beta
    Xbeta <- X %*% beta
    r <- y - Xbeta
    w <- as.vector(exp(-0.5*(tau*r)^2 ))
    
    
    W <- Diagonal(n=n, x = sqrt(as.vector(w)))
    Xtilde <- as.matrix(W%*%X)
    ytilde <- as.vector(W%*%y)
    beta <- as.vector(ncvfit(Xtilde, ytilde, init = beta_last, penalty=penalty,lambda = lambda,
                             max.iter = 100, warn = FALSE)$beta)
    
    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
  }
  
  return(list(beta=beta,iter=i))
}



#' L2E sparse regression with existing penalization methods
#'
#' \code{l2e_regression_sparse_ncv} performs robust sparse regression under the L2 criterion. Available penalties include lasso, MCP and SCAD.

#' @param y Response vector
#' @param X Design matrix
#' @param beta Initial vector of regression coefficients
#' @param tau Initial precision estimate
#' @param lambda Tuning parameter
#' @param penalty Available penalties include lasso, MCP and SCAD.
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Report the computing time
#'
#' @keywords internal
#'
l2e_regression_sparse_ncv <- function(y, X, beta, tau, lambda, penalty, max_iter=1e2,
                                      tol=1e-4, Show.Time=TRUE) {
  
  if (tau <= 0) stop("Entered non-positive initial tau")
  
  # initialize FISTA momentum parameters
  t_curr <- 1
  beta_prev <- beta
  beta_curr <- beta
  
  time <- proc.time()
  for (i in 1:max_iter) {
    
    # compute new momentum weight 
    t_next <- (1 + sqrt(1 + 4 * t_curr^2)) / 2
    momentum <- (t_curr - 1) / t_next
    
    # extrapolate beta
    beta_extrap <- beta_curr + momentum * (beta_curr - beta_prev)
    
    # save for next iteration
    beta_prev <- beta_curr
    t_curr <- t_next
    
    # update beta
    sol_beta <- update_beta_sparse_ncv(y, X, beta_extrap, tau, lambda, penalty, max_iter=1e2,tol=1e-4)
    beta_new <- sol_beta$beta
    
    # update tau
    r <- y - X%*%beta_new
    eta_last <- log(tau)  # get eta as in line 9
    res_eta <- update_eta_bktk(r,eta_last, tol=tol) # update eta as in line 10-12
    eta <- res_eta$eta
    tau <- max(exp(eta), 1e-10) # update tau as in line 13
    
    
    # Check for convergence
    A <- norm(as.matrix(beta_curr-beta_new),'f') < tol*(1 + norm(as.matrix(beta_curr),'f'))
    B <- abs(eta_last-eta) < tol*(1 + abs(eta_last))
    beta_curr <- beta_new
    if (A & B) break
    
    
  }
  if(Show.Time) print(proc.time() - time)
  
  return(list(beta=beta_curr,tau=tau, iter=i ))
  
}




#' Solution path of L2E sparse regression with existing penalization methods
#'
#' \code{L2E_sparse_ncv} computes the solution path of robust sparse regression under the L2 criterion. Available penalties include lasso, MCP and SCAD.
#'
#' @param y Response vector
#' @param X Design matrix
#' @param b Initial vector of regression coefficients, can be omitted
#' @param tau Initial precision estimate, can be omitted
#' @param lambdaSeq A decreasing sequence of values for the tuning parameter lambda, can be omitted
#' @param penalty Available penalties include lasso, MCP and SCAD.
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Report the computing time
#' @param refit Logical; whether to refit the model on the active set.
#' @return Returns a list object containing the estimates for beta (matrix) and
#' tau (vector) for each value of the tuning parameter lambda,
#' the run time (vector) for each lambda,
#' and the sequence of lambda used in the regression (vector)
#' @importFrom stats mad
#' @export
#' @examples
#' set.seed(12345)
#' n <- 100
#' tau <- 1
#' f <- matrix(c(rep(2,5), rep(0,45)), ncol = 1)
#' X <- X0 <- matrix(rnorm(n*50), nrow = n)
#' y <- y0 <- X0 %*% f + (1/tau)*rnorm(n)
#'
#' ## Clean Data
#' lambda <- 10^(-1)
#' sol <- L2E_sparse_ncv(y=y, X=X, lambdaSeq=lambda, penalty="SCAD")
#' r <- y - X %*% sol$Beta
#' ix <- which(abs(r) > 3/sol$Tau)
#' l2e_fit <- X %*% sol$Beta
#'
#' plot(y, l2e_fit, ylab='Predicted values', pch=16, cex=0.8)
#' points(y[ix], l2e_fit[ix], pch=16, col='blue', cex=0.8)
#'
#' ## Contaminated Data
#' i <- 1:5
#' y[i] <- 2 + y0[i]
#' X[i,] <- 2 + X0[i,]
#'
#' sol <- L2E_sparse_ncv(y=y, X=X, lambdaSeq=lambda, penalty="SCAD")
#' r <- y - X %*% sol$Beta
#' ix <- which(abs(r) > 3/sol$Tau)
#' l2e_fit <- X %*% sol$Beta
#'
#' plot(y, l2e_fit, ylab='Predicted values', pch=16, cex=0.8)
#' points(y[ix], l2e_fit[ix], pch=16, col='blue', cex=0.8)
#'
L2E_sparse_ncv <- function(y,X,b,tau,lambdaSeq,penalty="MCP",max_iter=1e2,tol=1e-4,Show.Time=TRUE,refit=FALSE) {
  
  if(missing(b)){
    b <- double(ncol(X))  # initial beta
  }
  
  if(missing(tau)){
    tau <- 1/mad(y)   # initial tau
  }
  
  if (tau <= 0) stop("Entered non-positive initial tau")
  
  if(missing(lambdaSeq)){
    lambdaSeq <- 10^seq(1, -4, length.out = 20)  # set a sequence of rho
  }
  
  
  Nlambda <- length(lambdaSeq)
  
  # save outputs
  Beta <- matrix(0, nrow = ncol(X), ncol = Nlambda)
  Tau <- double(Nlambda)
  
  time <- proc.time()
  for (j in 1:Nlambda) {
    
    res <- l2e_regression_sparse_ncv(y, X, b, tau, lambda=lambdaSeq[j], penalty=penalty,
                                     max_iter=max_iter, tol=tol, Show.Time = FALSE)
    Beta[, j] <- b <- res$beta
    Tau[j] <- tau <- res$tau
    
    if(refit){
      ## Get active set to refine the raw L2E betahat
      activesetL2E <- which(res$beta != 0)
      
      # if all parameter are penalized to 0, then just return the raw values
      if(length(activesetL2E) > 0){
        
        X_refine <- X[,activesetL2E,drop=FALSE] # make new covariates with just active set
        
        # run pure L2E without penalty term (i.e. lambda = 0)
        res_refit <- l2e_regression_sparse_ncv(y, X_refine, rep(0, length(activesetL2E)), tau, lambda=0, penalty=penalty,
                                               max_iter=max_iter, tol=tol, Show.Time = FALSE)
        
        beta_refit <- numeric(length=length(res$beta))
        beta_refit[activesetL2E] <- res_refit$beta
        
        Beta[, j] <- b <- beta_refit
        Tau[j] <- tau <- res_refit$tau
        
      }
    }
    
  }
  runtime <-  proc.time() - time
  if(Show.Time) print(runtime)
  
  return(list(Beta=Beta, Tau=Tau, runtime=runtime, lambdaSeq=lambdaSeq))
  
}


#' Internal function for k-fold cross-validation
#'
#' \code{cv_fold_l2e_ncv} fits the L2E sparse regression model on a training fold 
#' and evaluates the robust loss on the hold-out fold.
#'
#' @param i Integer; the index of the current fold to be used as the validation set.
#' @param y Response vector for the full dataset.
#' @param X Design matrix for the full dataset.
#' @param fold A vector of integers indicating which fold each observation belongs to.
#' @param cv.args A list of arguments to be passed to \code{L2E_sparse_ncv}.
#' @param method Character; the method to compute the objective loss, 
#' either "mean" or "median" (default). "median" is recommended for robustness.
#'
#' @return A numeric vector of the computed loss for each value in the 
#' \code{lambdaSeq} provided in \code{cv.args}.
#'
#' @keywords internal
cv_fold_l2e_ncv <- function(i, y, X, fold, cv.args, method="median") {
  cv.args$y <- y[fold!=i]
  cv.args$X <- X[fold!=i, , drop=FALSE]
  fit.i <- do.call("L2E_sparse_ncv", cv.args)
  
  # data in hold-out
  y_out <- y[fold==i]
  X_out <- X[fold==i, , drop=FALSE]
  
  
  L <- length(fit.i$lambdaSeq)
  loss <- double(L)
  
  for (l in 1:L) {
    bhat <- fit.i$Beta[, l]  # get the estimated beta with the l-th k
    Xbeta <- X_out %*% bhat
    r <- y_out - Xbeta
    tauhat <- fit.i$Tau[l]
    
    loss[l] <- objective_tau(tau = tauhat, r = r, method=method) ### use median instead of mean to account for outliers
  }
  
  return(loss)
}


#' Cross validation for L2E sparse regression with existing penalization methods
#'
#' \code{CV_L2E_sparse_ncv} performs k-fold cross-validation for robust sparse regression under the L2 criterion.
#'  Available penalties include lasso, MCP and SCAD.
#'
#' @param y Response vector
#' @param X Design matrix
#' @param beta0 Initial vector of regression coefficients, can be omitted
#' @param tau0 Initial precision estimate, can be omitted
#' @param lambdaSeq A decreasing sequence of tuning parameter lambda, can be omitted
#' @param penalty Available penalties include lasso, MCP and SCAD.
#' @param nfolds The number of cross-validation folds. Default is 5.
#' @param seed Users can set the seed of the random number generator to obtain reproducible results.
#' @param method Median or mean to compute the objective
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param trace Whether to trace the progress of the cross-validation
#' @return Returns a list object containing the mean and standard error of the cross-validation error -- CVE and CVSE -- for each value of k (vectors),
#' the index of the lambda with the minimum CVE and the lambda value itself (scalars),
#' the index of the lambda value with the 1SE CVE and the lambda value itself (scalars),
#' the sequence of lambda used in the regression (vector), and
#' a vector listing which fold each element of y was assigned to
#' 
#' @importFrom stats mad
#' @importFrom stats sd
#' @export
#' @examples
#' ## Completes in 20 seconds
#'
#' set.seed(12345)
#' n <- 100
#' tau <- 1
#' f <- matrix(c(rep(2,5), rep(0,45)), ncol = 1)
#' X <- X0 <- matrix(rnorm(n*50), nrow = n)
#' y <- y0 <- X0 %*% f + (1/tau)*rnorm(n)
#'
#' ## Clean Data
#' lambda <- 10^seq(-1, -2, length.out=20)
#' # (not run)
#' # cv <- CV_L2E_sparse_ncv(y=y, X=X, lambdaSeq=lambda, penalty="SCAD", seed=1234, nfolds=2)
#' # (lambda_min <- cv$lambda.min)
#'
#' # sol <- L2E_sparse_ncv(y=y, X=X, lambdaSeq=lambda_min, penalty="SCAD")
#' # r <- y - X %*% sol$Beta
#' # ix <- which(abs(r) > 3/sol$Tau)
#' # l2e_fit <- X %*% sol$Beta
#'
#' # plot(y, l2e_fit, ylab='Predicted values', pch=16, cex=0.8)
#' # points(y[ix], l2e_fit[ix], pch=16, col='blue', cex=0.8)
#'
#' ## Contaminated Data
#' i <- 1:5
#' y[i] <- 2 + y0[i]
#' X[i,] <- 2 + X0[i,]
#'
#' # (not run)
#' # cv <- CV_L2E_sparse_ncv(y=y, X=X, lambdaSeq=lambda, penalty="SCAD", seed=1234, nfolds=2)
#' # (lambda_min <- cv$lambda.min)
#'
#' # sol <- L2E_sparse_ncv(y=y, X=X, lambdaSeq=lambda_min, penalty="SCAD")
#' # r <- y - X %*% sol$Beta
#' # ix <- which(abs(r) > 3/sol$Tau)
#' # l2e_fit <- X %*% sol$Beta
#'
#' # plot(y, l2e_fit, ylab='Predicted values', pch=16, cex=0.8)
#' # points(y[ix], l2e_fit[ix], pch=16, col='blue', cex=0.8)
#'
CV_L2E_sparse_ncv <- function(y, X, beta0, tau0, lambdaSeq,  penalty="MCP", nfolds=5, seed=1234, method="median",
                              max_iter=1e2, tol=1e-4, trace=TRUE) {
  
  
  if(missing(lambdaSeq)){
    lambdaSeq <- 10^seq(1, -4, length.out = 20)  # set a sequence of lambda
  }
  
  if(missing(beta0)){
    beta0 <- double(ncol(X))  # initial beta
  }
  
  if(missing(tau0)){
    tau0 <- 1/mad(y)  # initial tau
  }
  
  if (tau0 <= 0) stop("Entered non-positive tau0")
  
  
  # Set up folds
  if (!missing(seed)) set.seed(seed)
  n <- length(y)
  fold <- sample(1:n %% nfolds)
  fold[fold==0] <- nfolds
  
  
  # Do cross-validation
  
  cv.args <- list()
  cv.args$b <- beta0
  cv.args$tau <- tau0
  cv.args$lambdaSeq <- lambdaSeq
  cv.args$penalty <- penalty
  cv.args$max_iter <- max_iter
  cv.args$tol <- tol
  cv.args$Show.Time <- FALSE
  cv.args$refit <- FALSE
  
  Loss <- matrix(0, nrow = nfolds, ncol = length(lambdaSeq))
  for (i in 1:nfolds) {
    if (trace) cat("Starting CV fold #", i, sep="","\n")
    res <- cv_fold_l2e_ncv(i, y, X, fold, cv.args, method=method)
    Loss[i, ] <- res
  }
  
  # Return
  cve <- apply(Loss, 2, mean)
  cvse <- apply(Loss, 2, sd)/sqrt(nfolds)
  min <- which.min(round(cve,8))
  
  
  # find the lambda.1se
  for (i in min:1) {
    if(cve[i]>cve[min]+cvse[min])
      break
  }
  
  if(min==1){
    lambda.1se <- lambdaSeq[1]
    min_1se <- 1
  }else{
    lambda.1se <- lambdaSeq[i+1]
    min_1se <- i+1
  }
  
  
  
  return(list(cve=cve, cvse=cvse, min=min, lambda.min=lambdaSeq[min], min_1se=min_1se, lambda.1se=lambda.1se,
              lambdaSeq=lambdaSeq,  fold=fold))
  
}



#' Eta update using Newton's method with backtracking
#'
#' \code{update_eta_bktk} updates the precision parameter tau = e^eta for L2E regression using Newton's method
#'
#' @param r Vector of residual
#' @param eta Initial estimate of eta
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @return Returns a list object containing the new estimate for eta (scalar),
#' the number of iterations (scalar) the update step utilized,
#' the eta and objective function solution paths (vectors), and
#' the first and second derivatives calculated via Newton's method (vectors)
#'
#'
#' @keywords internal
update_eta_bktk <- function(r, eta, max_iter=1e2, tol=1e-10) {
  
  
  n <- length(r)
  r_sq <- r^2
  stepsize <- 1
  stepsize_shrinkage <- 0.9
  
  
  first_derivative_seq <- double(max_iter)
  second_derivative_seq <- double(max_iter)
  Eta <- double(max_iter)
  Obj <- double(max_iter)
  
  if ( exp(eta)<=0 ) print("non-positive initial tau in update_eta_bktk")
  
  for (i in 1:max_iter) {
    
    eta_last <- eta
    tau_last <- exp(eta_last)  # avoid computing tau_last in the following
    
    # some elements for computing the derivatives
    v1 <- exp(-0.5*tau_last^2*r_sq) # the w_i's
    v2 <- r_sq*v1
    
    if(all(round(v1,2)==1)) break # Check if v1 is close to 1
    
    first_derivative <- tau_last/(2*sqrt(pi)) - tau_last*sqrt(2/pi)*mean(v1)+
      tau_last^3*sqrt(2/pi)*mean(v2)
    first_derivative_seq[i] <- first_derivative
    
    second_derivative <- tau_last/(2*sqrt(pi))+ 4*tau_last^3*sqrt(2/pi)*mean(v2)
    second_derivative_seq[i] <- second_derivative
    
    
    lam <-  first_derivative^2/second_derivative
    if(lam < tol) break
    
    ### backtracking
    dd <- -first_derivative/second_derivative
    f1 <- objective(eta_last + stepsize*dd, r)
    f0 <- objective(eta_last, r)
    
    
    while (f1>f0-0.5*stepsize*lam) {
      stepsize <- stepsize_shrinkage*stepsize
      f1 <- objective(eta_last + stepsize*dd, r)
    }
    
    eta <- eta_last + stepsize*dd
    stepsize <- 1
    
    Eta[i] <- eta
    Obj[i] <- objective(eta, r)
  }
  
  if(i>max_iter) i=i-1
  
  return(list(eta=as.numeric(eta),iter=i, Eta=Eta[1:i], Obj=Obj[1:i],
              first_derivative_seq = first_derivative_seq[1:i],
              second_derivative_seq = second_derivative_seq[1:i]))
  
}



#' Objective function of the L2E regression - tau
#'
#' \code{objective_tau} computes the objective of the L2E regression in terms of tau
#'
#' @param tau The current estimate of tau
#' @param r Vector of residuals
#' @param method Mean or median
#' @return Returns the output of the objective function (scalar)
#' @importFrom stats median
#' @keywords internal
objective_tau <- function(tau, r, method="mean"){
  
  v1 <- exp(-0.5*tau^2*r^2)
  
  s1 <- tau/(2*sqrt(pi))
  
  if(method=="mean"){
    s2 <- tau* sqrt(2/pi)*mean(v1)
  }else{
    s2 <- tau* sqrt(2/pi)*median(v1)
  }
  
  return(s1-s2)
}


#' Objective function of the L2E regression - eta
#'
#' \code{objective} computes the objective of the L2E regression in terms of eta
#'
#' @param eta The current estimate of eta
#' @param r Vector of residuals
#' @param method Mean or median
#' @return Returns the output of the objective function (scalar)
#' @importFrom stats median
#' @keywords internal
objective <- function(eta, r, method="mean"){
  
  v1 <- exp(-0.5*exp(2*eta)*r^2)
  
  s1 <- exp(eta)/(2*sqrt(pi))
  
  if(method=="mean"){
    s2 <- exp(eta)* sqrt(2/pi)*mean(v1)
  }else{
    s2 <- exp(eta)* sqrt(2/pi)*median(v1)
  }
  
  return(s1-s2)
}



