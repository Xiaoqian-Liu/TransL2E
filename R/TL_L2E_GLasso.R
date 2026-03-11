#' Transfer Learning for Robust L2E Group Sparse Regression
#'
#' \code{TL_L2E_glasso} performs robust transfer learning by detecting 
#' informative source datasets, co-learning on fused data, and optionally 
#' applying a debiasing step.
#'
#' @param y_t Response vector of the target dataset.
#' @param X_t Design matrix of the target dataset.
#' @param source_list A list of lists, where each element contains \code{y} (response) 
#' and \code{x} (design matrix) for a source dataset.
#' @param group Group indicator vector for the predictors.
#' @param nfolds Number of cross-validation folds. Default is 5.
#' @param penalty Penalty for sparsity (e.g., "grLasso", "grMCP", "grSCAD").
#' @param max_iter Maximum iterations for the optimization. Default is 100.
#' @param tol Relative tolerance for convergence. Default is 1e-04.
#' @param core Number of CPU cores for parallel computing during source detection.
#' @param seed Seed for reproducibility. Default is 123.
#' @param debias Logical; if TRUE, performs a debiasing step on the target residuals.
#' @param lambda_detect Character; which lambda to use for source detection ("lambda.1se" or "lambda.min").
#' @param lambda_colearn Character; which lambda to use for co-learning ("lambda.1se" or "lambda.min").
#' @param lambda_debias Character; which lambda to use for debiasing ("lambda.1se" or "lambda.min").
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Coefficients}: The final estimated beta coefficients.
#'   \item \code{Selected_Source}: A list of indices and data selected from each source.
#' }
#' 
#' @details This function implements a three-step transfer learning workflow:
#' \enumerate{
#'   \item \strong{Source Detection}: Uses parallel processing to evaluate each source dataset's 
#'   similarity to the target using Hellinger distance and density-ratio weights.
#'   \item \strong{Co-learning}: Aggregates the target data with selected source data to 
#'   refine the group lasso estimates.
#'   \item \strong{Debiasing}: If enabled, fits a sparse L2E model on the residuals to 
#'   correct for transfer bias.
#' }
#'
#' 
#'
#' @importFrom parallel makeCluster stopCluster 
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom doRNG %dorng% registerDoRNG
#' @export
#' @examples 
#' 
#' set.seed(123)
#' sigma_snr = 1
#' 
#' n_t = 200 # target sample size
#' n_s = 400 # source sample size
#' 
#' groups <- 10      # Number of groups
#' p_per_group <- 10  # Number of predictors per group
#' p <- groups * p_per_group  # Total number of predictors
#' 
#' group <- rep(1:groups, each = p_per_group)
#' 
#' 
#' 
#' K = 2 # number of source
#' rho = 0.5 # level of covariate shift
#' Sigma = source_cov(p, rho)
#' 
#' beta_t <- rep(0, p)
#' beta_t[group == 1] <- sample(c(-1,1), p_per_group, replace = TRUE, prob = c(0.5, 0.5))
#' beta_t[group == 3] <- sample(c(-1,1), p_per_group, replace = TRUE, prob = c(0.5, 0.5))
#' 
#' sigma_ms = 0.2 # model shift level
#' p_s = 0.3*p
#' 
#' out_prop = 0.2
#' out_sig = 2
#' 
#' 
#' # simulate the target dataset
#' X_t = matrix(rnorm(n_t * p), nrow = n_t, ncol = p) # target covariates
#' y_t = X_t %*% beta_t + rnorm(n_t, mean = 0, sd = sigma_snr) # target response
#'     
#'     
#' # 10% outler to the samples
#' y_t[1:(n_t * 0.1)] = y_t[1:(n_t * 0.1)] + out_sig*max(y_t)
#' 
#' source_data = source_cre(n_s = n_s, p = p, p_s = p_s, 
#'                       snr = sigma_snr, sigma_s = sigma_ms, Sigma_S = Sigma,
#'                       K = K, out_prop = out_prop, out_sig = out_sig, beta_t = beta_t)
#'                       
#'                       
#' TL_with_debias = TL_L2E_glasso(y_t, X_t, source_list = source_data, group = group, 
#'                               debias = TRUE, lambda_detect = "lambda.1se",
#'                               lambda_colearn="lambda.1se", lambda_debias="lambda.1se")
#' 
#' 
#'
#' # the fitted error
#' sum((beta_t - TL_with_debias$Coefficients)^2)
#' 
TL_L2E_glasso <- function(y_t, X_t, source_list, group,
                          nfolds = 5, penalty = "grLasso", 
                          max_iter = 100, tol = 1e-04, core = 1, seed = 123,
                          debias = TRUE, lambda_detect = "lambda.1se",
                          lambda_colearn="lambda.1se", lambda_debias="lambda.min"){
  
  
  # Basic input checks
  if(nrow(X_t) != nrow(y_t)){
    stop("The number of rows in X_t must match the number of rows in y_t.")
  }
  
  n_t <- nrow(y_t)  # Target sample size
  
  message("Step 1: Source Detection")
  
  # set up the number of clusters
  num_cores <- core
  
  cl <- makeCluster(num_cores)
  
  registerDoParallel(cl)
  
  # Stop the cluster on exit
  on.exit({
    try(stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  registerDoRNG(seed)
  
  loop_results <- foreach(i = seq_along(source_list),
                          .packages = c("TransL2E", "grpreg", "glmnet", "stats", "MASS", "Matrix"),
                          .export = c(
                            "y_t","X_t","source_list","n_t",
                            "penalty","group","max_iter","tol","nfolds","lambda_detect"
                          )) %dorng% {
                            
    select_source <- utils::getFromNamespace("select_source", "TransL2E")
    
    # start parallel computing
    current_source <- source_list[[i]]
    
    # Define number of samples in the current source
    ns_cur <- nrow(current_source$y)
    
    # Combine target and current source data
    y_cur <- rbind(y_t, current_source$y)
    X_cur <- rbind(X_t, current_source$x)
    
    # Fit the L2E model using cross-validation for tuning
    l2e_cur <- CV_L2E_group(y_cur, X_cur, penalty = penalty, group = group,
                                 max_iter = max_iter, tol = tol, 
                                 nfolds = nfolds, trace = FALSE)
    
    if(lambda_detect == "lambda.1se"){
      lambda_val <- l2e_cur$lambda.1se
    } else {
      lambda_val <- l2e_cur$lambda.min
    }
    
    l2e_cur_opt <- L2E_group(y_cur, X_cur, group = group, lambdaSeq = lambda_val, 
                                  penalty = penalty, Show.Time = FALSE)
    
    l2e_fit <- X_cur %*% l2e_cur_opt$Beta
    
    # Compute residuals and weights
    r_cur <- y_cur - l2e_fit
    w_cur <- as.vector(exp(-0.5 * (l2e_cur_opt$Tau * r_cur)^2))
    
    # Select source points
    index_cur <- select_source(w_cur[1:n_t], w_cur[(n_t + 1):(n_t + ns_cur)])
    
    # Prepare the data to be returned for this iteration
    y_selected <- current_source$y[index_cur, , drop = FALSE]
    x_selected <- current_source$x[index_cur,]
    
    # Return a list containing all the data we need.
    # This list will become one element in the final 'loop_results'.
    return(list(
      y = y_selected,
      x = x_selected,
      index = index_cur
    ))
    
  }
  
  
  selected_source <- loop_results
  
  # obtain the final fused x and y
  y_fuse <- y_t
  X_fuse <- X_t
  
  # Use lapply to extract the 'y' and 'x' components from each result,
  # then use do.call(rbind, ...) to stack them all into a big matrix.
  
  all_selected_y <- do.call(rbind, lapply(selected_source, `[[`, "y"))
  all_selected_x <- do.call(rbind, lapply(selected_source, `[[`, "x"))
  
  # Add the selected data to the target data
  y_fuse <- rbind(y_fuse, all_selected_y)
  X_fuse <- rbind(X_fuse, all_selected_x)
  
  
  message("Step 2: Co-learning")
  
  # Fit the L2E model on the fused dataset
  l2e_cl <- CV_L2E_group(y_fuse, X_fuse, penalty = penalty,  group = group,
                              max_iter = max_iter, tol = tol, nfolds = nfolds, trace = FALSE)
  
  if(lambda_colearn=="lambda.1se"){
    lambda_colearn = l2e_cl$lambda.1se
  }else{
    lambda_colearn = l2e_cl$lambda.min
  }
  
  l2e_cl_opt <- L2E_group(y_fuse, X_fuse, group = group,
                          lambdaSeq = lambda_colearn, penalty = penalty,
                               Show.Time = FALSE)
  
  # Estimated beta coefficients
  final_beta <- l2e_cl_opt$Beta
  
 # Step 3: Debias Step
  if(debias){
    message("Step 3: Debiasing")
    
    # Compute residuals on target data
    res_y <- y_t - X_t %*% final_beta
    
    # # Fit the L2E model on the residuals
    l2e_db <- CV_L2E_sparse_ncv(res_y, X_t, penalty = 'lasso',
                                # beta0 = as.vector(coef(cv.glmnet(x = X_t, y = res_y, alpha = 1,
                                #                               intercept = FALSE), s = "lambda.min"))[-1],
                                max_iter = max_iter, tol = tol, nfolds = nfolds, trace = FALSE)
    
    
    if(lambda_debias=="lambda.1se"){
      lambda_debias = l2e_db$lambda.1se
    }else{
      lambda_debias = l2e_db$lambda.min
    }
    
    l2e_db_opt <- L2E_sparse_ncv(res_y, X_t, lambdaSeq = lambda_debias, penalty = 'lasso',
                                 Show.Time = FALSE)
    
    
    # Update the final beta coefficients with the debias step
    final_beta <- final_beta + l2e_db_opt$Beta
  }
  
  return(list(Coefficients = final_beta, Selected_Source = selected_source))
}