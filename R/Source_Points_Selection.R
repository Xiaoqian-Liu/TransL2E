#' a function to select the source data points that have high similarity level with target set
#'
#' \code{select_source} identifies source data points that are likely to be 
#' informative for the target task by using a density-ratio based rejection 
#' sampling approach.
#'
#' @param target Numeric vector of weights/values from the target domain.
#' @param source Numeric vector of weights/values from the source domain.
#'
#' @return A vector of indices for the source data points that were accepted.
#' 
#' @details The function estimates the density of both domains and calculates 
#' a ratio. Points are accepted with a probability influenced by the 
#' Hellinger distance between distributions and the density ratio 
#' \eqn{f_T(x)/f_S(x)}. A small constant (\eqn{10^{-8}}) is added to the 
#' denominator for numerical stability.
#'
#' @importFrom stats density approx runif
#' @keywords internal


select_source <- function(target, source) {
  
  # size of source and target sample
  n_s = length(source)
  n_t = length(target)
  
  # # 1. Kernel density estimates for source & target
  dens_S <- density(source, from = 0, to = 1)
  dens_T <- density(target, from = 0, to = 1)
  
  # We'll approximate f_S(s_j) and f_T(s_j) at each source weight s_j
  ratio_vec <- numeric(n_s)
  
  for (j in 1:n_s) {
    s_val <- source[j]
    
    # Approximate the densities at s_val
    fS_val <- approx(dens_S$x, dens_S$y, xout = s_val, rule = 2)$y
    fT_val <- approx(dens_T$x, dens_T$y, xout = s_val, rule = 2)$y
    
    ratio_vec[j] <- fT_val * s_val / (fS_val + 10^-8)
  }
  
  # compute hellinger distance
  hd = hellinger_distance(target, source)
  
  # For the candidates, acceptance prob
  accept_prob <- exp(-sqrt(hd)) * pmin(1, ratio_vec)
  
  # Accept based on this probability
  u <- runif(n_s)
  source_index <- which(u <= accept_prob)
  
  return(source_index)
  
}

#' Compute Hellinger Distance between two samples
#'
#' \code{hellinger_distance} calculates the Hellinger distance between the 
#' kernel density estimates of the target and source vectors.
#'
#' @param target Numeric vector from the target distribution.
#' @param source Numeric vector from the source distribution.
#'
#' @return A scalar (numeric) representing the Hellinger distance, 
#' ranging from 0 to 1.
#'
#' 

#' @importFrom stats density approx
#' @keywords internal
hellinger_distance <- function(target, source) {
  
  # Estimate densities
  density_target <- density(target, n = 1000, from = 0, to = 1)
  density_source <- density(source, n = 1000, from = 0, to = 1)
  
  # Align densities on the same grid
  target_y <- density_target$y / sum(density_target$y)  # Normalize target density
  source_y <- approx(density_source$x, density_source$y, xout = density_target$x)$y
  source_y <- source_y / sum(source_y)  # Normalize source density
  
  return(sqrt(sum((sqrt(target_y) - sqrt(source_y))^2)) / sqrt(2))
}



#' Generate an AR(1) Covariance Matrix
#' 
#' Internal helper to simulate a p-by-p covariance matrix with an AR(1) structure.
#' 
#' @param p Integer. The number of features (dimensions).
#' @param rho Numeric. The correlation coefficient, typically between 0 and 1.
#' @return A p-by-p numeric matrix where element (i,j) is rho^|i-j|.
#' @export
source_cov = function(p, rho){

  source_covariance = matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      source_covariance[i,j] = rho^abs(i-j)
    }
  }
  return(source_covariance)
}


#' Simulate Multiple Source Datasets
#' 
#' Internal helper to create a list of source datasets with covariate shift, 
#' model shift, and response outliers.
#' 
#' @param n_s Integer. Sample size for each source dataset.
#' @param p Integer. Total number of features.
#' @param p_s Integer. Number of features affected by the model shift (delta_s).
#' @param snr Numeric. Standard deviation of the Gaussian noise in the response.
#' @param sigma_s Numeric. Standard deviation of the model shift (beta_s - beta_t).
#' @param Sigma_S Matrix. The covariance matrix for the design matrix X.
#' @param K Integer. The number of source datasets to generate.
#' @param out_prop Numeric. Proportion of samples in each source to be corrupted by outliers (0 to 1).
#' @param out_sig Numeric. Magnitude of the outlier shift relative to the maximum Y value.
#' @param beta_t Numeric vector. The ground truth coefficients for the target domain.
#' 
#' @return A list of length K, where each element is a list containing:
#' \item{x}{A design matrix of size n_s by p.}
#' \item{y}{A response vector of length n_s.}
#' 
#' @importFrom stats rnorm
#' @importFrom MASS mvrnorm
#' @export
source_cre <- function(n_s, p, p_s, snr, sigma_s, Sigma_S, K, out_prop, out_sig, beta_t) {

  # Generate a list of K source datasets
  source_list <- lapply(seq_len(K), function(i) {
    # Generate design matrix for the current source
    X <- mvrnorm(n = n_s, mu = rep(0, p), Sigma = Sigma_S)

    # Create a shift vector for beta
    delta_s <- rep(0, p)
    delta_s[1:p_s] <- rnorm(p_s, mean = 0, sd = sigma_s)

    # Adjust beta_t with the current shift to obtain beta_s
    beta_s <- beta_t + delta_s

    # Generate response vector Y with additive noise
    Y <- X %*% beta_s + rnorm(n_s, mean = 0, sd = snr)

    # Introduce outliers to a proportion of the responses
    n_out <- floor(n_s * out_prop)
    Y[1:n_out] <- Y[1:n_out] + out_sig*max(Y)

    # Introduce outliers to a proportion of the covariates
    #X[1:n_out, ] = X[1:n_out, ] + out_sig

    # Return the source dataset in the desired list format
    list(x = X, y = Y)
  })

  return(source_list)
}


