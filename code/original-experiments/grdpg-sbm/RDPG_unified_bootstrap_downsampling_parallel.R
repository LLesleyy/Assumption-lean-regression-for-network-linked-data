
###########################################
#    RDPG: sparse case SOLUTION           #
###########################################

## FIXED Rho Experiment first.

# Required Libraries
library(methods)
library(Rcpp)
library(RcppArmadillo)
library(pracma)
library(igraph)
library(mvtnorm)
library(MASS)
library(abind)
library(doParallel)

# Rcpp Setup: generating RDPG.
setClass("sbm",representation(B_mat="matrix",pi_vec="numeric",nu_seq = "function"))
create_simulation_function =  function(x,return_P_mat =FALSE) UseMethod("create_simulation_function")
create_simulation_function.sbm = function(sbm,return_P_mat =FALSE){
  
  cppFunction( 'arma::mat my_function_1( int n, double nu, NumericMatrix B_mat,
               NumericVector pi_vec ) {
               arma::mat A_mat(n,n,arma::fill::zeros);
               int num_class = pi_vec.size();
               IntegerVector my_classes = seq_len(num_class);
               
               IntegerVector class_assignment = Rcpp::sample(my_classes,n,true,pi_vec);
               
               for (int ii = 0; ii < n; ii++)
               {
               for(int jj = ii+1; jj < n; jj++) 
               {
               int x_temp = class_assignment[ii]-1;
               int y_temp = class_assignment[jj]-1;
               
               double P_temp = nu * B_mat(x_temp, y_temp);
               A_mat(ii,jj) = R::rbinom(1,P_temp);
               }
               }
               A_mat = arma::symmatu(A_mat);
               return(A_mat);
               }
               ',depends = "RcppArmadillo")
  cppFunction( 'List my_function_2( int n, double nu, NumericMatrix B_mat,
               NumericVector pi_vec ) {
               arma::mat A_mat(n,n,arma::fill::zeros);
               arma::mat P_mat(n,n,arma::fill::zeros);
               int num_class = pi_vec.size();
               IntegerVector my_classes = seq_len(num_class);
               
               IntegerVector class_assignment = Rcpp::sample(my_classes,n,true,pi_vec);
               
               for (int ii = 0; ii < n; ii++)
               {
               for(int jj = ii+1; jj < n; jj++) 
               {
               int x_temp = class_assignment[ii]-1;
               int y_temp = class_assignment[jj]-1;
               
               double P_temp = nu * B_mat(x_temp, y_temp);
               P_mat(ii,jj) = P_temp;
               A_mat(ii,jj) = R::rbinom(1,P_temp);
               }
               }
               A_mat = arma::symmatu(A_mat);
               P_mat = arma::symmatu(P_mat);
               List ret;
               ret["A_mat"] = A_mat;
               ret["P_mat"] = P_mat;
               ret["latent_pos"] = class_assignment;
               return ret;
               }
               ',depends = "RcppArmadillo")
  
  if(return_P_mat){  
    output = function(n) {
      return( my_function_2(n=n,nu=sbm@nu_seq(n),B_mat = sbm@B_mat, pi_vec = sbm@pi_vec))
    }
    return(output)
  }
  else{
    output = function(n) {
      return( my_function_1(n=n,nu=sbm@nu_seq(n),B_mat = sbm@B_mat, pi_vec = sbm@pi_vec))
    }
    return(output)
    
  }
}

# Parameters
N_graph = 5000
# only cut down to two sparsity cases.
sparsity = c(-0.65, -0.55)
N_sim = 300 # 0.95 \pm 1.96 \times 0.0117 \approx (0.927, 0.973)
B_resample = 2000  # unified for both bootstrap and down-sampling
K = 3

# True coefficient vector
beta_coeff_case3 <- c(1, 2, 1, 2, 1)

# SBM Model
B_mat = matrix(c(0.8, 0.2, 0.1, 0.2, 0.7, 0.15, 0.1, 0.15, 0.9), nrow = K, byrow = TRUE)
pi_vec = c(0.65, 0.25, 0.1)

# Latent eigen-decomposition.
eigen_B = eigen(B_mat)
ordering_B <- order(abs(eigen_B$values), decreasing = TRUE)
U = eigen_B$vectors[, ordering_B]
S_sqrt_inv = diag(sqrt(abs(eigen_B$values[ordering_B])))

# Parallel
registerDoParallel(cores = 8)

for (k in 1:length(sparsity)) {
  
  cat("Sparsity case:", sparsity[k], "\n")
  nu <- N_graph^sparsity[k]
  m_size <- round(5 * N_graph * nu)
  cat("Down-sample size:", m_size, "\n")
  
  all_downsampling <- list()
  all_bootstrap <- list()
  
  my.sbm_test = new("sbm", B_mat = B_mat, pi_vec = pi_vec, nu_seq = function(x) return(x^sparsity[k]))
  simulate_sbm_test = create_simulation_function(my.sbm_test, return_P_mat = TRUE)
  
  for (sim in 1:N_sim) {
    
    cat("  Simulation:", sim, "/", N_sim, "\n")
    
    # sample observed graph
    sample_graph = simulate_sbm_test(N_graph)
    A_obs = sample_graph$A_mat
    Z_latent = sample_graph$latent_pos
    edge_density <- sum(A_obs) / (N_graph * (N_graph - 1))
    mean_kernel <- drop(crossprod(my.sbm_test@pi_vec, my.sbm_test@B_mat %*% my.sbm_test@pi_vec))
    # Diagnostic estimate of nu using the known simulation kernel.
    nu_obv <- edge_density / mean_kernel
    
    # latent position construction.
    binary_matrix = matrix(0, ncol = K, nrow = length(Z_latent))
    binary_matrix[cbind(seq_along(Z_latent), Z_latent)] = 1
    X_latent = binary_matrix %*% U %*% S_sqrt_inv
    # Clean positions use unscaled B; only raw ASE needs rho normalization.
    X_latent_scaled <- X_latent
    
    # estimated ASE construction.
    eigen_decomp = eigen(A_obs)
    top_indices = order(abs(eigen_decomp$values), decreasing = TRUE)[1:K]
    X_estimated = eigen_decomp$vectors[, top_indices] %*% diag(sqrt(abs(eigen_decomp$values[top_indices])))
    X_estimated_scaled = nu^(-0.5) * X_estimated
    
    # Y generating process.
    cov_matrix = matrix(c(1, 0.3, 0.3, 1), 2)
    bivariate_cov = MASS::mvrnorm(N_graph, mu = c(1, 2), Sigma = cov_matrix)
    X_design_latent = cbind(bivariate_cov, X_latent_scaled)
    X_design_obs = cbind(bivariate_cov, X_estimated_scaled)
    y_obs = X_design_latent %*% beta_coeff_case3 + rnorm(N_graph)
    
    # Compute OLS estimator
    beta_obs <- solve(t(X_design_obs) %*% X_design_obs) %*% t(X_design_obs) %*% matrix(y_obs, nrow = N_graph)
    
    # one-time computation for observed statistics.
    observed_product <- X_design_obs * matrix(y_obs, nrow = N_graph, ncol = ncol(X_design_obs))
    observed_linear <- colMeans(observed_product)
    observed_quad <- crossprod(X_design_obs) / N_graph
    
    ###########################################
    # re-sampling starts: parallel computing. #
    ###########################################
    
    # Attention: if also bootstrap for nu_obv, consider the bootstrap framework.
    
    # Initialize empty matrices for this simulation
    beta_ds_mat <- matrix(NA, nrow = length(beta_coeff_case3), ncol = B_resample)
    beta_boot_mat <- matrix(NA, nrow = length(beta_coeff_case3), ncol = B_resample)
    
    # Parallel foreach
    results = foreach(b = 1:B_resample, .packages = c('base')) %dopar% {
      
      ### Down-sampling ###
      down_idx = sample(N_graph, m_size)
      X_d = X_design_obs[down_idx, , drop = FALSE]
      y_d = y_obs[down_idx]
      beta_ds_b = solve(crossprod(X_d) / m_size, colMeans(X_d * y_d))
      
      #### Bootstrap ###
      epsilon <- rnorm(N_graph, 0, 1)
      linear_res_centered <- sweep(observed_product, 2, observed_linear, `-`)
      perturb_linear_sum <- linear_res_centered * matrix(epsilon, nrow = length(epsilon), ncol = ncol(linear_res_centered))
      
      perturb_squared_sum <- matrix(0, nrow = ncol(X_design_obs), ncol = ncol(X_design_obs))
      for (i in 1:N_graph) {
        xi <- X_design_obs[i, ]
        perturb_squared_sum <- perturb_squared_sum + epsilon[i] * (tcrossprod(xi) - observed_quad)
      }
      
      avg_xy <- observed_linear + colMeans(perturb_linear_sum)
      avg_xx <- observed_quad + perturb_squared_sum / N_graph
      beta_boot_b <- solve(avg_xx, avg_xy)
      
      list(beta_ds = beta_ds_b, beta_boot = beta_boot_b)
    }
    
    # After foreach, unpack results
    for (b in 1:B_resample) {
      beta_ds_mat[, b] <- results[[b]]$beta_ds
      beta_boot_mat[, b] <- results[[b]]$beta_boot
    }
    
    # Now you can save:
    all_downsampling[[sim]] = beta_ds_mat
    all_bootstrap[[sim]] = beta_boot_mat
    
  }
  
  saveRDS(all_downsampling, file = paste0("0425_RDPG_downsampling_sparsity_", k, ".rds"))
  saveRDS(all_bootstrap, file = paste0("0425_RDPG_bootstrap_sparsity_", k, ".rds"))
  cat("Finished sparsity case", sparsity[k], "\n\n")
}



################################
#       coverage results       #
################################

library(ggplot2)

# Setup
sparsity_levels <- sparsity
beta_indices <- 1:2   # β₁ and β₂
true_beta <- c(1,2)
alpha <- 0.05
z_crit <- qnorm(1 - alpha/2)

# Initialize results
coverage_table <- data.frame()

for (k in 1:2) {  
  # sparsity 1 and 2
  
  # Load the down-sampling and bootstrap results
  downsampling_samples <- readRDS(paste0("0425_RDPG_downsampling_sparsity_", k, ".rds"))
  bootstrap_samples <- readRDS(paste0("0425_RDPG_bootstrap_sparsity_", k, ".rds"))
  
  num_sim <- length(downsampling_samples)
  
  # Initialize counters
  coverage_quantile_ds <- rep(0, 2)
  coverage_normal_ds <- rep(0, 2)
  coverage_quantile_boot <- rep(0, 2)
  
  for (sim in 1:num_sim) {
    beta_ds_mat <- downsampling_samples[[sim]]
    beta_boot_mat <- bootstrap_samples[[sim]]
    
    for (j in beta_indices) {
      
      beta_ds_j <- beta_ds_mat[j, ]
      beta_boot_j <- beta_boot_mat[j, ]
      
      true_val <- true_beta[j]
      
      ### Down-sampling CIs
      # Quantile CI
      ci_q_ds <- quantile(beta_ds_j, probs = c(alpha/2, 1 - alpha/2))
      if (true_val >= ci_q_ds[1] && true_val <= ci_q_ds[2]) {
        coverage_quantile_ds[j] <- coverage_quantile_ds[j] + 1
      }
      # Normal CI
      mean_ds <- mean(beta_ds_j)
      sd_ds <- sd(beta_ds_j)
      ci_n_ds <- mean_ds + c(-1, 1) * z_crit * sd_ds
      if (true_val >= ci_n_ds[1] && true_val <= ci_n_ds[2]) {
        coverage_normal_ds[j] <- coverage_normal_ds[j] + 1
      }
      
      ### Bootstrap CIs
      # Centered bootstrap samples
      centered_boot <- beta_boot_j - mean(beta_boot_j)
      # Compute quantiles
      q_centered <- quantile(centered_boot, probs = c(alpha/2, 1 - alpha/2))
      # Shift back (reverse the order)
      ci_bootstrap <- mean(beta_boot_j) - rev(q_centered)
      # Check coverage
      if (true_val >= ci_bootstrap[1] && true_val <= ci_bootstrap[2]) {
        coverage_quantile_boot[j] <- coverage_quantile_boot[j] + 1
      }
      }
    }
  
  # Calculate proportions
  coverage_quantile_ds <- coverage_quantile_ds / num_sim
  coverage_normal_ds <- coverage_normal_ds / num_sim
  coverage_quantile_boot <- coverage_quantile_boot / num_sim
  
  # Build the result table
  temp_table <- data.frame(
    Sparsity = rep(sparsity_levels[k], 2),
    Beta = paste0("beta_", beta_indices),
    True_Value = true_beta,
    Downsample_Quantile_CI = round(coverage_quantile_ds, 3),
    Downsample_Normal_CI = round(coverage_normal_ds, 3),
    Bootstrap_Quantile_CI = round(coverage_quantile_boot, 3)
  )
  
  coverage_table <- rbind(coverage_table, temp_table)
}

print(coverage_table)
