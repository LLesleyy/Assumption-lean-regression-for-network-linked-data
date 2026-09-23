
# libraries.
library(methods)
library(Rcpp)
library(RcppArmadillo)
library(pracma)
library(igraph)
library(mvtnorm)
library(MASS)
library(abind)

set.seed(26)

# generating RDPG.
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


########################
# simulation settings. #
########################

# smaller scale due to network noise.
N_graph=2000;
sparsity<-c(-0.75,-0.25)
k=3;
nu<-N_graph^(sparsity[k])

# Number of repetitions
N_sim <- 200
# Number of bootstrap samples
B_bootstrap <- 1000  


# RDPG Model.
K <- 3
# different SBM test cases
my.sbm_test = new("sbm",B_mat= matrix(c(0.8, 0.2, 0.1,
                                        0.2, 0.7, 0.15,
                                        0.1, 0.15, 0.9), nrow = K, byrow = TRUE), 
                  pi_vec=c(0.65,0.25,0.1),
                  nu_seq = function(x) return(x^(sparsity[k])))
# my.sbm_test1 = new("sbm",B_mat= matrix(c(0.02,0.03,0.03,0.01),nrow=2), 
#                  pi_vec=c(0.2,0.8),
#                  nu_seq=function(x)return(1))




#############################################################
# low rank structure and spectral decomposition on global B #
#############################################################

simulate_sbm_test=create_simulation_function(my.sbm_test,return_P_mat = TRUE)
B_mat= matrix(c(0.8, 0.2, 0.1,
                0.2, 0.7, 0.15,
                0.1, 0.15, 0.9), nrow = K, byrow = TRUE)

# Perform eigenvalue decomposition on B_mat
eigen_B <- eigen(B_mat)
# Extract the eigen-vectors and eigenvalues
U <- eigen_B$vectors
S <- eigen_B$values
# Order the eigenvalues in descending order of magnitude
order_indices <- order(abs(S), decreasing = TRUE)
# Reorder the eigenvalues and corresponding eigen-vectors
S <- S[order_indices]
U <- U[, order_indices]
# Construct the eigenvalue embedding by taking the top components
# Use the absolute values and the square root of the eigenvalues
S_sqrt_inv <- diag(sqrt(abs(S)))


######################################
#    prepare regression effects      #
######################################

#beta_coeff_case1 <- c(1, 2, 1, 0.5, 0.1) # strong network effect
beta_coeff_case1 <- c(1, 2, 1, 2, 1)     # strong network effect


#beta_coeff_case4 <- c(1, 2, 0, 0.5, 0.1) # medium network effect
beta_coeff_case2 <- c(1, 2, 0, 0.5, 0.01)   # medium network effect


beta_coeff_case3 <- c(1, 2, 0, 0.01, 0)    # weak network effect
#beta_coeff_case8 <- c(1, 2, 0, 0.01, 0.01) # weak network effect


############################################
########### Repetition Starts ##############
############################################

# 1. bootstrap coverage of conventional covariate will be discovered.
# 2. testing power of existence of network effect will be discovered.

# Store rejection results
rejections <- numeric(N_sim)
# Store all bootstrap samples (dimensions: beta_coeff x B_bootstrap x N_sim)
all_bootstrap_samples_RDPG <- list()

# Estimated rho_n will be used along with this experiment.
for (sim in 1:N_sim) {
  
  cat("Starting repetition", sim, "of", N_sim, "under sparsity",sparsity[k], "\n")
  
  # Generate the observed network data.
  sample_graph <- simulate_sbm_test(N_graph)
  A_obs <- sample_graph[['A_mat']]
  Z_latent_test <- sample_graph[['latent_pos']]
  P_generate <- sample_graph[['P_mat']]
  
  # estimated sparsity level.
  edge_density <- sum(A_obs) / (N_graph * (N_graph - 1))
  mean_kernel <- drop(crossprod(my.sbm_test@pi_vec, my.sbm_test@B_mat %*% my.sbm_test@pi_vec))
  # Diagnostic estimate of nu using the known simulation kernel.
  nu_obv <- edge_density / mean_kernel
  
  # Construct latent positions for current reliazation.
  binary_matrix <- matrix(0, ncol = K, nrow = length(Z_latent_test))
  binary_matrix[cbind(seq_along(Z_latent_test), Z_latent_test)] <- 1
  X_latent <- binary_matrix %*% U %*% S_sqrt_inv
  # Clean positions use unscaled B; only raw ASE needs rho normalization.
  X_latent_scaled <- X_latent
  
  # Compute ASE normalized by the known graph-generation sparsity factor.
  eigen_decomposition <- eigen(A_obs)
  sorted_indices <- order(abs(eigen_decomposition$values), decreasing = TRUE)
  sorted_eigenvalues <- eigen_decomposition$values[sorted_indices]
  sorted_eigenvectors <- eigen_decomposition$vectors[, sorted_indices]
  top_eigenvalues <- sorted_eigenvalues[1:K]
  top_eigenvectors <- sorted_eigenvectors[, 1:K]
  sqrt_lambda <- diag(sqrt(abs(top_eigenvalues)))
  X_estimated <- top_eigenvectors %*% sqrt_lambda
  X_estimated_scaled <- nu^(-0.5) * X_estimated
  #X_estimated_scaled <- nu_obv^(-0.5) * X_estimated
  
  # in case we want to use the estiamted rho_n and also bootstrap for it later.
  # If using nu_obv later, jointly bootstrap its degree projection:
  # g1_hat_nu <- rowSums(A_obs) / ((N_graph - 1) * mean_kernel) - nu_obv
  
  
  #######################################################
  # Data generating process and solve for observed beta #
  #######################################################
  
  # Generate two weakly correlated bivariate conventional covariates
  cov_matrix <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)  
  bivariate_covariates <- MASS::mvrnorm(N_graph, mu = c(1, 2), Sigma = cov_matrix)  
  
  # Generate observed response variable
  X_design_latent <- cbind(bivariate_covariates, X_latent_scaled)
  X_design_observed <- cbind(bivariate_covariates, X_estimated_scaled)
  y_obs <- X_design_latent %*% matrix(beta_coeff_case1, ncol = 1) + rnorm(N_graph, 0, 1)
  
  # Compute OLS estimator
  beta_obs <- solve(t(X_design_observed) %*% X_design_observed) %*% t(X_design_observed) %*% matrix(y_obs, nrow = N_graph)
  
  ##############################
  # Bootstrap procedure starts #
  ##############################
  beta_boot <- matrix(0, nrow = length(beta_coeff_case1), ncol = B_bootstrap)  # Store beta for each bootstrap
  
  for (B in 1:B_bootstrap) {
    
    cat("Starting bootstrap", B, "of", B_bootstrap, "\n")
    
    # Attention: if also bootstrap for nu_obv, consider the bootstrap framework
    
    epsilon <- rnorm(N_graph, 0, 1)  # Gaussian multiplier
    
    # Generate two observed statistics
    observed_product <- matrix(0, nrow = N_graph, ncol = length(beta_coeff_case1))
    for (ii in 1:N_graph) {
      observed_product[ii, ] <- X_design_observed[ii, ] * y_obs[ii]  
    }
    observed_linear <- colMeans(observed_product)
    observed_quad <- (1 / N_graph) * Reduce(`+`, lapply(1:N_graph, function(i) outer(X_design_observed[i, ], X_design_observed[i, ])))
    
    # Perturbation for bootstrap
    xi_boot_perturb_linear <- matrix(0, nrow = N_graph, ncol = length(beta_coeff_case1))
    xi_boot_pertub_quad <- vector("list", N_graph)
    
    for (kk in 1:N_graph) {
      xi_obs <- X_design_observed[kk, ]
      xi2_obs_mat <- outer(xi_obs, xi_obs)
      xi_boot_perturb_linear[kk, ] <- epsilon[kk] * (xi_obs * y_obs[kk] - observed_linear)
      xi_boot_pertub_quad[[kk]] <- epsilon[kk] * (xi2_obs_mat - observed_quad)
    }
    
    # Generate bootstrap sample
    average_x_boot <- observed_linear + colMeans(xi_boot_perturb_linear)
    average_x2_boot <- observed_quad + apply(abind(xi_boot_pertub_quad, along = 3), c(1, 2), mean)
    
    # Compute bootstrapped beta coefficients
    beta_boot[, B] <- solve(average_x2_boot) %*% matrix(average_x_boot, ncol = 1)
  }
  
  # Store bootstrap samples for this simulation
  all_bootstrap_samples_RDPG[[sim]] <- beta_boot  # Store all bootstrap samples for this iteration
  
  
  ############################################################
  # coverage and power computation along with the experiment #
  ############################################################
  
  # Compute test statistics for bootstrap
  Tn_boot_H1 <- numeric(B_bootstrap)
  for (B in 1:B_bootstrap) {
    Tn_boot_H1[B] = N_graph * sum((beta_boot[3:5, B] - beta_obs[3:5])^2)
  }
  Tn_obs_H1 = N_graph * sum(beta_obs[3:5]^2)
  
  # Upper-tail bootstrap test, as in the manuscript's network-effect corollary.
  # The +1 adjustment avoids reporting zero from finitely many bootstrap draws.
  p_value <- (1 + sum(Tn_boot_H1 >= Tn_obs_H1)) / (1 + length(Tn_boot_H1))
  critical_value <- unname(quantile(Tn_boot_H1, 0.95))
  
  # Store rejection decision
  rejections[sim] <- (Tn_obs_H1 > critical_value)
  
  print(paste("Finished Running iteration:", sim))
}

# Compute empirical power
power_estimate <- mean(rejections)
print(paste("Estimated power:", power_estimate))

# save all simulation results.
saveRDS(all_bootstrap_samples_RDPG, file = paste0("0415_RDPG_bootstrap_samples_sparsity_", sparsity[k], ".rds"))



