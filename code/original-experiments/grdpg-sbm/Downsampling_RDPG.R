
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

N_graph=3000;
sparsity<-c(-0.75,-0.7,-0.6)
k=1;
nu<-N_graph^(sparsity[k])


# RDPG Model.
K <- 3
# different SBM test cases
my.sbm_test = new("sbm",B_mat= matrix(c(0.8, 0.2, 0.1,
                                        0.2, 0.7, 0.15,
                                        0.1, 0.15, 0.9), nrow = K, byrow = TRUE), 
                  pi_vec=c(0.65,0.25,0.1),
                  nu_seq=function(x)return(x^(sparsity[k])))
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

beta_coeff_case3 <- c(1, 2, 0, 0.01, 0)    # weak network effect
#beta_coeff_case8 <- c(1, 2, 0, 0.01, 0.01) # weak network effect

############################################
########### Repetition Starts ##############
############################################

# 1. bootstrap coverage of conventional covariate will be discovered.
# 2. testing power of existence of network effect will be discovered.

# Number of repetitions
N_sim <- 300
# Number down-sampling replicates.
D_sampling <- 2000  
# Down-sampling_proportion.
m=0.2
# Down-sample_size.
m_size=N_graph*m
# Store all bootstrap samples (dimensions: beta_coeff x B_bootstrap x N_sim)
all_downsampling_samples_RDPG <- list()

# Estimated rho_n will be used along with this experiment.
for (sim in 1:N_sim) {
  
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
  # down-sampling procedure starts #
  ##############################

  beta_ds_RDPG <- matrix(0, nrow = length(beta_coeff_case1), ncol = D_sampling)
  
  for (D in 1:D_sampling) {
    down_sample_index <- sample(N_graph, m_size)
    
    # Subset design matrix and response
    X_down <- X_design_observed[down_sample_index, , drop = FALSE]
    y_down <- y_obs[down_sample_index]
    
    # Compute down-sampled sufficient statistics
    average_x2_ds <- crossprod(X_down) / m_size
    average_xy_ds <- colMeans(X_down * y_down)
    
    # Solve for down-sampled beta
    beta_ds_RDPG[, D] <- solve(average_x2_ds, average_xy_ds)
  }
  
  # Store ds samples for this simulation
  all_downsampling_samples_RDPG[[sim]] <- beta_ds_RDPG  # Store all bootstrap samples for this iteration
  
}
  
# save all simulation results.
saveRDS(all_downsampling_samples_RDPG, file = paste0("0420_RDPG_Downsampling_samples_sparsity_", k, ".rds"))




# the coverage in the sparse regime of RDPG
# Load the saved downsampling results
all_downsampling_samples_RDPG <- readRDS("0420_RDPG_Downsampling_samples_sparsity_1.rds")

# Set parameters
beta_true <- c(1, 2, 1, 2, 1)
alpha <- 0.05
z_crit <- qnorm(1 - alpha / 2)
num_sim <- length(all_downsampling_samples_RDPG)
num_beta <- length(beta_true)

# Initialize vectors to store coverage counts
coverage_quantile <- rep(0, num_beta)
coverage_normal <- rep(0, num_beta)

# Loop through simulations
for (sim in 1:num_sim) {
  beta_mat <- all_downsampling_samples_RDPG[[sim]]
  
  for (j in 1:num_beta) {
    beta_j <- beta_mat[j, ]
    true_val <- beta_true[j]
    
    # Quantile CI
    ci_q <- quantile(beta_j, probs = c(alpha / 2, 1 - alpha / 2))
    if (true_val >= ci_q[1] && true_val <= ci_q[2]) {
      coverage_quantile[j] <- coverage_quantile[j] + 1
    }
    
    # Normal CI
    mean_j <- mean(beta_j)
    sd_j <- sd(beta_j)
    ci_n <- mean_j + c(-1, 1) * z_crit * sd_j
    if (true_val >= ci_n[1] && true_val <= ci_n[2]) {
      coverage_normal[j] <- coverage_normal[j] + 1
    }
  }
}

# Convert to proportions
coverage_quantile <- coverage_quantile / num_sim
coverage_normal <- coverage_normal / num_sim

# Combine results
coverage_df <- data.frame(
  Beta_Index = 1:num_beta,
  True_Beta = beta_true,
  Quantile_CI_Coverage = round(coverage_quantile, 3),
  Normal_CI_Coverage = round(coverage_normal, 3)
)

print(coverage_df)

#> print(coverage_result)
#Beta True_Value Quantile_CI_Coverage Normal_CI_Coverage
#1 beta_1          1                0.850              0.843
#2 beta_2          2                0.617              0.617
