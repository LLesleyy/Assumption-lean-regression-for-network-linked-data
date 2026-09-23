# Load necessary libraries
library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)

set.seed(26)

# Load the C++ functions
sourceCpp("optimized_functions.cpp", rebuild = TRUE, verbose = TRUE)

# Simulation parameters
num_approx<-10^4; # as large as possible 
graphon_size <- 1000;
sparsity <- c(-0.75, -0.6, -0.55, -0.45, -0.4, -0.25)
k <- 6;nu <- graphon_size^(sparsity[k])
beta0<-1;beta1<-20; 
# beta2<-3;beta3<-2;beta_vec = matrix(c(beta2, beta3), ncol = 1);


# Initialize storage for beta_hat, beta_hat_par_cor, and beta_true
target_approx_results <- data.frame(
  beta_hat = numeric(num_approx),
  beta_hat_valid = numeric(num_approx)
)

for (iteration in 1:num_approx) {
  
  print(paste("Iteration:", iteration))
  
  # Generate X_vec
  X_vec <- rmvnorm(graphon_size, mean = c(0, 20), sigma = matrix(c(1, 0.4, 0.4, 1), nrow = 2))
  latent_positions <- X_vec[, 1]
  X1 <- X_vec[, 2]
  
  # Generate P matrix and adjacency matrix
  P_mat_temp <- generate_P_mat(latent_positions)
  A_mat <- simulate_from_P_mat(nu * P_mat_temp)
  P_mat <- nu * P_mat_temp
  di <- colSums(A_mat) # node-degree
  
  # Compute f_xi_Numerator and f_xi_Denominator in parallel
  f_xi_Numerator <- sapply(latent_positions, compute_neighborhood_avg_v2)
  f_xi_Denominator <- sapply(latent_positions, compute_kernel_avg_v2)
  # Compute f_xi
  f_xi <- f_xi_Numerator / f_xi_Denominator
  
  # Calculate neighborhood averages
  neighborhood_avg <- (A_mat %*% X1) / rowSums(A_mat)
  neighborhood_avg[di == 0] <- 0
  
  # Generate independent y
  error <- rnorm(graphon_size, 0, 0.5)
  y_gen <- beta0 + beta1 * f_xi + error # still under indepdent y generating process
  
  # Pair network data
  gen_data <- data.frame(y_gen = y_gen, x1_neighbourhood_ave = neighborhood_avg)
  colnames(gen_data) = c('y_gen', 'x1_neighbourhood_ave')
  
  # Compute beta_hat
  yzhat_mean <- mean(gen_data$y_gen * gen_data$x1_neighbourhood_ave)
  y_mean <- mean(gen_data$y_gen)
  zhat_mean <- mean(gen_data$x1_neighbourhood_ave)
  zhat2_mean <- mean(gen_data$x1_neighbourhood_ave^2)
  beta_hat <- (yzhat_mean - y_mean * zhat_mean) / (zhat2_mean - zhat_mean^2)
  
  # beta_hat for calling functions
  beta_hat_valid<-lm(y_gen~x1_neighbourhood_ave,data=gen_data)$coefficients[2]
  
  # Store results
  target_approx_results[iteration, ] <- c(beta_hat, beta_hat_valid)
}


