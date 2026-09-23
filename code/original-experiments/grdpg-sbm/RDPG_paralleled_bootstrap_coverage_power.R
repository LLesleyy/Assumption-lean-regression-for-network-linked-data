# Required libraries
library(methods)
library(Rcpp)
library(RcppArmadillo)
library(pracma)
library(igraph)
library(mvtnorm)
library(MASS)
library(doParallel)

set.seed(26)

# Define SBM class
setClass("sbm", representation(B_mat = "matrix", pi_vec = "numeric", nu_seq = "function"))

create_simulation_function <- function(x, return_P_mat = FALSE) UseMethod("create_simulation_function")

create_simulation_function.sbm <- function(sbm, return_P_mat = FALSE) {
  cppFunction('List my_function_2(int n, double nu, NumericMatrix B_mat, NumericVector pi_vec) {
    arma::mat A_mat(n,n,arma::fill::zeros);
    arma::mat P_mat(n,n,arma::fill::zeros);
    int num_class = pi_vec.size();
    IntegerVector my_classes = seq_len(num_class);
    IntegerVector class_assignment = Rcpp::sample(my_classes, n, true, pi_vec);
    for (int ii = 0; ii < n; ii++) {
      for (int jj = ii+1; jj < n; jj++) {
        int x_temp = class_assignment[ii]-1;
        int y_temp = class_assignment[jj]-1;
        double P_temp = nu * B_mat(x_temp, y_temp);
        P_mat(ii,jj) = P_temp;
        A_mat(ii,jj) = R::rbinom(1, P_temp);
      }
    }
    A_mat = arma::symmatu(A_mat);
    P_mat = arma::symmatu(P_mat);
    List ret;
    ret["A_mat"] = A_mat;
    ret["P_mat"] = P_mat;
    ret["latent_pos"] = class_assignment;
    return ret;
  }', depends = "RcppArmadillo")
  
  function(n) {
    my_function_2(n = n, nu = sbm@nu_seq(n), B_mat = sbm@B_mat, pi_vec = sbm@pi_vec)
  }
}

# Simulation settings
N_graph <- 2000
sparsity <- c(-0.75, -0.25)
k <- 2
nu <- N_graph^(sparsity[k])
N_sim <- 200
B_bootstrap <- 1000
K <- 3

my.sbm_test <- new("sbm", B_mat = matrix(c(0.8, 0.2, 0.1,
                                           0.2, 0.7, 0.15,
                                           0.1, 0.15, 0.9), nrow = K, byrow = TRUE),
                   pi_vec = c(0.65, 0.25, 0.1),
                   nu_seq = function(x) return(x^(sparsity[k])))

simulate_sbm_test <- create_simulation_function(my.sbm_test, return_P_mat = TRUE)
eigen_B <- eigen(my.sbm_test@B_mat)
U <- eigen_B$vectors[, order(abs(eigen_B$values), decreasing = TRUE)]
S <- sort(abs(eigen_B$values), decreasing = TRUE)
S_sqrt_inv <- diag(sqrt(S))
beta_coeff_case1 <- c(1, 2, 1, 2, 1)

all_bootstrap_samples_RDPG <- list()
rejections <- numeric(N_sim)

registerDoParallel(cores = 8)

for (sim in 1:N_sim) {
  cat("Starting repetition", sim, "of", N_sim, "under sparsity", sparsity[k], "\n")
  sample_graph <- simulate_sbm_test(N_graph)
  A_obs <- sample_graph[['A_mat']]
  Z_latent <- sample_graph[['latent_pos']]
  edge_density <- sum(A_obs) / (N_graph * (N_graph - 1))
  mean_kernel <- drop(crossprod(my.sbm_test@pi_vec, my.sbm_test@B_mat %*% my.sbm_test@pi_vec))
  # Diagnostic estimate of nu using the known simulation kernel.
  nu_obv <- edge_density / mean_kernel
  
  binary_matrix <- matrix(0, ncol = K, nrow = N_graph)
  binary_matrix[cbind(seq_along(Z_latent), Z_latent)] <- 1
  X_latent <- binary_matrix %*% U %*% S_sqrt_inv
  # Clean positions use unscaled B; only raw ASE needs rho normalization.
  X_latent_scaled <- X_latent
  
  eig_A <- eigen(A_obs)
  top_indices <- order(abs(eig_A$values), decreasing = TRUE)[1:K]
  X_estimated <- eig_A$vectors[, top_indices] %*% diag(sqrt(abs(eig_A$values[top_indices])))
  X_estimated_scaled <- nu^(-0.5) * X_estimated
  
  cov_matrix <- matrix(c(1, 0.3, 0.3, 1), 2)
  bivar_cov <- MASS::mvrnorm(N_graph, mu = c(1, 2), Sigma = cov_matrix)
  
  X_design_latent <- cbind(bivar_cov, X_latent_scaled)
  X_design_obs <- cbind(bivar_cov, X_estimated_scaled)
  y_obs <- X_design_latent %*% matrix(beta_coeff_case1, ncol = 1) + rnorm(N_graph)
  
  beta_obs <- solve(t(X_design_obs) %*% X_design_obs) %*% t(X_design_obs) %*% matrix(y_obs)
  
  # observed_product <- sweep(X_design_obs, 1, y_obs, `*`)  # n x p
  observed_product <- X_design_obs * matrix(y_obs, nrow = N_graph, ncol = ncol(X_design_obs))
  observed_linear <- colMeans(observed_product)
  observed_quad <- crossprod(X_design_obs) / N_graph
  
  beta_boot <- foreach(B = 1:B_bootstrap, .combine = 'cbind') %dopar% {
    epsilon <- rnorm(N_graph)
    
    # Step 2: Subtract observed_linear (mean of x_i * y_i across samples)
    xi1_centered <- sweep(observed_product, 2, observed_linear, `-`)  # n x p
    
    # Step 3: Multiply each row by epsilon
    xi1 <- sweep(xi1_centered, 1, epsilon, `*`)  # n x p
    
    xi2_sum <- matrix(0, nrow = 5, ncol = 5)
    for (i in 1:N_graph) {
      xi <- X_design_obs[i, ]
      xi2_sum <- xi2_sum + epsilon[i] * (tcrossprod(xi) - observed_quad)
    }
    avg_x <- observed_linear + colMeans(xi1)
    avg_xx <- observed_quad + xi2_sum / N_graph
    solve(avg_xx, avg_x)
  }
  
  all_bootstrap_samples_RDPG[[sim]] <- beta_boot
  Tn_boot <- colSums((beta_boot[3:5, ] - beta_obs[3:5])^2) * N_graph
  Tn_obs <- sum(beta_obs[3:5]^2) * N_graph
  # Upper-tail test; use the manuscript's 95th-percentile rejection rule.
  p_val <- (1 + sum(Tn_boot >= Tn_obs)) / (1 + length(Tn_boot))
  critical_value <- unname(quantile(Tn_boot, 0.95))
  rejections[sim] <- (Tn_obs > critical_value)
  print(paste("Finished Running iteration:", sim))
}

power_estimate <- mean(rejections)
print(paste("Estimated power:", power_estimate))
saveRDS(all_bootstrap_samples_RDPG, file = paste0("0415_parallel_RDPG_bootstrap_samples_sparsity_", sparsity[k], ".rds"))
