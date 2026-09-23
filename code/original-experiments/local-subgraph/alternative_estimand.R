

## pure approximation scheme.
## a separate experiment only targeting the population quantity.

library(methods)
library(Rcpp)
library(RcppArmadillo)
library(pracma)
library(igraph)
library(mvtnorm)
library(Matrix)
library(combinat)

## algorithm for searching the Rooted two-star/size-4-wheel and Peripheral two-star/size-4-wheel
# Set object classes for SBM, Graphon, LSM
setClass("graphon", representation(w_function = "character", nu_seq = "function"))
create_simulation_function = function(x, return_P_mat = FALSE) UseMethod("create_simulation_function")

create_simulation_function.graphon = function(graphon, return_P_mat = FALSE) {
  
  cppFunction('arma::mat my_function_1(int n, double nu) {
               arma::mat A_mat(n, n, arma::fill::zeros);
               NumericVector latent_positions = runif(n, 0, 1);
               for (int ii = 0; ii < n; ii++) {
                 for (int jj = ii + 1; jj < n; jj++) {
                   double P_temp = nu * w_function(latent_positions[ii], latent_positions[jj]);
                   A_mat(ii, jj) = R::rbinom(1, P_temp);
                 }
               }
               A_mat = arma::symmatu(A_mat);
               return(A_mat);
}', depends = "RcppArmadillo", includes = graphon@w_function)
  
  cppFunction('List my_function_2(int n, double nu) {
               arma::mat A_mat(n, n, arma::fill::zeros);
               arma::mat P_mat(n, n, arma::fill::zeros);
               NumericVector latent_positions = runif(n, 0, 1);
               for (int ii = 0; ii < n; ii++) {
                 for (int jj = ii + 1; jj < n; jj++) {
                   double P_temp = nu * w_function(latent_positions[ii], latent_positions[jj]);
                   A_mat(ii, jj) = R::rbinom(1, P_temp);
                   P_mat(ii, jj) = P_temp;
                 }
               }
               A_mat = arma::symmatu(A_mat);
               P_mat = arma::symmatu(P_mat);
               List ret;
               ret["A_mat"] = A_mat;
               ret["P_mat"] = P_mat;
               ret["latent_positions"] = latent_positions;
               return ret;
}', depends = "RcppArmadillo", includes = graphon@w_function)
  
  if (return_P_mat) {
    output = function(n) {
      return(my_function_2(n = n, nu = graphon@nu_seq(n)))
    }
    return(output)
  } else {
    output = function(n) {
      return(my_function_1(n = n, nu = graphon@nu_seq(n)))
    }
    return(output)
  }
}


#Comments: rank 2.
graphon.22 = 'double w_function(double u, double v) {
return(0.5/0.5*(u+v));
}'

# specify kernel function.
w_function.22 <- function(u, v) {
  return(0.5/0.5*(u+v))
}

# normalization constant.
integral2(w_function.22, 0, 1, 0, 1)$Q

# define f(xi) function.
f_xi_func.22_tri <- function(xi) {
  integrate_function <- function(u, v) {
    return(w_function.22(xi, u) * w_function.22(u, v)* w_function.22(v, xi))
  }
  as.numeric(integral2(integrate_function, 0, 1, 0, 1)$Q)
}


# Set parameters
sample_size<-5*10^4
beta0 <- 1; beta1 <- 20; beta2 <- 3; beta3 <- 2
beta_vec <- matrix(c(beta2, beta3), ncol = 1)
r <- 3
sparsity_levels <- c(-0.75,-0.25)  # Different sparsity levels to iterate over
graphon_size <- 1000  # Fixed graph size for all experiments

# Loop over each sparsity level
for (s in sparsity_levels) {
  
  cat("Running experiments for sparsity level:", s, "\n")
  
  # Initialize containers for OLS statistics
  y2_ds <- vector("numeric", sample_size)
  y2z_ds <- vector("numeric", sample_size)
  z_ds <- vector("numeric", sample_size)
  z_2nd_ds <- vector("numeric", sample_size)
  
  # Repeat the experiment for sample_size times
  for (ii in 1:sample_size) {
    cat("Running simulation repeat =", ii, "\n")
    
    # Update nu_seq function with current sparsity
    my_graphon_new <- new("graphon", w_function = graphon.22, nu_seq = function(x) return(x^(s)))
    my_graphon_simulate <- create_simulation_function(my_graphon_new, return_P_mat = TRUE)
    
    sample_i <- my_graphon_simulate(graphon_size)
    A_mat <- sample_i[['A_mat']]
    P_mat <- sample_i[['P_mat']]
    latent_positions <- sample_i[['latent_positions']]
    
    di <- colSums(A_mat)
    
    # Calculate f(xi) based on latent positions
    f_xi_tri <- sapply(latent_positions, f_xi_func.22_tri)
    
    # Calculate zi_hat based on adjacency matrix
    graph_obj_original <- graph_from_adjacency_matrix(A_mat, mode = "undirected", diag = TRUE)
    nu_obv <- sum(A_mat) / (2 * choose(graphon_size, 2))
    
    # Observed triangle count
    rooted_tri <- count_triangles(graph_obj_original)
    rooted_tri.density <- rooted_tri / (nu_obv^3 * choose(graphon_size - 1, r - 1))
    
    # Generate independent y
    error <- rnorm(graphon_size, 0, 1)
    y.response_original_2 <- beta0 + beta1 * f_xi_tri + error
    
    # Calculate OLS statistics
    y2_ds[ii] <- mean(y.response_original_2)
    y2z_ds[ii] <- mean(y.response_original_2 * rooted_tri.density)
    z_ds[ii] <- mean(rooted_tri.density)
    z_2nd_ds[ii] <- mean(rooted_tri.density^2)
  }
  
  # Compute the mean of the OLS statistics for the current sparsity level
  mean_y2_ds <- mean(y2_ds)
  mean_y2z_ds <- mean(y2z_ds)
  mean_z_ds <- mean(z_ds)
  mean_z_2nd_ds <- mean(z_2nd_ds)
  
  # Store results in a data frame
  results_df <- data.frame(
    Sparsity = s,
    Mean_y2_ds = mean_y2_ds,
    Mean_y2z_ds = mean_y2z_ds,
    Mean_z_ds = mean_z_ds,
    Mean_z_2nd_ds = mean_z_2nd_ds
  )
  
  # Save the results to a CSV file for each sparsity level
  file_name <- paste0("population_mean_results_", s, ".csv")
  write.csv(results_df, file_name, row.names = FALSE)
  cat("Saved results to", file_name, "\n")
}




