

# necessary libraries
library(methods)
library(Rcpp)
library(RcppArmadillo)
library(pracma)
library(igraph)
library(mvtnorm)
library(Matrix)
library(combinat)

## algorithm for search the Rooted two-star/size-4-wheel and Peripheral two-star/size-4-wheel
# Set object classes for SBM, Graphon, LSM
setClass("graphon",representation(w_function="character",nu_seq = "function"))
create_simulation_function =  function(x,return_P_mat =FALSE) UseMethod("create_simulation_function")

create_simulation_function.graphon = function(graphon,return_P_mat =FALSE) {
  
  cppFunction( 'arma::mat my_function_1( int n, double nu) {
               arma::mat A_mat(n,n,arma::fill::zeros);
               NumericVector latent_positions = runif(n,0,1);
               for (int ii = 0; ii < n; ii++)
               {
               for(int jj = ii+1; jj< n; jj++) 
               {
               double P_temp = nu * w_function(latent_positions[ii],latent_positions[jj]);
               A_mat(ii,jj) =  R::rbinom(1,P_temp);
               }
               }
               A_mat = arma::symmatu(A_mat);
               return(A_mat);
}
',depends = "RcppArmadillo", includes = graphon@w_function)
  
  cppFunction( 'List my_function_2( int n, double nu) {
               arma::mat A_mat(n,n,arma::fill::zeros);
               arma::mat P_mat(n,n,arma::fill::zeros);
               NumericVector latent_positions = runif(n,0,1);
               for (int ii = 0; ii < n; ii++)
               {
               for(int jj = ii+1; jj< n; jj++) 
               {
               double P_temp = nu * w_function(latent_positions[ii],latent_positions[jj]);
               A_mat(ii,jj) =  R::rbinom(1,P_temp);
               P_mat(ii,jj) = P_temp;
               }
               }
               A_mat = arma::symmatu(A_mat);
               P_mat = arma::symmatu(P_mat);
               List ret;
               ret["A_mat"] = A_mat;
               ret["P_mat"] = P_mat;
               ret["latent_positions"] = latent_positions;
               return ret;
}
',depends = "RcppArmadillo", includes = graphon@w_function) 
  
  
  if(return_P_mat){  
    output = function(n) {
      return( my_function_2(n=n,nu=graphon@nu_seq(n)))
    }
    return(output)
  }
  else{
    output = function(n) {
      return( my_function_1(n=n,nu=graphon@nu_seq(n)))
      
      return(output)
    }
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

f_xi_func.22_ts <- function(xi) {
  integrate_function <- function(u, v) {
    return(w_function.22(xi, u) * w_function.22(xi, v))
  }
  as.numeric(integral2(integrate_function, 0, 1, 0, 1)$Q)
}



# Simulation parameters
num_repeats <- 200;
bootstrap_size <- 2000;

beta0<-1;beta1<-20;beta2<-3;beta3<-2;
beta_vec = matrix(c(beta2, beta3), ncol = 1);

graphon_size <- 4000;r<-3;

sparsity<-c(-0.75,-0.6,-0.55,-0.45,-0.4,-0.25)
k=6;
nu<-graphon.size^(sparsity[k]);


# Containers for storing results
downsampling_results <- vector("list", num_repeats)

# Repeat the entire process for `num_repeats` times
for (repeat_idx in 1:num_repeats) {
  cat("Running repeat", repeat_idx, "of", num_repeats, "\n")
  
  ## Simulation settings
  my_graphon_new = new("graphon", w_function = graphon.22, nu_seq = function(x) return(x^(sparsity[k])))
  my_graphon_simulate = create_simulation_function(my_graphon_new, return_P_mat = TRUE)
  
  sample_i = my_graphon_simulate(graphon_size)
  A_mat = sample_i[['A_mat']]
  P_mat = sample_i[['P_mat']]
  latent_positions = sample_i[['latent_positions']]
  
  di = colSums(A_mat)
  
  X_vec = rmvnorm(graphon_size, mean = c(1, 3), sigma = matrix(c(1, 0.6, 0.6, 4), nrow = 2, ncol = 2, byrow = T))
  
  # Calculate f(xi) based on latent positions
  ts_results <- sapply(latent_positions, f_xi_func.22_ts)
  tri_results <- sapply(latent_positions, f_xi_func.22_tri)
  
  # Perform element-wise division
  f_xi <- ts_results / tri_results

  # Calculate zi_hat based on adjacency matrix
  graph_obj_original = graph_from_adjacency_matrix(A_mat, mode = "undirected", diag = TRUE)
  nu_obv = sum(A_mat) / (2 * choose(graphon_size, 2))
  
  # Observed rooted two-star count
  # try with true nu first.
  rooted_ts = choose(di, 2)
  rooted_ts.density = rooted_ts / (nu^2 * choose(graphon_size - 1, r - 1))

  rooted_tri = count_triangles(graph_obj_original)
  rooted_tri.density <- rooted_tri/(nu^3*choose(graphon.size-1,r-1))
  rooted_trans <- rooted_ts.density/rooted_tri.density
  
  # need to adjust for the nodewise squared measure
  # which will be rooted size-4 wheel/rooted bi-triangles
  ts_2nd_corrected <- choose(di, 4)/(nu^4 * choose(graphon.size, 2 * r - 2))
  tri_2nd_correted <- 
  
  rooted_trans_2nd_corrected <-
  
  # Generate independent y
  error = rnorm(graphon_size, 0, 0.5)
  y.response_original = beta0 + beta1 * f_xi + X_vec %*% beta_vec + error
  
  # Pair network data (y_i, \hat{z_i})
  original_data = as.data.frame(cbind(y.response_original, rooted_trans, X_vec))
  colnames(original_data) = c('y.response_original', 'local_trans_density', 'x1', 'x2')
  
  ## downsampling for the linear statistics and squared statistics
  # set down sampling size
  m = 0.7*graphon_size
  
  #m2=graphon_size
  B = bootstrap_size
  
  z_ds <- vector("numeric", B)
  x1_ds <- vector("numeric", B)
  x2_ds <- vector("numeric", B)
  y_ds <- vector("numeric", B)
  
  yz_ds <- vector("numeric", B)
  yx1_ds <- vector("numeric", B)
  yx2_ds <- vector("numeric", B)
  x1z_ds <- vector("numeric", B)
  x2z_ds <- vector("numeric", B)
  x1x2_ds <- vector("numeric", B)
  
  x1_2nd_ds <- vector("numeric", B)
  x2_2nd_ds <- vector("numeric", B)
  z_2nd_ds <- vector("numeric", B)
  
  beta_ds <- matrix(0, nrow = B, ncol = 4)
  
  for (ii in 1:B) {
    # down-sampling index
    ds_index <- sample(1:graphon_size, m)
    
    z_ds[ii] <- mean(original_data[ds_index, 2])
    x1_ds[ii] <- mean(original_data[ds_index, 3])
    x2_ds[ii] <- mean(original_data[ds_index, 4])
    y_ds[ii] <- mean(original_data[ds_index, 1])
    
    yz_ds[ii] <- mean(original_data[ds_index, 1] * original_data[ds_index, 2])
    yx1_ds[ii] <- mean(original_data[ds_index, 1] * original_data[ds_index, 3])
    yx2_ds[ii] <- mean(original_data[ds_index, 1] * original_data[ds_index, 4])
    x1z_ds[ii] <- mean(original_data[ds_index, 3] * original_data[ds_index, 2])
    x2z_ds[ii] <- mean(original_data[ds_index, 4] * original_data[ds_index, 2])
    x1x2_ds[ii] <- mean(original_data[ds_index, 3] * original_data[ds_index, 4])
    
    x1_2nd_ds[ii] <- mean(original_data[ds_index, 3]^2)
    x2_2nd_ds[ii] <- mean(original_data[ds_index, 4]^2)
    z_2nd_ds[ii] <- mean(original_data[ds_index, 2]^2)
    
    ds_X_XT <- matrix(c(
      1, z_ds[ii], x1_ds[ii], x2_ds[ii],
      z_ds[ii], z_2nd_ds[ii], x1z_ds[ii], x2z_ds[ii],
      x1_ds[ii], x1z_ds[ii], x1_2nd_ds[ii], x1x2_ds[ii],
      x2_ds[ii], x2z_ds[ii], x1x2_ds[ii], x2_2nd_ds[ii]
    ), nrow = 4, ncol = 4)
    
    # column vector 4*1
    ds_X_y <- matrix(c(y_ds[ii], yz_ds[ii], yx1_ds[ii], yx2_ds[ii]), ncol = 1)
    # a beta vector
    beta_ds[ii, ] <- solve(ds_X_XT) %*% ds_X_y
  }
  
  # Store results for this repeat
  downsampling_results[[repeat_idx]] <- beta_ds
  
  # Save the results to a file
  save(downsampling_results, file = paste0("1215_downsampling_sparsity_", k, ".rds"))
  
}

cat("Downsampling process completed and results saved.\n")
