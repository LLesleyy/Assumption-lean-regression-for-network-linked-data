

# necessary libraries
library(methods)
library(Rcpp)
library(RcppArmadillo)
library(pracma)
library(igraph)
library(mvtnorm)
library(Matrix)
library(combinat)

set.seed(26)


# note: the probability matrix is held as an input instead of any particular symmetric kernel or w function.
# with a gaussian probability matrix specified separately later.

# Function to simulate from a probability matrix
simulate_from_P_mat = function(P_mat){
  cppFunction( 'arma::mat my_function( int n, NumericMatrix P_mat
  ) {
               arma::mat A_mat(n,n,arma::fill::zeros);
               
               for (int ii = 0; ii < n; ii++)
               {
               for(int jj = ii+1; jj < n; jj++) 
               {
               
               double P_temp = P_mat(ii,jj);
               if( P_temp < 0) {
               P_temp = 0;
               }
               
               A_mat(ii,jj) = R::rbinom(1,P_temp);
               }
               }
               A_mat = arma::symmatu(A_mat);
               return(A_mat);
}
',depends = "RcppArmadillo")
  
  output = my_function(n= nrow(P_mat),P_mat=P_mat)
  return(output)
}


# Function to simulate from a specific Gaussian latent space model
generate_P_mat = function(latent_vec){
  n = length(latent_vec)
  P_mat_out = matrix(0,nrow=n,ncol=n)
  for(ii in 1:(n-1))
    for(jj in 1:(n-ii)) {
      P_mat_out[ii,ii+jj] = exp(-(latent_vec[ii]-latent_vec[ii+jj])^2/10)
    }
  P_mat_out[lower.tri(P_mat_out)] = t(P_mat_out)[lower.tri(P_mat_out)]
  return(P_mat_out)
}

compute_neighborhood_avg <- function(xi_i) {
  # Parameters for the bivariate normal distribution
  # mean=c(1,3,0), sigma = matrix(c(1,0.6 ,0.3, 0.6,4,-0.4,0.3,-0.4,1),nrow=3,ncol=3,byrow=T)
  mu <- c(1, 0)
  sigma <- matrix(c(1, 0.3, 0.3, 1),nrow=2,ncol=2,byrow=T)
  
  # Define the integrand, xi_i is given.
  integrand <- function(xi_j, x_j) {
    exp(-((xi_i - xi_j)^2) / 10) * x_j * dmvnorm(c(xi_j, x_j), mean = mu, sigma = sigma)
  }
  
  # Numerical integration
  result <- integral2(fun = Vectorize(function(xi_j, x_j) integrand(xi_j, x_j)),  # Vectorize the integrand
                      xmin = -6 * sqrt(sigma[1, 1]), xmax = 6 * sqrt(sigma[1, 1]),  # Approximate -Inf, Inf
                      ymin = -6 * sqrt(sigma[2, 2]), ymax = 6 * sqrt(sigma[2, 2])   # Approximate -Inf, Inf
                      )$Q
  return(result)
}



compute_kernel_avg <- function(xi_i) {
  # marginal normal from the joint distribution.
  # mean=c(1,3,0), sigma = matrix(c(1,0.6 ,0.3, 0.6,4,-0.4,0.3,-0.4,1),nrow=3,ncol=3,byrow=T)

  # Define the integrand for the expectation
  integrand <- function(xi_j) {
    exp(-((xi_i - xi_j)^2) / 10) * dnorm(xi_j, mean = 0, sd = 1)
  }
  
  # Compute the integral numerically
  result <- integrate(integrand, lower = 0 - 6 * 1, upper = 0 + 6 * 1)$value
  return(result)
}






# check validity 
# Simulation parameters
num_repeats <- 1;
bootstrap_size <- 2000;

beta0<-1;beta1<-20;beta2<-3;beta3<-2;
beta_vec = matrix(c(beta2, beta3), ncol = 1);

graphon_size <- 3500;r<-3;

sparsity<-c(-0.75,-0.6,-0.55,-0.45,-0.4,-0.25)
k=6;
nu<-graphon.size^(sparsity[k]);


# Containers for storing results
downsampling_results <- vector("list", num_repeats)

# Repeat the entire process for `num_repeats` times
for (repeat_idx in 1:num_repeats) {
  cat("Running repeat", repeat_idx, "of", num_repeats, "\n")
  
  X_vec = rmvnorm(graphon_size,mean=c(1,3,0), sigma = matrix(c(1,0.6 ,0.3, 0.6,4,-0.4,0.3,-0.4,1),nrow=3,ncol=3,byrow=T))
  latent_positions = X_vec[,3]
  X1<-X_vec[,1]
  X2<-X_vec[,2]
  X_covariate<-X_vec[,c(1,2)]

  # One dimension of the multivariate Normal is the latent position
  P_mat_temp = generate_P_mat(latent_positions)
  
  # for an observed graph sample, the link probability should be scaled by sparsity level in advance. the input parameter are graph size specific.
  # Simulate from Gaussian latent space model
  A_mat = simulate_from_P_mat(nu*P_mat_temp)
  P_mat = nu*P_mat_temp
  
  di = colSums(A_mat)
  
  # Calculate f(xi) based on latent positions
  f_xi_Numerator = sapply(latent_positions, compute_neighborhood_avg)
  f_xi_Denominator = sapply(latent_positions, compute_kernel_avg)
  f_xi = f_xi_Numerator/f_xi_Denominator
  
  # Calculate zi_hat based on adjacency matrix
  graph_obj_original = graph_from_adjacency_matrix(A_mat, mode = "undirected", diag = TRUE)
  nu_obv = sum(A_mat) / (2 * choose(graphon_size, 2))
  
  # Compute the neighborhood average
  neighborhood_avg <- (A_mat %*% X1) / rowSums(A_mat)
  # Handle nodes with no neighbors (degree = 0)
  neighborhood_avg[di == 0] <- 0  
  
  # compute for the node-wise adjusted squared terms.
  rooted_ts<-choose(di, 2)
  # matrix multiplication yield neighborhood summation
  Ax <- A_mat %*% X1 
  Ax_squared <- Ax^2
  # Subtract the diagonal contribution (self-loops removed)
  diagonal_contribution <- rowSums((A_mat %*% (X1^2)))
  # Compute the weighted rooted two-star for each node
  weighted_two_star <- Ax_squared - diagonal_contribution
  # corrected node-wise squared measure.
  nei_avg_2nd<-weighted_two_star/rooted_ts
  
  
  # Generate independent y
  error = rnorm(graphon_size, 0, 0.5) 
  # bivariate covariate, local average of one of them is also included.
  y.response_original = beta0 + beta1 * f_xi  + X_covariate %*% beta_vec + error
  
  # Pair network data (y_i, \hat{z_i})
  original_data = as.data.frame(cbind(y.response_original, neighborhood_avg, X_covariate))
  colnames(original_data) = c('y.response_original', 'x1_neighbourhood_ave', 'x1', 'x2')
  
  
  # to check if bias correction works for each node.
  corrected_node_measure<-(original_data$y.response_original*original_data$x1_neighbourhood_ave)/nei_avg_2nd
  noiseless_node_measure<-(original_data$y.response_original*f_xi)/f_xi^2
  non_corrected_measure<-(original_data$y.response_original*original_data$x1_neighbourhood_ave)/(original_data$x1_neighbourhood_ave^2)
  
  hist((corrected_node_measure-noiseless_node_measure))
  hist((non_corrected_measure-noiseless_node_measure))
  
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
  
  z_2nd_ds_adj<- vector("numeric", B)
  
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
    z_2nd_ds[ii]<-mean(original_data[ds_index, 2]^2)
    
    z_2nd_ds_adj[ii] <- mean(nei_avg_2nd[ds_index])
    
    ds_X_XT <- matrix(c(
      1, z_ds[ii], x1_ds[ii], x2_ds[ii],
      z_ds[ii], z_2nd_ds_adj[ii], x1z_ds[ii], x2z_ds[ii],
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
  save(downsampling_results, file = paste0("1215_downsampling_nei_avg_sparsity_", k, ".rds"))
}

cat("Downsampling process completed and results saved.\n")


