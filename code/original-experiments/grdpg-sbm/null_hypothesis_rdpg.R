
library(MASS)
library(abind)

# bootstrap of GRDPG.
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

# size of experiment.
N_graph=3500;

sparsity<-c(-0.4,-0.3,-0.2,-0.1)
k=3;
nu<-N_graph^(sparsity[k])

# Number of blocks
K <- 3


# different SBM test cases
my.sbm_test = new("sbm",B_mat= matrix(c(0.8, 0.2, 0.1,
                                        0.2, 0.7, 0.15,
                                        0.1, 0.15, 0.9), nrow = K, byrow = TRUE), 
                  pi_vec=c(0.65,0.25,0.1),
                  nu_seq = function(x) return(x^(sparsity[k])))

# another test case in from the reference paper.
# my.sbm_test1 = new("sbm",B_mat= matrix(c(0.02,0.03,0.03,0.01),nrow=2), 
#                  pi_vec=c(0.2,0.8),
#                  nu_seq=function(x)return(1))

simulate_sbm_test=create_simulation_function(my.sbm_test,return_P_mat = TRUE)

# observed graph.
sample_graph=simulate_sbm_test(N_graph)
A_obs=sample_graph[['A_mat']]
Z_latent_test=sample_graph[['latent_pos']]
P_generate=sample_graph[['P_mat']]

#############################################################
# construct latent positions:exactly the same way of ASE.   #
#############################################################
# keep all spectral dimensions, no truncation.
# B_mat=matrix(c(0.6,0.2,0.2,0.2),nrow=2)
# B_mat=matrix(c(0.02,0.03,0.03,0.01),nrow=2)
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

# Convert assignment vector into a one-hot matrix.
binary_matrix <- matrix(0, ncol = K, nrow = length(Z_latent_test))
binary_matrix[cbind(seq_along(Z_latent_test), Z_latent_test)] <- 1

# recover the underlying latent position matrix.
X_latent=binary_matrix  %*% U %*%S_sqrt_inv
# for regression, this is the xi_i after normalization that will not degenerate.
# Clean positions use unscaled B; only raw ASE needs rho normalization.
X_latent_scaled <- X_latent

#####################
# construct  ASE.   #
#####################
# generate observed ASE.
eigen_decomposition <- eigen(A_obs)
# Sort eigenvalues and corresponding eigen-vectors.
sorted_indices <- order(abs(eigen_decomposition$values), decreasing = TRUE)
sorted_eigenvalues <- eigen_decomposition$values[sorted_indices]
sorted_eigenvectors <- eigen_decomposition$vectors[, sorted_indices]
# Select the top K eigenvalues and corresponding eigenvectors
top_eigenvalues <- sorted_eigenvalues[1:K]
top_eigenvectors <- sorted_eigenvectors[, 1:K]
# Generate square root of the diagonal matrix
sqrt_lambda <- diag(sqrt(abs(top_eigenvalues)))
# Perform spectral embedding by multiplying eigenvectors with the square root of eigenvalues
# ASE also need to be normalized by sparsity accordingly.

X_estimated <- top_eigenvectors %*% sqrt_lambda
X_estimated_scaled <- nu^(-0.5)*X_estimated

################################################
#### data generating process for regression ####
################################################

beta_coeff <- c(1, 2, 0, 0, 0)

# Generate two weakly correlated bivariate conventional covariates
cov_matrix <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)  # Covariance matrix with weak correlation
bivariate_covariates <- MASS::mvrnorm(N_graph, mu = c(1, 2), Sigma = cov_matrix)  # Generate bivariate data

# Combine the latent features, intercept, and bivariate covariates
X_design_latent <- cbind(bivariate_covariates, X_latent_scaled)
X_design_observed<-cbind(bivariate_covariates, X_estimated_scaled)

y_obs <- X_design_latent %*% matrix(beta_coeff, ncol = 1) + rnorm(N_graph, 0, 1)

beta_obs<-solve(t(X_design_observed)%*%X_design_observed)%*%t(X_design_observed)%*%matrix(y_obs,nrow = N_graph)



######################################################
### bootstrap procedure for the testing statistics ###
######################################################
B_bootstrap <- 1000
beta_boot <- matrix(0, nrow = length(beta_coeff), ncol = B_bootstrap) 

for (B in 1:B_bootstrap) {
  epsilon <- rnorm(N_graph, 0, 1)  # Gaussian multiplier
  
  # Generate two observed statistics
  observed_product <- matrix(0, nrow = N_graph, ncol = length(beta_coeff))
  for (ii in 1:N_graph) {
    observed_product[ii, ] <- X_design_observed[ii, ] * y_obs[ii]  # Entry-wise multiplication
  }
  observed_linear <- colMeans(observed_product)
  observed_quad <- (1 / N_graph) * Reduce(`+`, lapply(1:N_graph, function(i) outer(X_design_observed[i, ], X_design_observed[i, ])))
  
  
  # Perturbation for both linear and quadratic terms
  xi_boot_perturb_linear <- matrix(0, nrow = N_graph, ncol = length(beta_coeff))
  xi_boot_pertub_quad <- vector("list", N_graph)
  
  for (kk in 1:N_graph) {
    xi_obs <- X_design_observed[kk, ]
    xi2_obs_mat <- outer(xi_obs, xi_obs)
    xi_boot_perturb_linear[kk, ] <- epsilon[kk] * (xi_obs * y_obs[kk] - observed_linear)
    xi_boot_pertub_quad[[kk]] <- epsilon[kk] * (xi2_obs_mat - observed_quad)
  }
  
  # Generate current bootstrap sample
  average_x_boot <- observed_linear + colMeans(xi_boot_perturb_linear)
  average_x2_boot <- observed_quad + apply(abind(xi_boot_pertub_quad, along = 3), c(1, 2), mean)
  
  # Calculate bootstrapped beta coefficients
  beta_boot[,B] <- solve(average_x2_boot) %*% matrix(average_x_boot, ncol = 1)
}


Tn_boot_H0=numeric(B_bootstrap)
for (B in 1:B_bootstrap){
  # only pick the network part coefficients.
  Tn_boot_H0[B]= N_graph*sum((beta_boot[3:5,B]-beta_obs[3:5])^2)
}
Tn_obs_H0=N_graph*sum(beta_obs[3:5]^2)

# Upper-tail bootstrap p-value for the squared-norm network-effect statistic.
p_value <- (1 + sum(Tn_boot_H0 >= Tn_obs_H0)) / (1 + length(Tn_boot_H0))

# Print the p-value
print(p_value)
