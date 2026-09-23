

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
# Rmk: conditional expectation function of the rooted two star
# spread of the kernel
#hist(sapply(runif(1000,0,1),f_xi_func.22_tri))
#hist(sapply(runif(1000,0,1),f_xi_func.22_ts))


# the function for searching all the non-centered two star, 
# get the count weighted by the root node, perform the pre-calculation for each node.
weighted_non_rooted_two_star_count <- function(A, y, x1, x2) {
  # Convert A to sparse matrix
  A_sparse <- Matrix(A, sparse = TRUE)
  
  # Calculate degrees and degrees - 1
  degrees <- rowSums(A_sparse)
  degrees_minus_1 <- degrees - 1
  
  # Precompute weighted adjacency matrices
  weighted_A_y <- A_sparse %*% Diagonal(x = y)
  weighted_A_x1 <- A_sparse %*% Diagonal(x = x1)
  weighted_A_x2 <- A_sparse %*% Diagonal(x = x2)
  
  # Compute weighted counts using matrix-vector multiplications
  scaling_matrix <- matrix(degrees_minus_1, nrow = n, ncol = n, byrow = TRUE)
  weighted_counts_yz <- rowSums(weighted_A_y * scaling_matrix)
  weighted_counts_x1z <- rowSums(weighted_A_x1 * scaling_matrix)
  weighted_counts_x2z <- rowSums(weighted_A_x2 * scaling_matrix)
  
  return(list(weighted_counts_yz = weighted_counts_yz,
              weighted_counts_x1z = weighted_counts_x1z,
              weighted_counts_x2z = weighted_counts_x2z))
}

find_size_4_non_rooted_wheels_count_sparse <- function(A) {
  # Convert A to a sparse matrix
  A_sparse <- Matrix(A, sparse = TRUE)
  
  # Calculate degrees of all nodes
  degrees <- rowSums(A_sparse)
  
  # Precompute valid combinations
  valid_combinations <- ifelse((degrees - 1) >= 3, choose(degrees - 1, 3), 0)
  
  # Compute counts using matrix-vector multiplication
  counts2 <- A_sparse %*% valid_combinations
  
  return(as.numeric(counts2))
}


## simulation settings
#  the current sparsity level 
sparsity<-c(-0.75,-0.6,-0.55,-0.45,-0.4,-0.25)
k=5;
## Parameters for outer repetition
outer_repetitions <- 200  # Number of outer repetitions
bootstrap.size <- 2000    # Number of bootstrap iterations
graphon.size<-4000;
# configuration parameters for the simulation study.
nu<-graphon.size^(sparsity[k]);
r<-3;
beta0<-1;beta1<-20;beta2<-3;beta3<-2;






# Initialize containers to store the results of each repetition
OLS_all_bootstrap_samples <- list()
# now we have the whole beta vector to be inferred.
bootstraped_beta_linear_multiplier<-matrix(0,ncol = 4, nrow=bootstrap.size);

# store the bootstrapped statistics, overall 13 statistics involved to calculate the beta vector.
bootstraped_y<-vector( "numeric" , bootstrap.size);
bootstraped_yx1<-vector( "numeric" , bootstrap.size);
bootstraped_yx2<-vector( "numeric" , bootstrap.size);
bootstraped_yz<-vector( "numeric" , bootstrap.size);

bootstraped_z2<-vector( "numeric" , bootstrap.size);
bootstraped_x1_2nd<-vector( "numeric" , bootstrap.size);
bootstraped_x2_2nd<-vector( "numeric" , bootstrap.size);

bootstraped_z<-vector( "numeric" , bootstrap.size);
bootstraped_x1<-vector( "numeric" , bootstrap.size);
bootstraped_x2<-vector( "numeric" , bootstrap.size);

bootstraped_x1x2<-vector( "numeric" , bootstrap.size);
bootstraped_x1z<-vector( "numeric" , bootstrap.size);
bootstraped_x2z<-vector( "numeric" , bootstrap.size);

bootstraped_nu<-vector( "numeric" , bootstrap.size);


g1_hat_yz <- vector( "numeric" , graphon.size)
g1_hat_x1z <- vector( "numeric" , graphon.size)
g1_hat_x2z <- vector( "numeric" , graphon.size)
g1_hat_z <- vector( "numeric" , graphon.size)

# add an additional part for pre-computation of rho_n.
g1_hat_nu <- vector( "numeric" , graphon.size)
g1_hat_z2 <- vector( "numeric" , graphon.size)

for (outer in 1:outer_repetitions) {
  cat("Starting repetition", outer, "of", outer_repetitions, "\n")
  
  ## observed graph
  # PART 2. Generating Original Dataset for bootstrap 
  # (re-run the same process for testing bootstrap reliability on different original dataset for coverage).
  my_graphon_new=new("graphon",w_function = graphon.22,nu_seq=function(x)return(x^(sparsity[k])))
  my_graphon_simulate = create_simulation_function(my_graphon_new,return_P_mat = TRUE)
  
  # generate one data-set as our original data-set.
  sample_i<-my_graphon_simulate(graphon.size) 
  # store the adjacency matrix.
  A_mat<-sample_i[['A_mat']]
  di=colSums(A_mat)
  P_mat<-sample_i[['P_mat']]
  # store the latent positions.
  latent_positions<-sample_i[['latent_positions']]
  
  X_vec = rmvnorm(graphon.size,mean=c(1,3), sigma = matrix(c(1,0.6,0.6,4),nrow=2,ncol=2,byrow=T))
  beta_vec=matrix(c(beta2,beta3),ncol =1)
  
  # 1. Calculate f(xi) based on latent positions
  f_xi_ts <- sapply(latent_positions, f_xi_func.22_ts)
  
  # 2. calculate zi_hat based on adjacency matrix.
  # create graph from adjacency matrix.
  graph_obj_original=graph_from_adjacency_matrix(A_mat,mode = "undirected",diag = TRUE)
  
  # estimate sparsity level based on the observed graph: number of edges/all possible pairs choose(n,2),
  # we need to adjust for the normalizing constant
  nu_obv<-sum(A_mat)/(2*choose(graphon.size,2))
  
  # # built-in triangle count function.
  # triangle.count_original <- count_triangles(graph_obj_original)
  # sum(triangle.count_original)/3
  # z.triangle.density <- triangle.count/(nu_obv_sample[i]^3*choose(n-1,r-1))
  
  # observed rooted two-star density.
  # start with the true nu.
  rooted_ts<-choose(di,2)
  # rooted_ts.density <- rooted_ts/(nu_obv^2*choose(graphon.size-1,r-1))
  rooted_ts.density <- rooted_ts/(nu^2*choose(graphon.size-1,r-1))
  #summary(rooted_ts.density)
  #hist(rooted_ts.density)
  
  ## 3. generating independent y.
  error<-rnorm(graphon.size,0,1)
  y.response_original<-beta0+beta1*f_xi_ts+X_vec%*%beta_vec+error
  
  # pair network data (y_i,\hat(z_i));
  original_data<-as.data.frame(cbind(y.response_original,rooted_ts.density,X_vec))
  colnames(original_data)<-c('y.response_original','rooted_ts_density','x1','x2')
  
  ## precomputations
  ## Pre-computation of g1_hat for the all the first-order mean statistics involving the network covariate.
  # return a matrix that records pair-wise edge involved triangle counts.
  A1=Matrix(A_mat,sparse=T)
  A2= A1 %*% A1
  diag(A2)=0
  # A2 gives the pair-wise length-2 walk count.
  
  ## 2. calculate \hat_{g(i)} for those first order mean statistics involving network covariate.
  # original_data<-as.data.frame(cbind(y.response_original,z.triangle.density_original,X_vec))
  # prepare the observed mean statistics:
  yz_mean <- mean((original_data$y.response_original)*(original_data$rooted_ts_density))
  x1z_mean <- mean((original_data$rooted_ts_density)*(original_data[,3]))
  x2z_mean <- mean((original_data$rooted_ts_density)*(original_data[,4]))
  z_mean<-mean(original_data$rooted_ts_density)
  
  # pre-calculate the factor.
  factor1<-1/choose(graphon.size-1,r-1);
  
  # calculate the leading term of Hoeffding decomposition for i from 1 to n.
  weighted_nonroot_count<-weighted_non_rooted_two_star_count(A_mat,original_data$y.response_original,
                                                             original_data$x1,original_data$x2)
  
  # Vectorized calculations for g1_hat_yz, g1_hat_x1z, g1_hat_x2z, g1_hat_z, g1_hat_nu
  local_center_contri_yz <- original_data[, 1] * rooted_ts
  local_center_contri_x1z <- original_data[, 3] * rooted_ts
  local_center_contri_x2z <- original_data[, 4] * rooted_ts
  
  other_node_contri_yz <- weighted_nonroot_count[['weighted_counts_yz']]
  other_node_contri_x1z <- weighted_nonroot_count[['weighted_counts_x1z']]
  other_node_contri_x2z <- weighted_nonroot_count[['weighted_counts_x2z']]
  
  # Precompute repeated terms
  common_factor <- factor1 / nu^2
  
  # Calculate all g1_hat terms.
  # rowSums(A2) as the nodewise length-2 walks.
  g1_hat_yz <- common_factor * (local_center_contri_yz + other_node_contri_yz) / 3 - yz_mean
  g1_hat_x1z <- common_factor * (local_center_contri_x1z + other_node_contri_x1z) / 3 - x1z_mean
  g1_hat_x2z <- common_factor * (local_center_contri_x2z + other_node_contri_x2z) / 3 - x2z_mean
  g1_hat_z <- common_factor * (rooted_ts + rowSums(A2)) / 3 - z_mean
  g1_hat_nu <- rowSums(A_mat) / (graphon.size - 1) - nu_obv
  
 
  
  ## Pre-computation of g1_hat for the second statistic.
  # return the rooted size-4 wheel, with the first column as rooted node
  peripheral_counts<-find_size_4_non_rooted_wheels_count_sparse(A_mat)
  
  # calculate the rooted size-4 wheel for the each node as the rooted node.
  rooted_wheel<-choose(di,4)
  
  # squared statistic as the original averaged local squared terms.
  z2_mean<-mean((original_data[,2])^2)
  
  # RMK: there is a constant factor need to be adjusted.
  factor2<-1/choose(graphon.size-1,2*r-1-1);
  g1_hat_z2 <- (factor2 / (5 * nu^4)) * (rooted_wheel + peripheral_counts) - z2_mean
  
  # RMK: be careful about the constant factor/and also whether to include distinct count/isomorphic count.
  ## bootstrapping procedures.
  ## Bootstrapping two statistics using a product gaussian multiplier and calculate bootstrapped beta_hat.
  
  y_mean<-mean(original_data[,1])
  x1_mean<-mean(original_data[,3])
  x2_mean<-mean(original_data[,4])
  
  x1_2nd_mean<-mean(original_data[,3]*original_data[,3])
  x2_2nd_mean<-mean(original_data[,4]*original_data[,4])
  x1x2_mean<-mean(original_data[,3]*original_data[,4])
  yx1_mean<-mean(original_data[,1]*original_data[,3])
  yx2_mean<-mean(original_data[,1]*original_data[,4])
  
  # X1=rnorm(graphon.size,1,sqrt(1/2));X2=rnorm(graphon.size,1,sqrt(1/3));
  # W=X1*X2;
  # mean(W)
  # var(W)
  # mean((W-1)^3)
  
  # an approximated linear multiplier bootstrap. 
  # now in this bootstrap cycle, add a few usual multiplier bootstrap schemes for those mean statistics not involving the network covariate, both the first order and second order mean statistics.
  # as to check the bootstrap samples are comparable to the observed statistics.
  observed_X_y<-matrix(c(y_mean,yz_mean,yx1_mean,yx2_mean),ncol = 1)
  
  observed_X_XT<-matrix(c(1,z_mean,x1_mean,x2_mean,
                          z_mean,z2_mean,x1z_mean,x2z_mean,
                          x1_mean,x1z_mean,x1_2nd_mean,x1x2_mean,
                          x2_mean,x2z_mean,x1x2_mean,x2_2nd_mean),
                        nrow = 4,ncol = 4)
  
  for (ii in 1:bootstrap.size){
    cat("Starting bootstrap", ii, "of", bootstrap.size, "\n")
    X1=rnorm(graphon.size,1,sqrt(1/2));X2=rnorm(graphon.size,1,sqrt(1/3));
    W=X1*X2;
    
    #W=rnorm(graphon.size,0,1)
    bootstraped_nu[ii]<-nu_obv+(2/graphon.size)*sum((W-1)*g1_hat_nu)
    
    ## 1.U-statistic type linear multiplier bootstrap.
    # 1.1 first order ones, all of them having order 3 U-statistic structure.
    # (nu_obv/bootstraped_nu[ii])^3
    # (nu_obv/bootstraped_nu[ii])^6
    bootstraped_yz[ii]<-(yz_mean+(r/graphon.size)*sum((W-1)*g1_hat_yz))
    bootstraped_x1z[ii]<-(x1z_mean+(r/graphon.size)*sum((W-1)*g1_hat_x1z))
    bootstraped_x2z[ii]<-(x2z_mean+(r/graphon.size)*sum((W-1)*g1_hat_x2z))
    bootstraped_z[ii]<-(z_mean+(r/graphon.size)*sum((W-1)*g1_hat_z))
    
    # 1.2 second order one, of order 5 U-statistic structure.
    bootstraped_z2[ii]<-(z2_mean+((2*r-1)/graphon.size)*sum((W-1)*g1_hat_z2))
    
    
    ## 2.usual multiplier bootstrap.
    # 2.1 first order mean statistics.
    bootstraped_y[ii]<-y_mean+mean((W-1)*(original_data[,1]-y_mean))
    bootstraped_x1[ii]<-x1_mean+mean((W-1)*(original_data[,3]-x1_mean))
    bootstraped_x2[ii]<-x2_mean+mean((W-1)*(original_data[,4]-x2_mean))
    # 2.2 second order mean statistics.
    bootstraped_x1_2nd[ii]<-x1_2nd_mean+mean((W-1)*(original_data[,3]*original_data[,3]-x1_2nd_mean))
    bootstraped_x2_2nd[ii]<-x2_2nd_mean+mean((W-1)*(original_data[,4]*original_data[,4]-x2_2nd_mean))
    bootstraped_x1x2[ii]<-x1x2_mean+mean((W-1)*(original_data[,3]*original_data[,4]-x1x2_mean))
    bootstraped_yx1[ii]<-yx1_mean+mean((W-1)*(original_data[,1]*original_data[,3]-yx1_mean))
    bootstraped_yx2[ii]<-yx2_mean+mean((W-1)*(original_data[,1]*original_data[,4]-yx2_mean))
    
    ## align all the bootstrapped mean statistics 
    # squared matrix 4*4
    bootstrap_X_XT<-matrix(c(1,bootstraped_z[ii],bootstraped_x1[ii],bootstraped_x2[ii],
                             bootstraped_z[ii],bootstraped_z2[ii],bootstraped_x1z[ii],bootstraped_x2z[ii],
                             bootstraped_x1[ii],bootstraped_x1z[ii],bootstraped_x1_2nd[ii],bootstraped_x1x2[ii],
                             bootstraped_x2[ii],bootstraped_x2z[ii],bootstraped_x1x2[ii],bootstraped_x2_2nd[ii]),
                           nrow = 4,ncol = 4)
    # column vector 4*1
    bootstrap_X_y<-matrix(c(bootstraped_y[ii],bootstraped_yz[ii],
                            bootstraped_yx1[ii],bootstraped_yx2[ii]),
                          ncol = 1)
    
    # a beta vector.
    bootstraped_beta_linear_multiplier[ii,]<-solve(bootstrap_X_XT)%*%bootstrap_X_y
    
    # test the stability/validity of the delta method.
    # delta_sumtest_boot[ii]<-bootstraped_yz[ii]+bootstraped_z2[ii]
    # delta_protest_boot[ii]<-bootstraped_yz[ii]*bootstraped_z2[ii]
  }
  
  # store the bootstrap samples.
  ## Step 5: Save the bootstrap samples for this repetition
  OLS_all_bootstrap_samples[[outer]] <- bootstraped_beta_linear_multiplier
  
  ## Save the updated results after current repetition
  saveRDS(OLS_all_bootstrap_samples, file = paste0("1213_original_OLS_all_bootstrap_samples_sparsity_", k, ".rds"))
  
}

