
# necessary libraries
library(methods)
library(Rcpp)
library(RcppArmadillo)
library(pracma)
library(igraph)
library(mvtnorm)
library(Matrix)
library(combinat)

library(ggplot2)
library(gridExtra)
# a unified framework of bootstrap for corrected and non-corrected OLS together.
# leverage the same bootstrap scheme.

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

# graph kernel
graphon.22 = 'double w_function(double u, double v) {
return(0.5/0.5*(u+v));
}'

w_function.22 <- function(u, v) {
  return(0.5/0.5*(u+v))
}


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

# parameters
sparsity_levels<-c(-0.75,-0.7,-0.6,-0.4,-0.3,-0.25);
outer_repetitions <- 200  # Number of outer repetitions
bootstrap.size <- 2000    # Number of bootstrap iterations
graphon.size<-4500;
k=6;
r<-3;



#######################################
#      investigate associations       #
#######################################

# estimated rho_n is used in this pipeline.
  
  current_sparsity <- sparsity_levels[k]
  nu <- graphon.size ^ current_sparsity
  
  ## observed graph
  my_graphon_new <- new("graphon", w_function = graphon.22, nu_seq = function(x) return(x^current_sparsity))
  #my_graphon_new=new("graphon",w_function = graphon.22,nu_seq=function(x)return(x^(sparsity[k])))
  my_graphon_simulate = create_simulation_function(my_graphon_new,return_P_mat = TRUE)
  sample_i<-my_graphon_simulate(graphon.size) 
  
  A_mat<-sample_i[['A_mat']]
  di=colSums(A_mat)
  P_mat<-sample_i[['P_mat']]
  latent_positions<-sample_i[['latent_positions']]
  
  # 1. Calculate f(xi) based on latent positions
  f_xi_ts <- sapply(latent_positions, f_xi_func.22_ts)
  
  # 2. calculate zi_hat based on adjacency matrix.
  nu_obv<-sum(A_mat)/(2*choose(graphon.size,2))
  rooted_ts<-choose(di,2)
  rooted_ts.density <- rooted_ts/(nu_obv^2*choose(graphon.size-1,r-1))
  #rooted_ts.density <- rooted_ts/(nu^2*choose(graphon.size-1,r-1))
  
  # 3. Generate outcome Y under complicated nonlinear structure
  X_vec = rmvnorm(graphon.size,mean=c(1,3), sigma = matrix(c(1,0.6,0.6,4),nrow=2,ncol=2,byrow=T))
  beta_vec=matrix(c(beta2,beta3),ncol =1)
  #y_response<-beta0+beta1*f_xi_ts+X_vec%*%beta_vec+error
  
  # generate nonlinear response, with network part signal amplified.
  y_response <- log(1 + 5 * f_xi_ts * abs(X_vec[,1])) + 
                sin(0.5 * X_vec[,2]) * sqrt(5 * f_xi_ts) + rnorm(graphon.size, 0, 1)

  # REMOVE noise to highlight the structure
  y_response_clean <- log(1 + 5 * f_xi_ts * abs(X_vec[,1])) + 
                      sin(0.5 * X_vec[,2]) * sqrt(5 * f_xi_ts)
  
  # scatter plots.
  # Correct original data
  original_data <- data.frame(
    y_response = y_response,
    rooted_ts_density = rooted_ts.density,
    x1 = X_vec[,1],
    x2 = X_vec[,2]
  )
  
  data_plot_noisy <- data.frame(
    Y = y_response_clean,
    rooted_ts_density = rooted_ts.density
  )
  
  # For clean plots (latent signal f(xi))
  data_plot_clean <- data.frame(
    Y = y_response_clean,
    f_xi = f_xi_ts,
    X1 = X_vec[,1],
    X2 = X_vec[,2]
  )
  
  # Common Biometrika-like theme
  bka_theme <- theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, hjust = 0.5),   # no bold
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      axis.ticks = element_line(color = "black", size = 0.3),
      panel.border = element_rect(color = "black", fill = NA, size = 0.3)
    )
  
  # === Plots ===
  p1_clean <- ggplot(data_plot_clean, aes(x = f_xi, y = Y)) +
    geom_point(color = "steelblue4", alpha = 0.5, size = 0.8) +
    labs(x = expression(Z[i]), y = "Y", title = "Y vs True Network Covariate") +
    bka_theme
  
  p2_clean <- ggplot(data_plot_clean, aes(x = X1, y = Y)) +
    geom_point(color = "darkorange3", alpha = 0.5, size = 0.8) +
    labs(x = expression(X[i1]), y = "Y", title = "Y vs X1") +
    bka_theme
  
  p3_clean <- ggplot(data_plot_clean, aes(x = X2, y = Y)) +
    geom_point(color = "seagreen4", alpha = 0.5, size = 0.8) +
    labs(x = expression(X[i2]), y = "Y", title = "Y vs X2") +
    bka_theme
  
  p4_noisy <- ggplot(data_plot_noisy, aes(x = rooted_ts_density, y = Y)) +
    geom_point(color = "orchid4", alpha = 0.5, size = 0.8) +
    labs(x = expression(hat(Z)[i]), y = "Y", title = "Y vs Noisy Network Covariate") +
    bka_theme
  
  # === Panel ===
  final_panel_plot <- grid.arrange(p1_clean, p2_clean, p3_clean, p4_noisy, nrow = 2)
  print(final_panel_plot)
  
  ########################################################
  #      Approximate the clean projection parameter       #
  ########################################################

  # 1. Set parameters
  set.seed(123)
  N_sim <- 10^5 # Monte Carlo size
  
  # Storage for nodewise samples
  f_xi_vals_all <- numeric(N_sim)
  X1_all <- numeric(N_sim)
  X2_all <- numeric(N_sim)
  Y_all  <- numeric(N_sim)
  
  # === MC Simulation ===
  for (i in 1:N_sim) {
    if (i %% 10000 == 0) {
      cat("Progress:", i, "/", N_sim, "\n")
    }
    
    # Sample latent position
    xi_i <- runif(1)
    f_xi_i <- f_xi_func.22_ts(xi_i)
    f_xi_vals_all[i] <- f_xi_i
    
    # Sample conventional covariates X1, X2
    X_i <- rmvnorm(1, mean = c(1, 3), sigma = matrix(c(1, 0.6, 0.6, 4), 2))
    X1_all[i] <- X_i[1]
    X2_all[i] <- X_i[2]
    
    # Generate response Y
    Y_all[i] <- log(1 + 5 * f_xi_i * abs(X_i[1])) +
                sin(0.5 * X_i[2]) * sqrt(5 * f_xi_i) +
                rnorm(1, 0, 1)  # noise
  }
  
  # === Compute population moments ===
  y_mean      <- mean(Y_all)
  z_mean      <- mean(f_xi_vals_all)
  x1_mean     <- mean(X1_all)
  x2_mean     <- mean(X2_all)
  
  z2_mean     <- mean(f_xi_vals_all^2)
  x1_2nd_mean <- mean(X1_all^2)
  x2_2nd_mean <- mean(X2_all^2)
  x1x2_mean   <- mean(X1_all * X2_all)
  
  yz_mean     <- mean(Y_all * f_xi_vals_all)
  yx1_mean    <- mean(Y_all * X1_all)
  yx2_mean    <- mean(Y_all * X2_all)
  x1z_mean    <- mean(X1_all * f_xi_vals_all)
  x2z_mean    <- mean(X2_all * f_xi_vals_all)
  
  # === Arrange into OLS system ===
  observed_X_y <- matrix(c(y_mean, yz_mean, yx1_mean, yx2_mean), ncol = 1)
  
  observed_X_XT <- matrix(c(
    1,      z_mean,     x1_mean,     x2_mean,
    z_mean, z2_mean,    x1z_mean,    x2z_mean,
    x1_mean, x1z_mean,  x1_2nd_mean, x1x2_mean,
    x2_mean, x2z_mean,  x1x2_mean,   x2_2nd_mean
  ), nrow = 4, byrow = TRUE)
  
  # Solve for population beta
  population_beta <- solve(observed_X_XT, observed_X_y)
  rownames(population_beta) <- c("Intercept", "Z (f(xi))", "X1", "X2")
  colnames(population_beta) <- "Population Beta"
  
  # === Output ===
  cat("Estimated Population Beta Coefficients:\n")
  print(round(population_beta, 4))
  
  # === Save result if needed ===
  saveRDS(population_beta, file = "population_beta_projection_estimate.rds")
  
  
  
  ########################################################
  #      Approximate the noisy projection parameter      #
  ########################################################
  
  # Parameters
  set.seed(123)
  N_sim <- 10^5   # Monte Carlo size
  sparsity_levels <- c(-0.75, -0.7, -0.6, -0.4, -0.3, -0.25)
  graphon.size <- 4500
  r <- 3  # two-star rooted size
  
  # Storage for results
  population_beta_noisy_list <- list()
  
  # Define w function
  w_function.22 <- function(u, v) {
    return((u + v))
  }
  
  # Define f_xi for rooted 2-star
  f_xi_func.22_ts <- function(xi) {
    integrate_function <- function(u, v) {
      return(w_function.22(xi, u) * w_function.22(xi, v))
    }
    as.numeric(integral2(integrate_function, 0, 1, 0, 1)$Q)
  }
  
  # Iterate over sparsity levels
  for (k in seq_along(sparsity_levels)) {
    
    current_sparsity <- sparsity_levels[k]
    cat("==== Starting simulation for sparsity =", current_sparsity, "====\n")
    
    nu <- graphon.size^current_sparsity
    my_graphon_new <- new("graphon", w_function = graphon.22, nu_seq = function(x) return(x^current_sparsity))
    
    # Storage for each sparsity, we only keep samples for node 1 to approximate the population moment.
    rooted_ts_density_all <- numeric(N_sim)
    X1_all <- numeric(N_sim)
    X2_all <- numeric(N_sim)
    Y_all <- numeric(N_sim)
    
    for (i in 1:N_sim) {
      
      cat("current MC progress:", i, "/", N_sim, "\n")
      
      # simulation takes longer time.
      #if (i %% 1000 == 0) {
      #  cat("Progress for sparsity", current_sparsity, ":", i, "/", N_sim, "\n")
      #}
      
      ## observed graph
      my_graphon_simulate = create_simulation_function(my_graphon_new,return_P_mat = TRUE)
      sample_i<-my_graphon_simulate(graphon.size) 
      
      A_mat<-sample_i[['A_mat']]
      di=colSums(A_mat)
      P_mat<-sample_i[['P_mat']]
      latent_positions<-sample_i[['latent_positions']]
      nu_obv<-sum(A_mat)/(2*choose(graphon.size,2))
      
      # Latent position
      xi_i <- latent_positions[1]
      # two-star calculation
      f_xi_i <- f_xi_func.22_ts(xi_i)
      
      # Observed rooted two-star density (noisy Z_i)
      di_i=di[1]
      rooted_ts_density_i <- choose(di_i,2) / (nu_obv^2 * choose(graphon.size - 1, r - 1))
      rooted_ts_density_all[i] <- rooted_ts_density_i
      
      # Conventional covariates
      X_i <- rmvnorm(1, mean = c(1, 3), sigma = matrix(c(1, 0.6, 0.6, 4), 2))
      X1_all[i] <- X_i[1]
      X2_all[i] <- X_i[2]
      
      # Generate Y-with clean Z_i.
      Y_all[i] <- log(1 + 5 * f_xi_i * abs(X_i[1])) +
                  sin(0.5 * X_i[2]) * sqrt(5 * f_xi_i) +
                  rnorm(1, 0, 1)
    }
    
    # === Compute approximate population moments ===
    y_mean <- mean(Y_all)
    z_mean <- mean(rooted_ts_density_all)
    x1_mean <- mean(X1_all)
    x2_mean <- mean(X2_all)
    
    z2_mean <- mean(rooted_ts_density_all^2)
    x1_2nd_mean <- mean(X1_all^2)
    x2_2nd_mean <- mean(X2_all^2)
    x1x2_mean <- mean(X1_all * X2_all)
    
    yz_mean <- mean(Y_all * rooted_ts_density_all)
    yx1_mean <- mean(Y_all * X1_all)
    yx2_mean <- mean(Y_all * X2_all)
    x1z_mean <- mean(X1_all * rooted_ts_density_all)
    x2z_mean <- mean(X2_all * rooted_ts_density_all)
    
    # === Arrange into OLS system ===
    observed_X_y <- matrix(c(y_mean, yz_mean, yx1_mean, yx2_mean), ncol = 1)
    
    observed_X_XT <- matrix(c(
      1,      z_mean,     x1_mean,     x2_mean,
      z_mean, z2_mean,    x1z_mean,    x2z_mean,
      x1_mean, x1z_mean,  x1_2nd_mean, x1x2_mean,
      x2_mean, x2z_mean,  x1x2_mean,   x2_2nd_mean
    ), nrow = 4, byrow = TRUE)
    
    # Solve for population beta
    population_beta <- solve(observed_X_XT, observed_X_y)
    rownames(population_beta) <- c("Intercept", "Z (Rooted 2-star)", "X1", "X2")
    colnames(population_beta) <- paste0("Sparsity_", current_sparsity)
    
    # Save
    population_beta_noisy_list[[k]] <- population_beta
    
    # Save intermediate result
    saveRDS(population_beta_noisy_list, file = "population_beta_projection_estimates_noisy_rooted_ts.rds")
  }
  
  # === After loop, output all ===
  cat("Finished All Sparsity Levels!\n")
  print(population_beta_noisy_list)
  
  ########################################################
  #      Approximate the noisy projection parameter      #
  #      faster version with true sparsity parameter     #
  ########################################################
  
  # Libraries
  library(mvtnorm)
  library(foreach)
  library(doParallel)
  library(parallel)
  
  # Set global parameters
  set.seed(123)
  N_sim <- 10^5  # Monte Carlo iterations per sparsity
  graphon.size <- 4500
  r <- 3
  sparsity_levels <- c(-0.75, -0.7, -0.6, -0.4, -0.3, -0.25)
  
  # Define graphon kernel
  w_function.22 <- function(u, v) {
    return(u + v)
  }
  
  # Define f_xi for rooted 2-star (population version)
  f_xi_func.22_ts <- function(xi) {
    integrate_function <- function(u, v) {
      return(w_function.22(xi, u) * w_function.22(xi, v))
    }
    as.numeric(integral2(integrate_function, 0, 1, 0, 1)$Q)
  }
  
  # Prepare for parallel computation
  n_cores <- detectCores() - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Storage
  population_beta_noisy_list <- list()
  
  # Iterate over sparsity levels
  for (k in seq_along(sparsity_levels)) {
    
    current_sparsity <- sparsity_levels[k]
    cat("==== Starting simulation for sparsity =", current_sparsity, "====\n")
    
    nu <- graphon.size^current_sparsity
    
    # Monte Carlo simulation with parallel backend
    result <- foreach(i = 1:N_sim, .combine = rbind, .packages = c('mvtnorm', 'pracma')) %dopar% {
      
      # 1. Sample all latent positions
      xi_all <- runif(graphon.size)
      
      # 2. Compute connection probabilities for node 1
      probs <- pmin(nu * w_function.22(xi_all[1], xi_all[-1]), 1)
      
      # 3. Simulate Bernoulli connections for node 1
      A_1j <- rbinom(graphon.size - 1, 1, probs)
      
      # 4. Degree of node 1
      d_1 <- sum(A_1j)
      
      # 5. Compute observed rooted 2-star density
      rooted_2star_density <- if (d_1 >= 2) {
        choose(d_1, 2) / (nu^2 * choose(graphon.size - 1, r - 1))
      } else {
        0
      }
      
      # 6. Conventional covariates
      X_i <- rmvnorm(1, mean = c(1, 3), sigma = matrix(c(1, 0.6, 0.6, 4), 2))
      
      # 7. True latent signal (for Y generation)
      f_xi_1 <- f_xi_func.22_ts(xi_all[1])
      
      # 8. Generate Y using f_xi_1
      Y_i <- log(1 + 5 * f_xi_1 * abs(X_i[1])) +
             sin(0.5 * X_i[2]) * sqrt(5 * f_xi_1) +
             rnorm(1, 0, 1)
      
      # 9. return list for monte-carlo iterations.
      c(Y_i, rooted_2star_density, X_i[1], X_i[2])
      
    }
    
    # Result columns
    colnames(result) <- c("Y", "Z_hat", "X1", "X2")
    
    # Extract components
    Y_all <- result[, "Y"]
    Z_hat_all <- result[, "Z_hat"]
    X1_all <- result[, "X1"]
    X2_all <- result[, "X2"]
    
    # === Compute approximate population moments ===
    y_mean      <- mean(Y_all)
    z_mean      <- mean(Z_hat_all)
    x1_mean     <- mean(X1_all)
    x2_mean     <- mean(X2_all)
    
    z2_mean     <- mean(Z_hat_all^2)
    x1_2nd_mean <- mean(X1_all^2)
    x2_2nd_mean <- mean(X2_all^2)
    x1x2_mean   <- mean(X1_all * X2_all)
    
    yz_mean     <- mean(Y_all * Z_hat_all)
    yx1_mean    <- mean(Y_all * X1_all)
    yx2_mean    <- mean(Y_all * X2_all)
    x1z_mean    <- mean(X1_all * Z_hat_all)
    x2z_mean    <- mean(X2_all * Z_hat_all)
    
    # === Set up OLS ===
    observed_X_y <- matrix(c(y_mean, yz_mean, yx1_mean, yx2_mean), ncol = 1)
    
    observed_X_XT <- matrix(c(
      1,        z_mean,     x1_mean,     x2_mean,
      z_mean,   z2_mean,    x1z_mean,    x2z_mean,
      x1_mean,  x1z_mean,   x1_2nd_mean, x1x2_mean,
      x2_mean,  x2z_mean,   x1x2_mean,   x2_2nd_mean
    ), nrow = 4, byrow = TRUE)
    
    # Estimate projection beta
    population_beta <- solve(observed_X_XT, observed_X_y)
    rownames(population_beta) <- c("Intercept", "Z_hat", "X1", "X2")
    colnames(population_beta) <- paste0("Sparsity_", current_sparsity)
    
    # Save results
    population_beta_noisy_list[[k]] <- population_beta
    cat("==== Completed sparsity =", current_sparsity, "====\n")
    
    # Save after each sparsity level to avoid data loss
    saveRDS(population_beta_noisy_list, file = "population_beta_projection_estimates_noisy_rooted_ts_PARALLEL_FINAL.rds")
  }
  
  # Stop parallel cluster
  stopCluster(cl)
  
  # Output final results
  cat("==== Finished all sparsity levels! ====\n")
  print(population_beta_noisy_list)
  
  
  
  