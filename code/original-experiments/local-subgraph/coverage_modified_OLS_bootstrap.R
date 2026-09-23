# Load the required libraries
library(methods)
library(Rcpp)
library(RcppArmadillo)
library(pracma)
library(igraph)
library(mvtnorm)
library(Matrix)
library(combinat)

# True regression coefficients
true_beta <- c(1, 20, 3, 2)

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args)) args[1L] else "results"

# Load the saved bootstrap results for the dense case:n^(-0.6)
loaded_data<-readRDS(file.path(results_dir, "all_bootstrap_samples_sparsity_2.rds"))
# modified case, a dense regime

# Number of repetitions
num_repeats <- length(loaded_data)

# Initialize matrix to store CI coverage results
ci_coverage <- matrix(0, nrow = num_repeats, ncol = length(true_beta))

# Calculate CI coverage for each repeat
for (repeat_idx in 1:num_repeats) {
  bootstraped_beta_linear_multiplier <- loaded_data[[repeat_idx]]
  
  # Compute confidence intervals for each beta
  ci_lower <- apply(bootstraped_beta_linear_multiplier, 2, quantile, probs = 0.025)
  ci_upper <- apply(bootstraped_beta_linear_multiplier, 2, quantile, probs = 0.975)
  
  # Check if the true beta coefficients fall within the CI
  ci_coverage[repeat_idx, ] <- (true_beta >= ci_lower & true_beta <= ci_upper)
}

# Calculate CI coverage rate for each coefficient
ci_coverage_rate <- colMeans(ci_coverage)
names(ci_coverage_rate) <- paste0("Beta", 0:(length(true_beta) - 1))

# Print CI coverage results.
for (i in seq_along(ci_coverage_rate)) {
  print(paste0("Coverage for the dense case n^(-0.6), ", names(ci_coverage_rate)[i], ": ", ci_coverage_rate[i]))
}


###################################################################################
# Load the saved bootstrap results for the dense case:n^(-0.4)
loaded_data2<-readRDS(file.path(results_dir, "all_bootstrap_samples_sparsity_5.rds"))

# Number of repetitions
num_repeats <- length(loaded_data2)

# Initialize matrix to store CI coverage results
ci_coverage <- matrix(0, nrow = num_repeats, ncol = length(true_beta))

# Calculate CI coverage for each repeat
for (repeat_idx in 1:num_repeats) {
  bootstraped_beta_linear_multiplier <- loaded_data2[[repeat_idx]]
  
  # Compute confidence intervals for each beta
  ci_lower <- apply(bootstraped_beta_linear_multiplier, 2, quantile, probs = 0.025)
  ci_upper <- apply(bootstraped_beta_linear_multiplier, 2, quantile, probs = 0.975)
  
  # Check if the true beta coefficients fall within the CI
  ci_coverage[repeat_idx, ] <- (true_beta >= ci_lower & true_beta <= ci_upper)
}

# Calculate CI coverage rate for each coefficient
ci_coverage_rate <- colMeans(ci_coverage)
names(ci_coverage_rate) <- paste0("Beta", 0:(length(true_beta) - 1))

# Print CI coverage results with descriptive text
for (i in seq_along(ci_coverage_rate)) {
  print(paste0("Coverage for the dense case n^(-0.4), ", names(ci_coverage_rate)[i], ": ", ci_coverage_rate[i]))
}



###################################################################################
# Load the saved bootstrap results for the dense case:n^(-0.25)
loaded_data2<-readRDS(file.path(results_dir, "all_bootstrap_samples_sparsity_6.rds"))
# with n=3500

# Number of repetitions
num_repeats <- length(loaded_data2)

# Initialize matrix to store CI coverage results
ci_coverage <- matrix(0, nrow = num_repeats, ncol = length(true_beta))

# Calculate CI coverage for each repeat
for (repeat_idx in 1:num_repeats) {
  bootstraped_beta_linear_multiplier <- loaded_data2[[repeat_idx]]
  
  # Compute confidence intervals for each beta
  ci_lower <- apply(bootstraped_beta_linear_multiplier, 2, quantile, probs = 0.025)
  ci_upper <- apply(bootstraped_beta_linear_multiplier, 2, quantile, probs = 0.975)
  
  # Check if the true beta coefficients fall within the CI
  ci_coverage[repeat_idx, ] <- (true_beta >= ci_lower & true_beta <= ci_upper)
}

# Calculate CI coverage rate for each coefficient
ci_coverage_rate <- colMeans(ci_coverage)
names(ci_coverage_rate) <- paste0("Beta", 0:(length(true_beta) - 1))

# Print CI coverage results with descriptive text
for (i in seq_along(ci_coverage_rate)) {
  print(paste0("Coverage for the dense case n^(-0.25), ", names(ci_coverage_rate)[i], ": ", ci_coverage_rate[i]))
}


###################################################################################
# Load the saved bootstrap results for the dense case:n^(-0.25)
loaded_data_1213 <- readRDS(file.path(results_dir, "1213_all_bootstrap_samples_sparsity_6.rds"))
# with n=4000, 200:2000(4000)

# Number of repetitions
num_repeats <- length(loaded_data_1213)

# Initialize matrix to store CI coverage results
ci_coverage <- matrix(0, nrow = num_repeats, ncol = length(true_beta))

# Calculate CI coverage for each repeat
for (repeat_idx in 1:num_repeats) {
  bootstraped_beta_linear_multiplier <- loaded_data_1213[[repeat_idx]]
  
  # Compute confidence intervals for each beta
  ci_lower <- apply(bootstraped_beta_linear_multiplier, 2, quantile, probs = 0.025)
  ci_upper <- apply(bootstraped_beta_linear_multiplier, 2, quantile, probs = 0.975)
  
  # Check if the true beta coefficients fall within the CI
  ci_coverage[repeat_idx, ] <- (true_beta >= ci_lower & true_beta <= ci_upper)
}

# Calculate CI coverage rate for each coefficient
ci_coverage_rate <- colMeans(ci_coverage)
names(ci_coverage_rate) <- paste0("Beta", 0:(length(true_beta) - 1))

# Print CI coverage results with descriptive text and formatted percentages
for (i in seq_along(ci_coverage_rate)) {
  coverage_percentage <- sprintf("%.2f", ci_coverage_rate[i] * 100)
  print(paste0("Coverage for the modified case under dense case n^(-0.25), ", 
               names(ci_coverage_rate)[i], ": ", coverage_percentage, "%"))
}



###################################################################################
# Load the saved bootstrap results for the dense case:n^(-0.40)
loaded_data_1213 <- readRDS(file.path(results_dir, "1213_all_bootstrap_samples_sparsity_5.rds"))
# with n=4000, 200:2000(4000)

# Number of repetitions
num_repeats <- length(loaded_data_1213)

# Initialize matrix to store CI coverage results
ci_coverage <- matrix(0, nrow = num_repeats, ncol = length(true_beta))

# Calculate CI coverage for each repeat
for (repeat_idx in 1:num_repeats) {
  bootstraped_beta_linear_multiplier <- loaded_data_1213[[repeat_idx]]
  
  # Compute confidence intervals for each beta
  ci_lower <- apply(bootstraped_beta_linear_multiplier, 2, quantile, probs = 0.025)
  ci_upper <- apply(bootstraped_beta_linear_multiplier, 2, quantile, probs = 0.975)
  
  # Check if the true beta coefficients fall within the CI
  ci_coverage[repeat_idx, ] <- (true_beta >= ci_lower & true_beta <= ci_upper)
}

# Calculate CI coverage rate for each coefficient
ci_coverage_rate <- colMeans(ci_coverage)
names(ci_coverage_rate) <- paste0("Beta", 0:(length(true_beta) - 1))

# Print CI coverage results with descriptive text and formatted percentages
for (i in seq_along(ci_coverage_rate)) {
  coverage_percentage <- sprintf("%.2f", ci_coverage_rate[i] * 100)
  print(paste0("Coverage for the modified case under dense case n^(-0.40), ", 
               names(ci_coverage_rate)[i], ": ", coverage_percentage, "%"))
}



###################################################################################
# Load the saved bootstrap results for the dense case:n^(-0.25)
loaded_data_1213_original <- readRDS(file.path(results_dir, "1213_original_OLS_all_bootstrap_samples_sparsity_6.rds"))
# with n=4000, 200:2000(4000)

# Number of repetitions
num_repeats <- length(loaded_data_1213_original)

# Initialize matrix to store CI coverage results
ci_coverage <- matrix(0, nrow = num_repeats, ncol = length(true_beta))

# Calculate CI coverage for each repeat
for (repeat_idx in 1:num_repeats) {
  bootstraped_beta_linear_multiplier <- loaded_data_1213_original[[repeat_idx]]
  
  # Compute confidence intervals for each beta
  ci_lower <- apply(bootstraped_beta_linear_multiplier, 2, quantile, probs = 0.025)
  ci_upper <- apply(bootstraped_beta_linear_multiplier, 2, quantile, probs = 0.975)
  
  # Check if the true beta coefficients fall within the CI
  ci_coverage[repeat_idx, ] <- (true_beta >= ci_lower & true_beta <= ci_upper)
}

# Calculate CI coverage rate for each coefficient
ci_coverage_rate <- colMeans(ci_coverage)
names(ci_coverage_rate) <- paste0("Beta", 0:(length(true_beta) - 1))

# Print CI coverage results with descriptive text and formatted percentages
for (i in seq_along(ci_coverage_rate)) {
  coverage_percentage <- sprintf("%.2f", ci_coverage_rate[i] * 100)
  print(paste0("Coverage for the original OLS under dense case n^(-0.25), ", 
               names(ci_coverage_rate)[i], ": ", coverage_percentage, "%"))
}





###################################################################################
# Load the saved bootstrap results for the dense case:n^(-0.40)
loaded_data_1213_original <- readRDS(file.path(results_dir, "1213_original_OLS_all_bootstrap_samples_sparsity_5.rds"))
# with n=4000, 200:2000(4000)

# Number of repetitions
num_repeats <- length(loaded_data_1213_original)

# Initialize matrix to store CI coverage results
ci_coverage <- matrix(0, nrow = num_repeats, ncol = length(true_beta))

# Calculate CI coverage for each repeat
for (repeat_idx in 1:num_repeats) {
  bootstraped_beta_linear_multiplier <- loaded_data_1213_original[[repeat_idx]]
  
  # Compute confidence intervals for each beta
  ci_lower <- apply(bootstraped_beta_linear_multiplier, 2, quantile, probs = 0.025)
  ci_upper <- apply(bootstraped_beta_linear_multiplier, 2, quantile, probs = 0.975)
  
  # Check if the true beta coefficients fall within the CI
  ci_coverage[repeat_idx, ] <- (true_beta >= ci_lower & true_beta <= ci_upper)
}

# Calculate CI coverage rate for each coefficient
ci_coverage_rate <- colMeans(ci_coverage)
names(ci_coverage_rate) <- paste0("Beta", 0:(length(true_beta) - 1))

# Print CI coverage results with descriptive text and formatted percentages
for (i in seq_along(ci_coverage_rate)) {
  coverage_percentage <- sprintf("%.2f", ci_coverage_rate[i] * 100)
  print(paste0("Coverage for the original OLS under dense case n^(-0.40), ", 
               names(ci_coverage_rate)[i], ": ", coverage_percentage, "%"))
}
