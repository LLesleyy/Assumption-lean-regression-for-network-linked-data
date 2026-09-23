
# Set parameters for downsampling
m <- 0.5*nrow(all_eligible_data_combined)  # Number of samples to include in each partial summation
k <- 2000   # Number of downsampling iterations

# Prepare design matrix (X) and response vector (y)
X <- model.matrix(y_raw_wave2 ~ y_raw_wave1 + GRC_NUMERIC + RETURN + TREATMENT_FACTOR + 
                    NEIGHBORHOOD_AVG_TREATMENT + SCHID + Race_overall + GENDER + LIVEWB + node_degree.density, 
                  data = all_eligible_data_combined)
y <- all_eligible_data_combined$y_raw_wave2

# Calculate nodewise products for XtX and Xty
nodewise_XtX <- lapply(1:nrow(X), function(i) X[i, ] %*% t(X[i, ]))
nodewise_Xty <- lapply(1:nrow(X), function(i) X[i, ] * y[i])

# Initialize a matrix to store beta estimates from each iteration
beta_estimates <- matrix(0, nrow = k, ncol = ncol(X))

# Perform downsampling k times
set.seed(42)  # For reproducibility
for (i in 1:k) {
  # Randomly sample m out of n nodewise products
  sample_indices <- sample(1:length(nodewise_XtX), m, replace = TRUE)
  
  # Sum the selected nodewise products for this iteration
  XtX_sample_sum <- Reduce(`+`, nodewise_XtX[sample_indices])
  Xty_sample_sum <- Reduce(`+`, nodewise_Xty[sample_indices])
  
  # Calculate downsampled beta for this iteration
  beta_hat_downsampled <- solve(XtX_sample_sum) %*% Xty_sample_sum
  
  # Store the result
  beta_estimates[i, ] <- beta_hat_downsampled
}

# `beta_estimates` now holds all k downsampled beta estimates, one per row
beta_estimates  # Each row corresponds to the coefficient estimates for an iteration

variance_matrix<-var(beta_estimates)


# Compute the mean of downsampled beta estimates
beta_mean <- colMeans(beta_estimates)

# Calculate the standard errors from the downsampled variance matrix
standard_errors <- sqrt(diag(variance_matrix))

# Calculate t-statistics for each coefficient
t_statistics <- beta_mean / standard_errors

# Calculate p-values for each coefficient
# Assuming large sample sizes, using normal distribution for approximation
p_values <- 2 * (1 - pnorm(abs(t_statistics)))

# Variable names corresponding to each coefficient
variable_names <- colnames(X)

# Create a summary table for coefficients, standard errors, t-statistics, and p-values
summary_table <- data.frame(
  Coefficient = beta_mean,
  Std_Error = standard_errors,
  t_value = t_statistics,
  p_value = p_values,
  row.names = variable_names
)

# Display the summary table
summary_table





