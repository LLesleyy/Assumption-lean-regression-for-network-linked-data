
# filter out a few schools that 
# the large bootstrap variance is the key factor of preventing find significance.
# down-sampling helps reduce variance.
# sample splitting for model selection
# evidence of non-linearity as to compare the ratio the between the default std and our bootstrap variance
# count functional might be be appropriate for inference using the full samples, but downsampling should work for this case.
# this might not be a good example for using regressing with the count functional and test the bias-corrected estimator in the standard full sample regime,
# since the aggregation might lead to the unreasonable sparse graph.



# Required libraries
library(dplyr)  # For data manipulation
library(tidyr)  # For handling missing values
library(ggplot2) # For plotting
library(igraph)  # For spectral decomposition
library(AID)

# Load the packages
# library(sandwich)
# library(lmtest)


args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[1L] else file.path("data", "ProcessedData.Rda")
output_folder <- if (length(args) >= 2L) args[2L] else file.path("results", "real-data")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
load(data_file)

# Function to compute spectral embedding with sorted eigenvalues
compute_spectral_embedding <- function(A, num_dimensions) {
  # Compute the eigen decomposition of the adjacency matrix
  eig <- eigen(A, symmetric = TRUE)
  
  # Sort eigenvalues in descending order and rearrange eigen-vectors accordingly
  sorted_indices <- order(abs(eig$values), decreasing = TRUE)  # Sort by magnitude in descending order
  sorted_eigenvalues <- eig$values[sorted_indices]
  sorted_eigenvectors <- eig$vectors[, sorted_indices]
  
  # Select the top eigen-vectors and eigenvalues based on the sorted order
  selected_vectors <- sorted_eigenvectors[, 1:num_dimensions]
  selected_eigenvalues <- abs(sorted_eigenvalues[1:num_dimensions])^(1/2)
  
  # Apply the scaling to the selected eigenvectors
  spectral_embedding <- selected_vectors %*% diag(selected_eigenvalues)
  
  return(spectral_embedding)
}

# Alternative function for Truncated SVD embedding with sorted singular values
compute_truncated_svd_embedding <- function(A, num_dimensions) {
  # Perform Truncated SVD using svd()
  svd_result <- svd(A)
  
  # Sort singular values in descending order and rearrange singular vectors accordingly
  sorted_indices <- order(svd_result$d, decreasing = TRUE)  # Sort singular values in descending order
  sorted_singular_values <- svd_result$d[sorted_indices]
  sorted_u_vectors <- svd_result$u[, sorted_indices]
  
  # Select the top singular vectors and singular values based on the sorted order
  selected_vectors <- sorted_u_vectors[, 1:num_dimensions]
  selected_singular_values <- sorted_singular_values[1:num_dimensions]
  
  # Scale the singular vectors by the singular values
  svd_embedding <- selected_vectors %*% diag(selected_singular_values)
  
  return(svd_embedding)
}


# Correct function to map survey responses to numeric values, fully vectorized
map_survey_responses <- function(x) {
  levels <- c("(0) Almost nobody", "(1) A few people", "(2) About 25%", 
              "(3) About 50%", "(4) About 75%", "(5) Almost everyone")
  numeric_values <- 0:5
  mapped_values <- ifelse(is.na(x), NA, numeric_values[match(x, levels)])  # Use ifelse to handle NAs properly
  return(mapped_values)
}

map_obs_responses <- function(x) {
  levels <- c("(0) Never", "(1) 1-2 times/month", "(2) About 1 time/week", 
              "(3) 2-3 times/week", "(4) Every day")
  numeric_values <- 0:4
  mapped_values <- ifelse(is.na(x), NA, numeric_values[match(x, levels)])  # Use ifelse to handle NAs properly
  return(mapped_values)
}

# Function to apply different powers and visualize the results
check_normality_of_transformed_y <- function(data, powers = c(1, 1.2, 0.5, 2, -1)) {
  for (power in powers) {
    # Apply the power transformation
    transformed_y <- data$y_raw_wave2^power
    
    # Create a histogram and Q-Q plot for the transformed data
    p1 <- ggplot(data.frame(transformed_y), aes(x = transformed_y)) +
      geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black") +
      ggtitle(paste("Histogram of y^", power)) +
      xlab(paste("y^", power)) + theme_minimal()
    
    p2 <- ggplot(data.frame(transformed_y), aes(sample = transformed_y)) +
      stat_qq() + stat_qq_line() +
      ggtitle(paste("Q-Q Plot of y^", power)) +
      theme_minimal()

  
    # View the test results
    cat("Box-Cox Test for y^", power, ":\n")
    print(boxcoxnc(y = transformed_y))
    print(p1)  # Display histogram
    print(p2)  # Display Q-Q plot
  }
}




# Initialize a vector to store which schools are treated
treated_schools <- vector("logical", length = length(lcc.list))
# Identify indices of treated schools
for (i in 1:length(lcc.list)) {
  df <- lcc.list[[i]][["df"]]  # Extract the dataframe for demographic features
  treated_schools[i] <- any(grepl("Treatment", df[["SCHTREAT"]]))  # Check if the school is treated
}

treated_school_indices <- which(treated_schools)  # Indices of schools that are treated
cat("Number of treated schools:", length(treated_school_indices), "\n")





# Initialize list to store filtered data for treated schools
filtered_data_combined <- list()
# Loop over only the treated schools
for (i in treated_school_indices) {
  
  # Extract the school data for the current treated school
  school_data <- lcc.list[[i]]
  df <- school_data[["df"]]  # Extract the dataframe for demographic features
  
  # Total number of students in the school before any filtering
  total_students <- nrow(df)
  
  ########## Step 1: Filter out students with missing demographic information
  ########## the FIRST place we filter out records
  eligible_students <- df %>%
    filter(!is.na(TREAT) & !is.na(SCHID) & !is.na(GENDER) & 
             !is.na(GRC) & !is.na(RETURN) & !is.na(LIVEWB) & 
             !(is.na(ETHW) & is.na(ETHB) & is.na(ETHH) & is.na(ETHA) & is.na(ETHC) & is.na(ETHSA) & is.na(ETHO)))
  
  ####################### Step 2: Convert survey responses from wave 2 to numeric values
  survey_cols_wave2 <- c("PN1W2", "PN2W2", "PN3W2", "PN4W2", "PN5W2", "PN6W2", 
                         "PN7W2", "PN8W2", "PN9W2", "PN10W2", "PN11W2", "PN12W2", "PN13W2")

  # some of the questions were about negative atmosphere. 
  # So we transform then by the 5− the original scores before combining with other questions to ensure all scores are measuring positive atmosphere.
  # PN3W2 PN4W2 PN5W2 PN6W2.  
  # the reference study only picks two of them.
  eligible_students <- eligible_students %>%
    mutate(across(all_of(survey_cols_wave2), ~ map_survey_responses(as.character(.)), .names = "numeric_{col}"))
  # check if scores are appropriately transformed.
  
  # Step 3: Filter out students with all NAs in wave 2 survey responses
  # The SECOND place we filter out records
  eligible_students <- eligible_students %>%
    rowwise() %>%
    filter(sum(is.na(c_across(starts_with("numeric_PN")))) < length(survey_cols_wave2))
  
  # Step 3*: Filter out students with any NAs in the wave 2 survey responses
  #eligible_students <- eligible_students %>%
  #  filter(if_all(starts_with("numeric_PN"), ~ !is.na(.)))
  
  # List of columns that represent the negative atmosphere (to be transformed)
   negative_atmosphere_cols <- c("numeric_PN3W2", "numeric_PN4W2", "numeric_PN5W2", "numeric_PN6W2")
  #negative_atmosphere_cols <- c("numeric_PN4W2","numeric_PN6W2")
  
  # Transform the scores for negative atmosphere questions (5 - original score)
  eligible_students <- eligible_students %>%
    mutate(across(all_of(negative_atmosphere_cols), ~ 5 - ., .names = "transformed_{col}"))
  
  # List of columns representing all survey questions, including the transformed negative atmosphere questions
  all_survey_cols <- c("numeric_PN1W2", "numeric_PN2W2", "transformed_numeric_PN3W2", "transformed_numeric_PN4W2",
                       "transformed_numeric_PN5W2", "transformed_numeric_PN6W2", "numeric_PN7W2", "numeric_PN8W2", 
                       "numeric_PN9W2", "numeric_PN10W2", "numeric_PN11W2", "numeric_PN12W2", "numeric_PN13W2")
  
  #all_survey_cols <- c("numeric_PN1W2", "numeric_PN2W2", "transformed_numeric_PN4W2","transformed_numeric_PN6W2", 
  #                     "numeric_PN3W2", "numeric_PN5W2", "numeric_PN7W2", "numeric_PN8W2", 
  #                     "numeric_PN9W2", "numeric_PN10W2", "numeric_PN11W2", "numeric_PN12W2", "numeric_PN13W2")
  
  # Calculate the average score across all available questions
  eligible_students <- eligible_students %>%
    rowwise() %>%
    mutate(y_raw_wave2 = mean(c_across(all_of(all_survey_cols)), na.rm = TRUE))  # Calculate average, excluding NA values
  
  ##################### repeat the same procedure for the wave-1 response.
  survey_cols_wave1 <- c("PN1", "PN2", "PN3", "PN4", "PN5", "PN6","PN7", "PN8", "PN9", "PN10", "PN11", "PN12")
  
  eligible_students <- eligible_students %>%
    mutate(across(all_of(survey_cols_wave1), ~ map_survey_responses(as.character(.)), .names = "wave1_numeric_{col}"))
  # check if scores are appropriately transformed.
  
  # Step 4: Filter out students with all NAs in wave 1 survey responses
  eligible_students <- eligible_students %>%
    rowwise() %>%
    filter(sum(is.na(c_across(starts_with("wave1_numeric_PN")))) < length(survey_cols_wave1))
  
  # List of columns that represent the negative atmosphere (to be transformed)
  negative_atmosphere_cols_wave1 <- c("wave1_numeric_PN3", "wave1_numeric_PN4", 
                                      "wave1_numeric_PN5", "wave1_numeric_PN6","wave1_numeric_PN9")

  # Transform the scores for negative atmosphere questions (5 - original score)
  eligible_students <- eligible_students %>%
    mutate(across(all_of(negative_atmosphere_cols_wave1), ~ 5 - ., .names = "transformed_{col}"))
  
  # List of columns representing all survey questions, including the transformed negative atmosphere questions
  all_survey_cols_wave1 <- c("wave1_numeric_PN1", "wave1_numeric_PN2", 
                             "transformed_wave1_numeric_PN3", "transformed_wave1_numeric_PN4",
                             "transformed_wave1_numeric_PN5", "transformed_wave1_numeric_PN6", 
                             "wave1_numeric_PN7", "wave1_numeric_PN8", 
                             "transformed_wave1_numeric_PN9", 
                             "wave1_numeric_PN10", "wave1_numeric_PN11", "wave1_numeric_PN12")
  
  # Calculate the average score across all available questions
  eligible_students <- eligible_students %>%
    rowwise() %>%
    mutate(y_raw_wave1 = mean(c_across(all_of(all_survey_cols_wave1)), na.rm = TRUE))  # Calculate average, excluding NA values
  
  
  
  ################### general info: Number of students remaining after filtering
  num_remaining <- nrow(eligible_students)
  
  # Count and display the number of students remove
  num_removed <- total_students - num_remaining
  cat("School", i, ": Total students =", total_students, ", Students removed =", num_removed, ", Students remaining =", num_remaining, "\n")
  
  # Apply the transformation: raise to the power of 1.2 for wave 2 data
  if (nrow(eligible_students) > 0) {
    eligible_students <- eligible_students %>%
      mutate(y_wave2 = y_raw_wave2^1.2,
             y_wave2_v2 = (y_raw_wave2-1)^0.9)
  } else {
    warning(paste("No valid non-NA values for transformation in school", i))
    next
  }
  
  ################## Step 5: Create a new binary variable TREAT_BINARY
  eligible_students <- eligible_students %>%
    mutate(TREAT_BINARY = ifelse(grepl("Treatment", TREAT), 1, 0))
  # and 0 by default will be baseline level.

  ################## Step 5': the treatment assignment and the participation of the 
  ################## workshop might differ, also process an alternative variable
  ################## whether this student attend any one of this intervention meeting.
  meeting_attendence_cols <- c("RTSM1", "RTSM2", "RTSM3", "RTSM4","RTSM5", "RTSM6", "RTSM7", "RTSM8", "RTSM9")
  eligible_students <- eligible_students %>%
    rowwise() %>%
    mutate(
      TREATMENT_FACTOR = ifelse(any(c_across(all_of(meeting_attendence_cols)) == "(1) Yes", na.rm = TRUE), 1, 0)
    ) %>%
    ungroup()  # Remove row-wise grouping after the operation
  
  ################# Correct extraction of GRC_NUMERIC from GRC
  eligible_students <- eligible_students %>%
    mutate(GRC_NUMERIC = case_when(
      GRC == "(0) 5th grade" ~ 5,
      GRC == "(1) 6th grade" ~ 6,
      GRC == "(2) 7th grade" ~ 7,
      GRC == "(3) 8th grade" ~ 8,
      TRUE ~ NA_real_  # If any unexpected values are found
    ))
  # Check the summary of the corrected GRC_NUMERIC
  summary(eligible_students$GRC_NUMERIC)
  
  ################ Step 6: Construct the combined adjacency matrix A_combined using A1 and A2
  A1 <- school_data[["A1"]]
  A2 <- school_data[["A2"]]
  
  # Use row indices to subset the adjacency matrices
  eligible_indices <- as.integer(rownames(eligible_students))
  A1 <- A1[eligible_indices, eligible_indices, drop = FALSE]
  A2 <- A2[eligible_indices, eligible_indices, drop = FALSE]
  
  # Construct the combined adjacency matrix
  A_combined <- pmax(A1, A2)  # Element-wise maximum to combine A1 and A2
  A_combined_min <- A1*A2 # Element-wise product to combine A1 and A2 as the minimum network
  
  # Step 7: Compute Spectral Embedding
  # Determine the number of dimensions to keep based on the eigen gap or other criteria
  num_dimensions_to_keep <- 2  # Example: keep the top 2 dimensions
  spectral_embedding <- compute_spectral_embedding(A_combined, num_dimensions_to_keep)
  truncated_SVD <- compute_truncated_svd_embedding(A_combined, num_dimensions_to_keep)
  
  # Combine the spectral embedding and SVD with eligible students data
  colnames(spectral_embedding) <- paste0("Spectral_", 1:num_dimensions_to_keep)
  eligible_students <- cbind(eligible_students, spectral_embedding)
  
  colnames(truncated_SVD) <- paste0("SVD_", 1:num_dimensions_to_keep)
  eligible_students <- cbind(eligible_students, truncated_SVD)
  
  
  # Also compute the local subgraph count as potential covariate.
  graph_size <- nrow(A_combined)
  di <- colSums(A_combined)
  nu_obv <- sum(A_combined) / (2 * choose(graph_size, 2))
  
  rooted_ts <- choose(di, 2)
  rooted_ts.density <- rooted_ts / (nu_obv^2 * choose(graph_size - 1, 2))
  eligible_students <- cbind(eligible_students, rooted_ts.density)
  
  node_degree.density<-di/(nu_obv*(graph_size - 1))
  eligible_students <- cbind(eligible_students, node_degree.density)
  
  # also construct the neighbourhood average of the treatment effect.
  # Extract TREATMENT_FACTOR as a vector
  treatment_vector <- eligible_students$TREATMENT_FACTOR
  # Calculate the sum of TREATMENT_FACTOR for each student's neighbors
  neighborhood_sum <- A_combined %*% treatment_vector
  # Calculate the number of neighbors for each student
  num_neighbors <- rowSums(A_combined)
  # Compute the neighborhood average, handling cases with no neighbors to avoid division by zero
  neighborhood_avg_treatment <- ifelse(num_neighbors > 0, neighborhood_sum / num_neighbors, 0)
  # Add the result to the eligible_students dataframe
  eligible_students <- eligible_students %>%
    mutate(NEIGHBORHOOD_AVG_TREATMENT = neighborhood_avg_treatment)
  
  

  # Append processed data to the filtered list for combined analysis
  filtered_data_combined[[i]] <- eligible_students
}

# Combine all eligible students across treated schools for combined analysis
all_eligible_data_combined <- bind_rows(filtered_data_combined)

# Merge ETHC and ETHSA and all other levels into OTHER and create one categorical variable.
all_eligible_data_combined <- all_eligible_data_combined %>%
  mutate(Race_overall = ifelse(!(ETHW == "(1) White ethnicity" | ETHB == "(1) Black ethnicity" | ETHH == "(1) Hispanic" | ETHA == "(1) Asian-American"), 
                             "Others", 
                             ifelse(ETHW == "(1) White ethnicity", "White", 
                                    ifelse(ETHB == "(1) Black ethnicity", "Black", 
                                           ifelse(ETHH == "(1) Hispanic", "Hispanic", "Asian")))))
# Reorder the Race_overall factor to make "Others" the baseline level
all_eligible_data_combined <- all_eligible_data_combined %>%
  mutate(Race_overall = relevel(factor(Race_overall), ref = "Others"))

# only select the relevant variables into the final data-set.
# GR is wave-1 grade info
# transformed_wave1_numeric_PN4W2
all_eligible_data_combined <- all_eligible_data_combined %>%
  select(SCHID, GRC, GR, GRC_NUMERIC, RETURN, TREAT, TREAT_BINARY, TREATMENT_FACTOR, NEIGHBORHOOD_AVG_TREATMENT,
         Race_overall, GENDER, LIVEWB, rooted_ts.density, node_degree.density,
         y_raw_wave2, y_wave2, y_wave2_v2, y_raw_wave1, Spectral_1, Spectral_2, SVD_1, SVD_2, 
         PN1W2, PN2W2, PN3W2, PN4W2, PN5W2, PN6W2, 
         PN7W2, PN8W2, PN9W2, PN10W2, PN11W2, PN12W2, PN13W2,
         PN1, PN2, PN3, PN4, PN5, PN6,PN7, PN8, PN9, PN10, PN11, PN12, OTHETH,
         starts_with("numeric_PN"),
         starts_with("wave1_numeric_PN"),
         starts_with("transformed_numeric_PN"),
         starts_with("transformed_wave1_numeric_PN"),
         starts_with("ETH"), # Retain race variables starting with ETH
         starts_with("RTS")) # Recording of attending the meeting.
# check all the variable names.
str(all_eligible_data_combined)

# Step 8: Check number of eligible treated schools
num_eligible_schools <- length(treated_school_indices)
cat("Number of eligible treated schools:", num_eligible_schools, "\n")

if (num_eligible_schools != 28) {
  stop("The number of eligible treated schools is not equal to 28. Please check the filtering criteria.")
}

# Print the total number of students included in the regression model for combined analysis
total_students <- nrow(all_eligible_data_combined)
cat("Total number of students included in the regression model for combined analysis:", total_students, "\n")


real_model_controled_us_PN <- lm(y_raw_wave2 ~ y_raw_wave1 + GRC_NUMERIC + RETURN + TREATMENT_FACTOR + NEIGHBORHOOD_AVG_TREATMENT + SCHID + Race_overall + GENDER + LIVEWB + node_degree.density, 
                           data = all_eligible_data_combined)
summary(real_model_controled_us_PN)
original_coefficients <- coef(real_model_controled_us_PN)

   
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

# downsample variance estimate.
variance_matrix<-var(beta_estimates)

# Calculate the standard errors from the downsampled variance matrix
standard_errors <- sqrt(diag(variance_matrix))

# Calculate t-statistics for each coefficient
t_statistics <- original_coefficients / standard_errors

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






































real_model_controled <- lm(y_wave2_v2 ~ GRC_NUMERIC + RETURN + TREATMENT_FACTOR + SCHID + Race_overall + GENDER + LIVEWB, 
                         data = all_eligible_data_combined)
summary(real_model_controled)


###################### following are the models for comparison.
# Perform Regression for wave 2 using spectral embedding and interaction term
# Benchmark OLS fitting without any network covariate over transformed friendliness score.
real_model_OLS_SCH <- lm(y_wave2 ~ GRC_NUMERIC + RETURN + TREAT_BINARY + SCHID + Race_overall + GENDER + LIVEWB, 
                   data = all_eligible_data_combined)
summary(real_model_OLS_SCH)

real_model_OLS <- lm(y_wave2 ~ GRC_NUMERIC + RETURN + TREAT_BINARY  + Race_overall + GENDER + LIVEWB, 
                   data = all_eligible_data_combined)
summary(real_model_OLS)

real_model_OLS_SCH_v2 <- lm(y_wave2_v2~ GRC_NUMERIC + RETURN + TREAT_BINARY  + SCHID + Race_overall + GENDER + LIVEWB, 
                     data = all_eligible_data_combined)
summary(real_model_OLS_SCH_v2)

# comparable fitting with spectral embedding.
real_model_1 <- lm(y_wave2 ~ GRC_NUMERIC + RETURN + TREAT_BINARY + SCHID+ Race_overall + GENDER + LIVEWB + Spectral_1 + Spectral_2, 
                   data = all_eligible_data_combined)
summary(real_model_1)
real_model_1_SVD <- lm(y_wave2 ~ GRC_NUMERIC + RETURN + TREAT_BINARY + SCHID+ Race_overall + GENDER + LIVEWB + SVD_1 + SVD_2, 
                   data = all_eligible_data_combined)
summary(real_model_1_SVD)

# our target model: without distribution assumption on y/intial rating control/rdpg embedding
real_model_1_raw <- lm(y_raw_wave2 ~ y_raw_wave1 + GRC_NUMERIC + RETURN + TREAT_BINARY + SCHID+ Race_overall + GENDER + LIVEWB + Spectral_1 + Spectral_2, 
                   data = all_eligible_data_combined)
summary(real_model_1_raw)



# interaction between grade and returning students
real_model_2 <- lm(y_raw_wave2 ~ GRC_NUMERIC * RETURN + TREAT_BINARY  + SCHID+ Race_overall + GENDER + LIVEWB + Spectral_1 + Spectral_2, 
                                      data = all_eligible_data_combined)
summary(real_model_2)








########################## chi-square test for the overall network effect over the spectral vector:
#  Compute the Sandwich Covariance Matrix
sandwich_cov <- vcovHC(real_model_1_raw, type = "HC")  # Heteroskedasticity-robust covariance matrix

# Extract the Submatrix for the Spectral Regressors (Spectral_1 and Spectral_2)
spectral_indices <- which(names(coef(real_model_4)) %in% c("Spectral_1", "Spectral_2"))
spectral_cov_submatrix <- sandwich_cov[spectral_indices, spectral_indices]

# Step 3: Compute the Quadratic Form
spectral_coefs <- coef(real_model_4)[spectral_indices]
quadratic_form <- t(spectral_coefs) %*% solve(spectral_cov_submatrix) %*% spectral_coefs

# Step 4: Perform the Chi-Square Test
# Degrees of freedom is equal to the number of spectral dimensions
df <- length(spectral_coefs)

# Compute the p-value
p_value <- pchisq(quadratic_form, df = df, lower.tail = FALSE)

# Output the test statistic and p-value
cat("Chi-square test statistic:", quadratic_form, "\n")
cat("Degrees of freedom:", df, "\n")
cat("P-value for the chi-square test:", p_value, "\n")

# Interpretation
if (p_value < 0.05) {
  cat("The overall network effect (captured by the spectral embedding) is significant.\n")
} else {
  cat("The overall network effect (captured by the spectral embedding) is not significant.\n")
}


# Extract the new standard errors from the sandwich covariance matrix
new_se <- sqrt(diag(sandwich_cov))

# Extract the coefficient estimates from the original model
coefficients_estimates <- coef(real_model_1_raw)

# Recalculate the t-values
new_t_values <- coefficients_estimates / new_se

# Recalculate the p-values based on the new t-values and degrees of freedom
df_residual <- df.residual(real_model_1_raw)  # Get the residual degrees of freedom
new_p_values <- 2 * pt(abs(new_t_values), df = df_residual, lower.tail = FALSE)

# Create a significance column based on the p-values
significance <- cut(new_p_values, 
                    breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), 
                    labels = c("***", "**", "*", ".", ""), 
                    right = FALSE)

# Create a new summary table with updated values and significance
results_table <- data.frame(
  Coefficient = coefficients_estimates,
  Std_Error = new_se,
  t_value = new_t_values,
  p_value = new_p_values,
  Significance = significance
)

# Print the updated regression table
print(results_table)









######################## start model selections

# Split the data-set.
# randomly evenly split the whole data-set into two halves.
rand_inx <- sample(1:total_students, total_students / 2, replace = TRUE)
sp1 <- all_eligible_data_combined[rand_inx, ] # SP1 for model selection

sp2_inx <- setdiff(1:total_students, rand_inx)
sp2 <- all_eligible_data_combined[sp2_inx, ]  # SP2 for final fitting and estimation


# Fit the full model
full_mode_sp1 <- lm(y_raw_wave2 ~ GRC_NUMERIC + RETURN + TREAT_BINARY + SCHID+ Race_overall + GENDER + LIVEWB + Spectral_1 + Spectral_2, 
                 data = sp1)
summary(full_mode_sp1)
# Calculate sandwich covariance matrix
sandwich_cov_matrix <- vcovHC(full_mode_sp1, type = "HC0")

# Obtain the sandwich-based p-values
coeftest_result <- coeftest(full_mode_sp1, sandwich_cov_matrix)
p_values <- coeftest_result[, "Pr(>|t|)"]
print(p_values)



# Identify which levels are associated with School (SCHID) and Race (Race_overall)
school_levels <- grep("^SCHID", names(p_values))
race_levels <- grep("^Race_overall", names(p_values))

# Apply Bonferroni correction
adjusted_school_pvals <- p.adjust(p_values[school_levels], method = "bonferroni")
adjusted_race_pvals <- p.adjust(p_values[race_levels], method = "bonferroni")

# Identify insignificant levels (adjusted p-value >= 0.05)
insignificant_schools <- names(p_values[school_levels])[adjusted_school_pvals >= 0.05]
insignificant_races <- names(p_values[race_levels])[adjusted_race_pvals >= 0.05]

print(insignificant_schools)
print(insignificant_races)


# Extract significant levels after Bonferroni correction
school_levels_to_keep <- names(p_values[school_levels])[adjusted_school_pvals < 0.05]
# Remove the "SCHID" prefix from the names of the levels to keep
school_levels_to_keep <- gsub("SCHID", "", school_levels_to_keep)
# Convert the significant levels to numeric (assuming SCHID is numeric in your dataset)
school_levels_to_keep <- as.numeric(school_levels_to_keep)
print(school_levels_to_keep)


sp1 <- sp1 %>%
  mutate(SCHID_reduced = ifelse(SCHID %in% school_levels_to_keep, SCHID, "sch_baseline"))

sp1 <- sp1 %>%
  mutate(SCHID_reduced = factor(SCHID_reduced, levels = c("sch_baseline", setdiff(unique(SCHID_reduced), "sch_baseline"))))

# Merge insignificant Race levels into reference level (Others)
#race_levels_to_keep <- names(p_values[race_levels])[adjusted_race_pvals < 0.05]
#sp1 <- sp1 %>%
#  mutate(Race_reduced = ifelse(Race_overall %in% race_levels_to_keep, Race_overall, "Others"))


# Refit the model with reduced categorical variables
current_model <- lm(y_raw_wave2 ~ GRC_NUMERIC + RETURN + TREAT_BINARY + SCHID_reduced + GENDER + LIVEWB, data = sp1)

  # Calculate sandwich covariance matrix
  sandwich_cov_matrix <- vcovHC(current_model, type = "HC0")
  # Obtain sandwich-based p-values
  coeftest_result <- coeftest(current_model, sandwich_cov_matrix)
  p_values <- coeftest_result[, "Pr(>|t|)"]
  # keep intercept and treatment (keep them in the model)
  p_values <- p_values[!names(p_values) %in% c("(Intercept)", "TREAT_BINARY")]
  print(p_values)
  # Check if the largest p-value is significant
  max_p <- max(p_values)
  result <- ifelse(max_p >= 0.05, "Greater than or equal to 0.05", "Less than 0.05")
  print(result)
  # Find the variable with the largest p-value
  variable_to_remove <- names(p_values)[which.max(p_values)]
  
current_model <- lm(y_raw_wave2 ~ GRC_NUMERIC + RETURN + TREAT_BINARY + SCHID_reduced  + LIVEWB, data = sp1)

# Calculate sandwich covariance matrix
sandwich_cov_matrix <- vcovHC(current_model, type = "HC0")
# Obtain sandwich-based p-values
coeftest_result <- coeftest(current_model, sandwich_cov_matrix)
p_values <- coeftest_result[, "Pr(>|t|)"]
# keep intercept and treatment (keep them in the model)
p_values <- p_values[!names(p_values) %in% c("(Intercept)", "TREAT_BINARY")]
print(p_values)
# Check if the largest p-value is significant
max_p <- max(p_values)
result <- ifelse(max_p >= 0.05, "Greater than or equal to 0.05", "Less than 0.05")
print(result)

# Print final model coefficients with sandwich estimator
sandwich_cov_matrix <- vcovHC(current_model, type = "HC0")
coeftest_result <- coeftest(current_model, sandwich_cov_matrix)

# Add significance markers
significance_levels <- ifelse(coeftest_result[, "Pr(>|t|)"] < 0.001, "***", 
                              ifelse(coeftest_result[, "Pr(>|t|)"] < 0.01, "**", 
                                     ifelse(coeftest_result[, "Pr(>|t|)"] < 0.05, "*", "")))

# Combine results with significance levels
final_result <- cbind(coeftest_result, Significance = significance_levels)

# Output the final result
print(final_result)



# Using the refined model terms from the final step of the backward elimination
sp2 <- sp2 %>%
  mutate(SCHID_reduced = ifelse(SCHID %in% school_levels_to_keep, SCHID, "sch_baseline"))
sp2 <- sp2 %>%
  mutate(SCHID_reduced = factor(SCHID_reduced, levels = c("sch_baseline", setdiff(unique(SCHID_reduced), "sch_baseline"))))

final_formula <- y_wave2 ~ GRC_NUMERIC + RETURN + TREAT_BINARY + SCHID_reduced + LIVEWB

# Fit the final model on sp2
final_model_sp2 <- lm(final_formula, data = sp2)

# Calculate sandwich covariance matrix for sp2
sandwich_cov_matrix_sp2 <- vcovHC(final_model_sp2, type = "HC0")

# Get the coefficient table with sandwich estimates
final_result_sp2 <- coeftest(final_model_sp2, sandwich_cov_matrix_sp2)

# Print the result
print(final_result_sp2)
