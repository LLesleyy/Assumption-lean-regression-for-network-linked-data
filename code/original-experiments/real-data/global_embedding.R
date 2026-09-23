
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
library(abind)
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

selected_school_indices<-c()

# processing for each school
num_dimensions_to_keep <- 2  # Example: keep the top 2 dimensions

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
  
  graph_size <- nrow(A_combined)
  di <- colSums(A_combined)
  nu_obv <- sum(A_combined) / (2 * choose(graph_size, 2))
  
  # Step 7: Compute Spectral Embedding
  # Determine the number of dimensions to keep based on the eigen gap or other criteria

  spectral_embedding <- compute_spectral_embedding(A_combined, num_dimensions_to_keep)
  spectral_embedding<- nu_obv^(-0.5)*spectral_embedding
  
  truncated_SVD <- compute_truncated_svd_embedding(A_combined, num_dimensions_to_keep)
  
  # Combine the spectral embedding and SVD with eligible students data
  colnames(spectral_embedding) <- paste0("Spectral_", 1:num_dimensions_to_keep)
  eligible_students <- cbind(eligible_students, spectral_embedding)
  
  colnames(truncated_SVD) <- paste0("SVD_", 1:num_dimensions_to_keep)
  eligible_students <- cbind(eligible_students, truncated_SVD)
  
  
  # Also compute the local subgraph count as potential covariate.
  rooted_ts <- choose(di, 2)
  rooted_ts.density <- rooted_ts / (nu_obv^2 * choose(graph_size - 1, 2))
  eligible_students <- cbind(eligible_students, rooted_ts.density, graph_size, nu_obv)
  
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

# Dimensions
d <- num_dimensions_to_keep
global_dim <- d * length(treated_school_indices)  # max possible dimension, will be reduced below
padded_embeddings_list <- list()

# Go through each treated school
for (k in treated_school_indices) {
  school_data <- filtered_data_combined[[k]]
  
  # Skip if NULL or missing spectral embedding
  if (is.null(school_data) || !any(grepl("^Spectral_", colnames(school_data)))) next
  
  spectral_data <- school_data %>% select(starts_with("Spectral_")) %>% as.matrix()
  n_students <- nrow(spectral_data)
  
  # Initialize padded matrix
  padded_matrix <- matrix(0, nrow = n_students, ncol = global_dim)
  
  # Fill the correct d-dimensional block
  col_range <- ((school_counter - 1) * d + 1):(school_counter * d)
  padded_matrix[, col_range] <- spectral_data
  
  padded_embeddings_list[[length(padded_embeddings_list) + 1]] <- padded_matrix
}

# Update final dimension based on valid schools only
global_dim <- d * (school_counter - 1)
global_spectral_embeddings <- do.call(rbind, padded_embeddings_list)

# Bind with student data
all_eligible_data_combined_with_embedding <- bind_cols(
  all_eligible_data_combined,
  as.data.frame(global_spectral_embeddings)
)

# Rename columns
colnames(all_eligible_data_combined_with_embedding)[
  (ncol(all_eligible_data_combined) + 1):ncol(all_eligible_data_combined_with_embedding)
] <- paste0("GlobalSpectral_", 1:global_dim)


# overwrite to preserve previous datset name.
all_eligible_data_combined<-all_eligible_data_combined_with_embedding


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
         y_raw_wave2, y_wave2, y_wave2_v2, y_raw_wave1, #Spectral_1, Spectral_2, SVD_1, SVD_2, 
         PN1W2, PN2W2, PN3W2, PN4W2, PN5W2, PN6W2, 
         PN7W2, PN8W2, PN9W2, PN10W2, PN11W2, PN12W2, PN13W2,
         PN1, PN2, PN3, PN4, PN5, PN6,PN7, PN8, PN9, PN10, PN11, PN12, OTHETH,
         starts_with("numeric_PN"),
         starts_with("wave1_numeric_PN"),
         starts_with("transformed_numeric_PN"),
         starts_with("transformed_wave1_numeric_PN"),
         starts_with("ETH"),         # Retain race variables starting with ETH
         starts_with("RTS"),         # Recording of attending the meeting.
         starts_with("Spectral_"),   # individual spectral embedding
         starts_with("SVD_"),        # individual SVD
         starts_with("GlobalSpectral_")) # global padded spectral embedding
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



################################################################
######################## Model Fitting #########################
################################################################

# Get all column names starting with "GlobalSpectral_"
global_spectral_vars <- grep("^GlobalSpectral_", names(all_eligible_data_combined), value = TRUE)

# Build formula dynamically
formula_full <- reformulate(
  c( "GRC_NUMERIC", "RETURN", "TREAT_BINARY", "SCHID", 
     "Race_overall", "GENDER", "LIVEWB", global_spectral_vars),
     response = "y_raw_wave2"
)

# Fit model
real_model_full1 <- lm(formula_full, data = all_eligible_data_combined)
summary(real_model_full1)


# refine the model by merging some schools.
# Create a vector of significant SCHID levels
significant_schools <- c("19", "27", "34", "35", "45", "49", "51", "56")
# Create a new variable SCHID_collapsed
all_eligible_data_combined$SCHID_collapsed <- ifelse(
  as.character(all_eligible_data_combined$SCHID) %in% significant_schools,
  as.character(all_eligible_data_combined$SCHID),
  "Other"
)

# Convert it to a factor
all_eligible_data_combined$SCHID_collapsed <- factor(all_eligible_data_combined$SCHID_collapsed)

# Replace SCHID with SCHID_collapsed in the variable list
vars_with_collapse <- c(
  "GRC_NUMERIC", "RETURN", "TREAT_BINARY", "SCHID_collapsed", 
  "Race_overall", "GENDER", "LIVEWB", global_spectral_vars
)

# Rebuild formula with SCHID_collapsed
formula_with_collapse <- reformulate(vars_with_collapse, response = "y_raw_wave2")

# Fit the refined model
Re_real_model_full1 <- lm(formula_with_collapse, data = all_eligible_data_combined)
summary(Re_real_model_full1)

##################################################
# bootstrap for the regression effect estimates. #
##################################################

# Prepare design matrix (X) and response vector (y)
# Use the same formula as in lm()
X_matrix <- model.matrix(formula_with_collapse, data = all_eligible_data_combined)

# Response vector remains the same
y_vec <- all_eligible_data_combined$y_raw_wave2

# Compute OLS estimator
model_beta_obs <- solve(t(X_matrix) %*% X_matrix) %*% t(X_matrix) %*% matrix(y_vec, nrow = total_students)
model_beta_obs_v2 <-Re_real_model_full1$coefficients


# Bootstrap procedure
B_bootstrap<-300;
model_beta_boot <- matrix(0, nrow = length(model_beta_obs), ncol = B_bootstrap)  # Store beta for each bootstrap


for (B in 1:B_bootstrap) {
  epsilon <- rnorm(total_students, 0, 1)  # Gaussian multiplier
  
  # Generate two observed statistics
  observed_product <- matrix(0, nrow = total_students, ncol = length(model_beta_obs))
  for (ii in 1:total_students) {
    observed_product[ii, ] <- X_matrix[ii, ] * y_vec[ii]  
  }
  observed_linear <- colMeans(observed_product)
  observed_quad <- (1 / total_students) * Reduce(`+`, lapply(1:total_students, function(i) outer(X_matrix[i, ], X_matrix[i, ])))
  
  # Perturbation for bootstrap
  xi_boot_perturb_linear <- matrix(0, nrow = total_students, ncol = length(model_beta_obs))
  xi_boot_pertub_quad <- vector("list", total_students)
  
  for (kk in 1:total_students) {
    xi_obs <- X_matrix[kk, ]
    xi2_obs_mat <- outer(xi_obs, xi_obs)
    xi_boot_perturb_linear[kk, ] <- epsilon[kk] * (xi_obs * y_vec[kk] - observed_linear)
    xi_boot_pertub_quad[[kk]] <- epsilon[kk] * (xi2_obs_mat - observed_quad)
  }
  
  # Generate bootstrap sample
  average_x_boot <- observed_linear + colMeans(xi_boot_perturb_linear)
  average_x2_boot <- observed_quad + apply(abind(xi_boot_pertub_quad, along = 3), c(1, 2), mean)
  
  # Compute bootstrapped beta coefficients
  model_beta_boot[, B] <- solve(average_x2_boot) %*% matrix(average_x_boot, ncol = 1)
  
  print(paste("Finished Running bootstrap iteration:", B))
}




B_bootstrap=300;
# detect existence of any overall network effect carried by the top spectral embeddings.
n_covaraite <- length(model_beta_obs)
model_Tn_boot <- numeric(B_bootstrap)
for (B in 1:B_bootstrap) {
  model_Tn_boot[B] = total_students * sum((model_beta_boot[(n_covaraite-55):n_covaraite, B] - model_beta_obs[(n_covaraite-55):n_covaraite])^2)
}
model_Tn_obs = total_students * sum(model_beta_obs[(n_covaraite-55):n_covaraite]^2)

# Compute two-sided p-value
p_upper_network <- mean(model_Tn_boot >= model_Tn_obs)  
p_lower_network <- mean(model_Tn_boot <= model_Tn_obs)  
p_value_network <- 2 * min(p_upper_network, p_lower_network)
p_value_network <- min(p_value_network, 1)
p_value_network





# re-generate the p-value based on the bootstrap variance estimates.
# Compute bootstrap variance and standard errors
beta_var_boot <- apply(model_beta_boot[,1:300], 1, var)  # Variance across bootstrap samples
beta_se_boot <- sqrt(beta_var_boot)  # Standard errors


# Compute updated p-values using bootstrap standard errors
z_scores <- abs(model_beta_obs_v2 / beta_se_boot)  # Compute z-scores
p_values_boot <- 2 * (1 - pnorm(z_scores))  # Two-tailed p-values

# Create a summary table
p_value_table <- data.frame(
  Coefficient = names(model_beta_obs_v2),
  Estimate = model_beta_obs_v2,
  Bootstrap_SE = beta_se_boot,
  P_Value = p_values_boot
)

# Print the updated p-value table
print(p_value_table)

