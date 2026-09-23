



# Required libraries
library(dplyr)  # For data manipulation
library(tidyr)  # For handling missing values
library(ggplot2) # For plotting
library(igraph)  # For spectral decomposition
library(AID)
library(Matrix)  # Load Matrix package for block diagonal construction

# Load the packages
#library(sandwich)
#library(lmtest)

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


# for an aggregated large graph.
A_combined_list <- list()

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
graph_summary <- data.frame(School = integer(), Graph_Size = integer(), Sparsity = numeric())

# Loop over only the treated schools
for (i in treated_school_indices) {
  
  # Extract the school data for the current treated school
  school_data <- lcc.list[[i]]
  df <- school_data[["df"]]  # Extract the data-frame for demographic features
  
  # Total number of students in the school before any filtering
  total_students <- nrow(df)
  
  ########## Step 1: Filter out students with missing demographic information
  ########## the FIRST place we filter out records
  eligible_students <- df %>%
    filter(!is.na(TREAT) & !is.na(SCHID) & !is.na(GENDER) & 
             !is.na(GRC) & !is.na(RETURN) & !is.na(LIVEWB) & 
             !(is.na(ETHW) & is.na(ETHB) & is.na(ETHH) & is.na(ETHA) & is.na(ETHC) & is.na(ETHSA) & is.na(ETHO)))
  
  
  #################### there is an other set of scores which are more relevant
  #################### -wave 1 DN scores
  obs_cols_wave1 <- c("DN1", "DN2", "DN3", "DN4", "DN5", "DN6", 
                      "DN7", "DN8", "DN9", "DN10", "DN11", "DN12", "DN13")
  
  # the rating is for the level of conflicts people observed at school
  # So we transform them by the 5− the original scores ；before combining with other questions to ensure all scores are measuring negative level of conflicts here.
  # DN1-DN8: negative behaviors  
  # DN9-DN13: positive actions against conflicts.  
  eligible_students <- eligible_students %>%
    mutate(across(all_of(obs_cols_wave1), ~ map_obs_responses(as.character(.)), .names = "numeric_{col}"))
  # check if scores are appropriately transformed.
  
  # Step 3: Filter out students with all NAs in wave 1 survey responses
  # The SECOND place we filter out records
  eligible_students <- eligible_students %>%
    rowwise() %>%
    filter(sum(is.na(c_across(starts_with("numeric_DN")))) < length(obs_cols_wave1))
  
  # ATTENTION
  # List of columns that represent the positive atmosphere (to be transformed)
  # cautious: here we are capturing the conflicts level
  positive_atmosphere_cols_wave1 <- c("numeric_DN9", "numeric_DN10", "numeric_DN11", "numeric_DN12", "numeric_DN13")
  
  # Transform the scores for negative atmosphere questions (5 - original score)
  eligible_students <- eligible_students %>%
    mutate(across(all_of(positive_atmosphere_cols_wave1), ~ 5 - ., .names = "transformed_{col}"))
  
  # List of columns representing all survey questions, including the transformed negative atmosphere questions
  all_obs_cols_wave1 <- c("numeric_DN1", "numeric_DN2", "numeric_DN3", "numeric_DN4",
                          "numeric_DN5", "numeric_DN6", "numeric_DN7", "numeric_DN8", 
                          "transformed_numeric_DN9", "transformed_numeric_DN10", "transformed_numeric_DN11", "transformed_numeric_DN12", "transformed_numeric_DN13")
  
  # Calculate the average score across all available questions
  eligible_students <- eligible_students %>%
    rowwise() %>%
    mutate(DN_raw_wave1 = mean(c_across(all_of(all_obs_cols_wave1)), na.rm = TRUE))  # Calculate average, excluding NA values
  
  
  
  
  #################### -wave 2 DN scores
  obs_cols_wave2 <- c("DN1W2", "DN2W2", "DN3W2", "DN4W2", "DN5W2", "DN6W2", "DN7W2", "DN8W2", 
                      "DN9W2", "DN10W2", "DN11W2", "DN12W2", "DN13W2", "DN14W2")
  
  # the rating is for the level of conflicts people observed at school
  # So we transform then by the 5− the original scores before combining with other questions to ensure all scores are measuring positive atmosphere.
  # DN1W2-DN8W2: negative.  
  # DN9W2-DN13W2: positive against conflicts. 
  # extra DN14W2 is also negative
  eligible_students <- eligible_students %>%
    mutate(across(all_of(obs_cols_wave2), ~ map_obs_responses(as.character(.)), .names = "wave2_numeric_{col}"))
  # check if scores are appropriately transformed.
  
  # Step 3: Filter out students with all NAs in wave 1 survey responses
  # The SECOND place we filter out records
  eligible_students <- eligible_students %>%
    rowwise() %>%
    filter(sum(is.na(c_across(starts_with("wave2_numeric_DN")))) < length(obs_cols_wave2))
  
  # List of columns that represent the positive atmosphere (to be transformed)
  # cautious: here we are capturing the conflicts level
  positive_atmosphere_cols_wave2 <- c("wave2_numeric_DN9W2", "wave2_numeric_DN10W2", "wave2_numeric_DN11W2", "wave2_numeric_DN12W2", "wave2_numeric_DN13W2")
  
  # Transform the scores for negative atmosphere questions (5 - original score)
  eligible_students <- eligible_students %>%
    mutate(across(all_of(positive_atmosphere_cols_wave2), ~ 5 - ., .names = "transformed_{col}"))
  
  # List of columns representing all survey questions, including the transformed negative atmosphere questions
  all_obs_cols_wave2 <- c("wave2_numeric_DN1W2", "wave2_numeric_DN2W2", "wave2_numeric_DN3W2", "wave2_numeric_DN4W2",
                          "wave2_numeric_DN5W2", "wave2_numeric_DN6W2", "wave2_numeric_DN7W2", "wave2_numeric_DN8W2", 
                          "transformed_wave2_numeric_DN9W2", "transformed_wave2_numeric_DN10W2", "transformed_wave2_numeric_DN11W2", "transformed_wave2_numeric_DN12W2", "transformed_wave2_numeric_DN13W2",
                          "wave2_numeric_DN14W2")
  
  # Calculate the average score across all available questions
  eligible_students <- eligible_students %>%
    rowwise() %>%
    mutate(DN_raw_wave2 = mean(c_across(all_of(all_obs_cols_wave2)), na.rm = TRUE))  # Calculate average, excluding NA values
  
  
  ################### general info: Number of students remaining after filtering
  num_remaining <- nrow(eligible_students)
  
  # Count and display the number of students remove
  num_removed <- total_students - num_remaining
  cat("School", i, ": Total students =", total_students, ", Students removed =", num_removed, ", Students remaining =", num_remaining, "\n")
  
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
    ungroup()  # Remove rowwise grouping after the operation
  
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
  
  A_combined_list[[i]] <- A_combined
  
  # Step 7: Compute Spectral Embedding
  # Determine the number of dimensions to keep based on the eigen gap or other criteria
  num_dimensions_to_keep <- 5  # Example: keep the top 2 dimensions
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
  
  # Store graph information for this school
  graph_summary <- rbind(graph_summary, data.frame(School = i, Graph_Size = graph_size, Sparsity = nu_obv))
  
  # Step 8: Extract and plot the eigenvalues for spectral embedding
  #eig <- eigen(A_combined, symmetric = TRUE)  # Compute eigen decomposition again for all eigenvalues
  #sorted_eigenvalues <- sort(abs(eig$values), decreasing = TRUE)  # Sort eigenvalues in descending order
  # Plot eigenvalues
  #eigen_plot <- ggplot(data.frame(Index = 1:length(sorted_eigenvalues), Eigenvalue = sorted_eigenvalues), aes(x = Index, y = Eigenvalue)) +
  #  geom_point() +
  #  geom_line() +
  #  ggtitle(paste("Eigenvalues for School", i)) +
  #  xlab("Index") +
  #  ylab("Eigenvalue") +
  #  theme_minimal()
  # Save the plot to the folder
  #ggsave(filename = paste0(output_folder, "/Eigenvalues_School_", i, ".png"), plot = eigen_plot)
  
  # Step 9: Extract and plot singular values for truncated SVD
  #svd_result <- svd(A_combined)
  #sorted_singular_values <- sort(svd_result$d, decreasing = TRUE)  # Sort singular values in descending order
  # Plot singular values
  #svd_plot <- ggplot(data.frame(Index = 1:length(sorted_singular_values), SingularValue = sorted_singular_values), aes(x = Index, y = SingularValue)) +
  #  geom_point() +
  #  geom_line() +
  #  ggtitle(paste("Singular Values for School", i)) +
  #  xlab("Index") +
  #  ylab("Singular Value") +
  #  theme_minimal()
  
  # Save the plot to the folder
  #ggsave(filename = paste0(output_folder, "/SingularValues_School_", i, ".png"), plot = svd_plot)
  
  
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
  select(SCHID, GRC, GR, GRC_NUMERIC, RETURN, TREAT, TREAT_BINARY, TREATMENT_FACTOR, NEIGHBORHOOD_AVG_TREATMENT, Race_overall, GENDER, LIVEWB, 
         DN_raw_wave1, DN_raw_wave2, 
         Spectral_1, Spectral_2, Spectral_3, Spectral_4, Spectral_5, 
         SVD_1, SVD_2, SVD_3, SVD_4, SVD_5, 
         OTHETH, 
         rooted_ts.density, node_degree.density,
         starts_with("numeric_DN"),
         starts_with("wave2_numeric_DN"),
         starts_with("transformed_numeric_DN"),
         starts_with("transformed_wave2_numeric_DN"),
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



#####################################################
# concatenate all adjacency into a aggregated graph #
#####################################################
# Construct final block diagonal adjacency matrix for the entire dataset
A_final <- bdiag(A_combined_list)
# Convert back to dense matrix if needed
A_final <- as.matrix(A_final)
# Check dimensions
cat("Final aggregated adjacency matrix size:", dim(A_final), "\n")
# Compute the number of nodes (size of A_final)
num_nodes <- nrow(A_final)
# Compute total number of edges (sum of A_final elements divided by 2 for undirected graphs)
total_edges <- sum(A_final) / 2
# Compute the sparsity parameter (graph density)
sparsity_nu <- total_edges / choose(num_nodes, 2)  # Equivalent to 2 * total_edges / (num_nodes * (num_nodes - 1))
# Print the results
cat("Estimated sparsity parameter (ν):", sparsity_nu, "\n")
