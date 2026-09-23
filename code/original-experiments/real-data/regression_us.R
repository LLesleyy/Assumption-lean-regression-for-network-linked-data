# Required libraries
library(dplyr)  # For data manipulation
library(tidyr)  # For handling missing values
library(ggplot2) # For plotting

# Correct function to map survey responses to numeric values, fully vectorized
map_survey_responses <- function(x) {
  levels <- c("(0) Almost nobody", "(1) A few people", "(2) About 25%", 
              "(3) About 50%", "(4) About 75%", "(5) Almost everyone")
  numeric_values <- 0:5
  mapped_values <- ifelse(is.na(x), NA, numeric_values[match(x, levels)])  # Use ifelse to handle NAs properly
  return(mapped_values)
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
  
  # Step 1: Filter out students with missing demographic information
  eligible_students <- df %>%
    #filter(ETHO == "(0) Race/ethnicity included") %>%
    filter(!is.na(TREAT) & !is.na(SCHID) & !is.na(GENDER) & 
             !is.na(GRC) & !is.na(RETURN) & !is.na(LIVEWB) & 
             !is.na(ETHW) & !is.na(ETHB) & !is.na(ETHH) & 
             !is.na(ETHA) & !is.na(ETHC) & !is.na(ETHSA))
  
  # Step 2: Convert survey responses from wave 2 to numeric values
  survey_cols_wave2 <- c("PN1W2", "PN2W2", "PN3W2", "PN4W2", "PN5W2", "PN6W2", 
                         "PN7W2", "PN8W2", "PN9W2", "PN10W2", "PN11W2", "PN12W2", "PN13W2")
  
  eligible_students <- eligible_students %>%
    mutate(across(all_of(survey_cols_wave2), ~ map_survey_responses(as.character(.)), .names = "numeric_{col}"))
  
  # Step 3: Filter out students with all NAs in wave 2 survey responses
  eligible_students <- eligible_students %>%
    rowwise() %>%
    filter(sum(is.na(c_across(starts_with("numeric_PN")))) < length(survey_cols_wave2))
  
  # Step 4: Compute the response variable 'y' as the average of available survey responses from wave 2
  eligible_students <- eligible_students %>%
    mutate(y_raw_wave2 = mean(c_across(starts_with("numeric_PN")), na.rm = TRUE))
  
  # Number of students remaining after filtering
  num_remaining <- nrow(eligible_students)
  
  # Count and display the number of students removed
  num_removed <- total_students - num_remaining
  cat("School", i, ": Total students =", total_students, ", Students removed =", num_removed, ", Students remaining =", num_remaining, "\n")
  
  # Apply the transformation: raise to the power of 1.2 for wave 2 data
  if (nrow(eligible_students) > 0) {
    eligible_students <- eligible_students %>%
      mutate(y_wave2 = y_raw_wave2^1.2)
  } else {
    warning(paste("No valid non-NA values for transformation in school", i))
    next
  }
  
  # Step 5: Create a new binary variable TREAT_BINARY
  eligible_students <- eligible_students %>%
    mutate(TREAT_BINARY = ifelse(grepl("Treatment", TREAT), 1, 0))
  
  # Step 6: Construct the combined adjacency matrix A_combined using A1 and A2
  A1 <- school_data[["A1"]]
  A2 <- school_data[["A2"]]
  
  # Use row indices to subset the adjacency matrices
  eligible_indices <- as.integer(rownames(eligible_students))
  A1 <- A1[eligible_indices, eligible_indices, drop = FALSE]
  A2 <- A2[eligible_indices, eligible_indices, drop = FALSE]
  
  # Construct the combined adjacency matrix
  A_combined <- pmax(A1, A2)  # Element-wise maximum to combine A1 and A2
  
  # Compute Local Rooted Two-Star Density for A_combined
  degree_vector <- rowSums(A_combined)
  rooted_two_star_density <- choose(degree_vector, 2) / choose(nrow(A_combined) - 1, 2)
  
  # Estimate Sparsity Level for A_combined
  num_edges <- sum(A_combined) / 2  # since A_combined is assumed to be undirected, divide by 2
  num_possible_edges <- choose(nrow(A_combined), 2)
  sparsity_level <- num_edges / num_possible_edges
  
  # Normalize rooted_two_star_density for A_combined
  normalized_rooted_two_star_density <- rooted_two_star_density / (sparsity_level^2)
  eligible_students$normalized_rooted_two_star_density <- normalized_rooted_two_star_density
  
  
  # Step 7: Calculate rooted triangle density
  rooted_triangles_density <- rowSums((A_combined %*% A_combined) * A_combined) / (2*choose(nrow(A_combined) - 1, 2))  # Calculate rooted triangles density for each node
  normalized_rooted_triangle_density <- rooted_triangles_density / (sparsity_level^3)  # Normalize by sparsity level cubed
  eligible_students$normalized_rooted_triangle_density <- normalized_rooted_triangle_density
  
  # Step 8: Calculate nodewise degree normalized by estimated sparsity level
  normalized_nodewise_degree <- degree_vector / ((nrow(A_combined) - 1)*sparsity_level)
  eligible_students$normalized_nodewise_degree <- normalized_nodewise_degree
  
  # Append processed data to the filtered list for combined analysis
  filtered_data_combined[[i]] <- eligible_students
}

# Combine all eligible students across treated schools for combined analysis
all_eligible_data_combined <- bind_rows(filtered_data_combined)

# Step 8: Check number of eligible treated schools
num_eligible_schools <- length(treated_school_indices)
cat("Number of eligible treated schools:", num_eligible_schools, "\n")

if (num_eligible_schools != 28) {
  stop("The number of eligible treated schools is not equal to 28. Please check the filtering criteria.")
}

# Print the total number of students included in the regression model for combined analysis
total_students_in_regression_combined <- nrow(all_eligible_data_combined)
cat("Total number of students included in the regression model for combined analysis:", total_students_in_regression_combined, "\n")

# Step 9: Calculate and display summary statistics for rooted_two_star_density from combined analysis
for (i in treated_school_indices) {
  school_name <- paste("School", i)
  eligible_students <- filtered_data_combined[[i]]
  
  if (!is.null(eligible_students) && nrow(eligible_students) > 0) {
    # Calculate summary statistics for combined analysis
    summary_stats <- summary(eligible_students$rooted_two_star_density_combined)
    sd_value <- sd(eligible_students$rooted_two_star_density_combined, na.rm = TRUE)
    cat("Summary statistics for normalized rooted_two_star_density_combined in", school_name, ":\n")
    print(summary_stats)
    cat("Standard Deviation:", sd_value, "\n\n")
    
    # Plot Histogram for the new response variable 'y_wave2'
    print(
      ggplot(eligible_students, aes(x = y_wave2)) +
        geom_histogram(bins = 30, fill = "purple", color = "white", alpha = 0.7) +
        ggtitle(paste("Histogram of New Response (y_wave2) for", school_name, "from Combined Analysis")) +
        xlab("Transformed Response (y_wave2)") +
        ylab("Frequency") +
        theme_minimal()
    )
  }
}

# Step 10: Define regressors and prepare for regression for combined analysis
regressors_combined_ts <- c("TREAT_BINARY", "SCHID", "ETHW", "ETHB", "ETHH", "ETHA", 
                         "ETHC", "ETHSA", "GRC", "GENDER", "RETURN", "LIVEWB", 
                         "normalized_rooted_two_star_density")

regressors_combined_ts_2 <- c("TREAT_BINARY", "ETHW", "ETHB", "ETHH", "ETHA", 
                         "ETHC", "ETHSA", "GRC", "GENDER", "RETURN", "LIVEWB", 
                         "normalized_rooted_two_star_density")

regressors_combined_tri_2 <- c("TREAT_BINARY", "ETHW", "ETHB", "ETHH", "ETHA", 
                           "ETHC", "ETHSA", "GRC", "GENDER", "RETURN", "LIVEWB", 
                           "normalized_rooted_triangle_density")
regressors_combined_di_2 <- c("TREAT_BINARY", "ETHW", "ETHB", "ETHH", "ETHA", 
                           "ETHC", "ETHSA", "GRC", "GENDER", "RETURN", "LIVEWB", 
                           "normalized_nodewise_degree")

regressors_combined_3 <- c("TREAT_BINARY", "ETHW", "ETHB", "ETHH", "ETHA", 
                           "ETHC", "ETHSA", "GRC", "GENDER", "RETURN", "LIVEWB")

# Perform Regression for the combined analysis
regression_model_combined_ts <- lm(y_wave2 ~ ., data = all_eligible_data_combined %>% select(all_of(c("y_wave2", regressors_combined_ts))))
summary(regression_model_combined_ts)

regression_model_combined_ts_2 <- lm(y_wave2 ~ ., data = all_eligible_data_combined %>% select(all_of(c("y_wave2", regressors_combined_ts_2))))
summary(regression_model_combined_ts_2)

regression_model_combined_tri_2 <- lm(y_wave2 ~ ., data = all_eligible_data_combined %>% select(all_of(c("y_wave2", regressors_combined_tri_2))))
summary(regression_model_combined_tri_2)

regression_model_combined_di_2 <- lm(y_wave2 ~ ., data = all_eligible_data_combined %>% select(all_of(c("y_wave2", regressors_combined_di_2))))
summary(regression_model_combined_di_2)

regression_model_combined_3 <- lm(y_wave2 ~ ., data = all_eligible_data_combined %>% select(all_of(c("y_wave2", regressors_combined_3))))
summary(regression_model_combined_3)

# variability of subgraph count covariate
hist(all_eligible_data_combined$normalized_rooted_two_star_density)
hist(all_eligible_data_combined$normalized_rooted_triangle_density)
hist(all_eligible_data_combined$normalized_nodewise_degree)


