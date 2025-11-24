################################################################################
# Merge Option C Orthogonal Localizer Beta Values with Behavioral and Demographic Data
# UPDATED: Uses OPTION C extraction (orthogonal localizer voxel selection)
# SIMPLIFIED: Ignores ROI column in behavioral file
################################################################################

library(dplyr)
library(tidyr)
library(readr)

################################################################################
# CONFIGURATION
################################################################################

# Input files
beta_file <- "C:/Users/anama/OneDrive/Documents/Diss1_DMPM_Paper/2025_Results/top25percentile_beta_values.csv"
behavioral_file <- "C:/Users/anama/OneDrive/Documents/Diss1_DMPM_Paper/Gram_ROIsep.csv"

# Output file
output_dir <- "C:/Users/anama/OneDrive/Documents/Diss1_DMPM_Paper/2025_Results"
output_file <- "merged_data_for_analysis.csv"

################################################################################
# LOAD DATA
################################################################################


cat("Loading orthogonal localizer beta values...\n")
if (!file.exists(beta_file)) {
  stop("ERROR: Option C beta file not found!\n",
       "Expected: ", beta_file, "\n",
       "Make sure you've run the extraction script.\n")
}
beta_data <- read_csv(beta_file, show_col_types = FALSE) %>%
  mutate(participant_id = as.character(participant_id))  # Convert to character

cat("Loading behavioral/demographic data...\n")
if (!file.exists(behavioral_file)) {
  stop("ERROR: Behavioral file not found!\n",
       "Expected: ", behavioral_file, "\n")
}
behav_data <- read_csv(behavioral_file, show_col_types = FALSE)

cat(paste("✓ Beta data:", nrow(beta_data), "observations\n"))
cat(paste("✓ Behavioral data:", nrow(behav_data), "observations\n\n"))

################################################################################
# PARSE ROI NAMES TO GET NETWORK AND HEMISPHERE
################################################################################

cat("Parsing ROI names to extract network and hemisphere...\n")
beta_data <- beta_data %>%
  mutate(
    # Extract network from ROI name
    network = case_when(
      grepl("hippocampus|parahippocampal", roi) ~ "DM",
      grepl("caudate|putamen", roi) ~ "PM",
      TRUE ~ NA_character_
    ),
    # Extract hemisphere from ROI name
    hemisphere = case_when(
      grepl("left", roi) ~ "left",
      grepl("right", roi) ~ "right",
      TRUE ~ NA_character_
    ),
    # Add structure name for reference
    structure = case_when(
      grepl("hippocampus", roi) ~ "hippocampus",
      grepl("parahippocampal", roi) ~ "parahippocampal",
      grepl("caudate", roi) ~ "caudate",
      grepl("putamen", roi) ~ "putamen",
      TRUE ~ NA_character_
    )
  )

cat("Network and hemisphere mapping verified:\n")
print(beta_data %>% distinct(roi, network, hemisphere, structure))
cat("\n")

################################################################################
# PREPARE BEHAVIORAL DATA
################################################################################

cat("Preparing behavioral data...\n")

# Rename columns and create condition labels
# IGNORING the 'roi' column from behavioral file
behav_data <- behav_data %>%
  rename(
    participant_id = sub,
    sex = female,  # Assuming 1 = female, 0 = male
    language_score = celf_sts,
    IQ_score = kbit2,
    condition_code = cond,
    accuracy = acc,
    reaction_time = rt
  ) %>%
  mutate(
    participant_id = as.character(participant_id),
    condition = case_when(
      condition_code == 0 ~ "correct_grammar",
      condition_code == 1 ~ "plurality_error",
      condition_code == 2 ~ "finiteness_error",
      TRUE ~ NA_character_
    )
  ) %>%
  # Remove the roi column from behavioral data - we'll use ROIs from beta data
  select(-any_of(c("roi", "roi_code")))

# Get participant-level variables (should be constant across conditions)
participant_info <- behav_data %>%
  group_by(participant_id) %>%
  summarise(
    sex = first(sex),
    handedness = first(handedness),
    language_score = first(language_score),
    IQ_score = first(IQ_score),
    .groups = "drop"
  )

# Get condition-level behavioral measures (one per participant per condition)
# These will be the same across all ROIs for each participant/condition
behavioral_measures <- behav_data %>%
  group_by(participant_id, condition) %>%
  summarise(
    accuracy = mean(accuracy, na.rm = TRUE),  # Average if multiple trials
    reaction_time = mean(reaction_time, na.rm = TRUE),
    .groups = "drop"
  )

cat(paste("  Unique participants in behavioral data:", 
          n_distinct(behav_data$participant_id), "\n"))
cat(paste("  Unique participants in beta data:", 
          n_distinct(beta_data$participant_id), "\n"))
cat(paste("  Behavioral measures per participant:", 
          nrow(behavioral_measures) / n_distinct(behavioral_measures$participant_id), "\n\n"))

################################################################################
# MERGE DATASETS
################################################################################

cat("Merging datasets...\n")

# Merge strategy:
# 1. Beta data has: participant_id, condition, roi, network, structure, hemisphere, beta_value, n_voxels_used
# 2. Behavioral has: participant_id, condition, accuracy, RT (NOT ROI-specific)
# 3. Demographics: participant_id, sex, handedness, language, IQ
# 
# Each beta observation gets matched with behavioral data for that participant/condition. Same behavioral values will apply to all ROIs.

merged_data <- beta_data %>%
  # Add behavioral measures (matched on participant and condition)
  left_join(
    behavioral_measures,
    by = c("participant_id", "condition")
  ) %>%
  # Add demographic info (matched on participant only)
  left_join(
    participant_info,
    by = "participant_id"
  )

cat(paste("✓ Merged data:", nrow(merged_data), "observations\n\n"))

################################################################################
# DATA VALIDATION
################################################################################

# Check for missing values
cat("Missing values by key column:\n")
missing_summary <- data.frame(
  column = names(merged_data),
  n_missing = colSums(is.na(merged_data)),
  pct_missing = round(100 * colSums(is.na(merged_data)) / nrow(merged_data), 1)
) %>%
  filter(n_missing > 0) %>%
  arrange(desc(n_missing))

if (nrow(missing_summary) > 0) {
  print(missing_summary, row.names = FALSE)
  
  if (any(missing_summary$column %in% c("beta_value", "accuracy", "reaction_time"))) {
    cat("\nWarning: Missing values in key variables!\n")
    cat("  This may indicate incomplete behavioral data for some participants/conditions.\n")
  }
} else {
  cat("✓ No missing values in merged data\n")
}
cat("\n")

# Check participant counts
cat("Participant coverage:\n")
cat(paste("  Participants in merged data:", 
          n_distinct(merged_data$participant_id), "\n"))
cat(paste("  Participants with complete data:", 
          merged_data %>% 
            filter(complete.cases(.)) %>% 
            pull(participant_id) %>% 
            n_distinct(), "\n"))
cat("\n")

# Check observations per participant
cat("Observations per participant (should be 3 conditions × 8 ROIs = 24 per participant):\n")
obs_per_participant <- merged_data %>%
  count(participant_id)

cat(paste("  Mean observations:", round(mean(obs_per_participant$n), 1), "\n"))
cat(paste("  Median:", median(obs_per_participant$n), "\n"))
cat(paste("  Range:", min(obs_per_participant$n), "-", max(obs_per_participant$n), "\n"))
cat(paste("  (Expected: 24 = 3 conditions × 8 ROIs)\n\n"))

# Participants with very low data
low_data_participants <- obs_per_participant %>%
  filter(n < 12) %>%  # Less than 50% of possible observations
  arrange(n)

if (nrow(low_data_participants) > 0) {
  cat("Participants with <12 observations (may want to exclude):\n")
  print(low_data_participants, row.names = FALSE)
  cat("\n")
}

# Check balance across conditions
cat("Observations by condition:\n")
condition_counts <- merged_data %>%
  count(condition) %>%
  arrange(desc(n))
print(condition_counts, row.names = FALSE)
cat("\n")

# Check balance across networks
cat("Observations by network:\n")
network_counts <- merged_data %>%
  count(network) %>%
  arrange(desc(n))
print(network_counts, row.names = FALSE)
cat("\n")

# Check balance across ROIs
cat("Observations by ROI:\n")
roi_counts <- merged_data %>%
  count(roi) %>%
  arrange(desc(n))
print(roi_counts, row.names = FALSE)
cat("\n")

# Check voxel counts
cat("Voxels used (mean per ROI selection):\n")
voxel_summary <- merged_data %>%
  group_by(roi) %>%
  summarise(
    mean_voxels = mean(n_voxels_used, na.rm = TRUE),
    min_voxels = min(n_voxels_used, na.rm = TRUE),
    max_voxels = max(n_voxels_used, na.rm = TRUE)
  ) %>%
  arrange(mean_voxels)
print(voxel_summary, row.names = FALSE)
cat("\n")

################################################################################
# DESCRIPTIVE STATISTICS
################################################################################

cat(strrep("=", 80), "\n")
cat("DESCRIPTIVE STATISTICS\n")
cat(strrep("=", 80), "\n\n")

# Participant demographics
cat("Participant Demographics:\n")
demo_summary <- participant_info %>%
  summarise(
    n_participants = n(),
    n_female = sum(sex == 1, na.rm = TRUE),
    pct_female = round(100 * n_female / n_participants, 1),
    mean_language = round(mean(language_score, na.rm = TRUE), 1),
    sd_language = round(sd(language_score, na.rm = TRUE), 1),
    mean_IQ = round(mean(IQ_score, na.rm = TRUE), 1),
    sd_IQ = round(sd(IQ_score, na.rm = TRUE), 1)
  )
print(demo_summary)
cat("\n")

# Beta values by network and condition
cat("Beta values by network and condition:\n")
beta_summary <- merged_data %>%
  group_by(network, condition) %>%
  summarise(
    n = n(),
    mean_beta = round(mean(beta_value, na.rm = TRUE), 3),
    sd_beta = round(sd(beta_value, na.rm = TRUE), 3),
    min_beta = round(min(beta_value, na.rm = TRUE), 3),
    max_beta = round(max(beta_value, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(network, condition)
print(beta_summary, n = Inf)
cat("\n")

# Accuracy by condition
cat("Accuracy by condition:\n")
acc_summary <- merged_data %>%
  group_by(condition) %>%
  summarise(
    n_obs = n(),
    n_with_acc = sum(!is.na(accuracy)),
    mean_acc = round(mean(accuracy, na.rm = TRUE), 3),
    sd_acc = round(sd(accuracy, na.rm = TRUE), 3),
    .groups = "drop"
  )
print(acc_summary)
cat("\n")

# Reaction time by condition
cat("Reaction time by condition:\n")
rt_summary <- merged_data %>%
  group_by(condition) %>%
  summarise(
    n_obs = n(),
    n_with_rt = sum(!is.na(reaction_time)),
    mean_rt = round(mean(reaction_time, na.rm = TRUE), 0),
    sd_rt = round(sd(reaction_time, na.rm = TRUE), 0),
    .groups = "drop"
  )
print(rt_summary)
cat("\n")

################################################################################
# ROI REPRESENTATION ANALYSIS
################################################################################

cat(strrep("=", 80), "\n")
cat("ROI REPRESENTATION IN FINAL DATASET\n")
cat(strrep("=", 80), "\n\n")

# Calculate coverage by ROI
roi_representation <- merged_data %>%
  count(roi) %>%
  mutate(
    max_possible = n_distinct(merged_data$participant_id) * 3,  # 3 conditions
    coverage_pct = round(100 * n / max_possible, 1)
  ) %>%
  arrange(coverage_pct)

cat("Note: Coverage shows observations per ROI (3 conditions per participant)\n\n")
print(roi_representation, row.names = FALSE)
cat("\n")

################################################################################
# SAVE MERGED DATA
################################################################################

# Save full merged dataset
output_path <- file.path(output_dir, output_file)
write_csv(merged_data, output_path)
cat(paste("✓ Merged data saved to:", output_file, "\n"))
cat(paste("  Location:", output_path, "\n"))
cat(paste("  Observations:", nrow(merged_data), "\n"))
cat(paste("  Participants:", n_distinct(merged_data$participant_id), "\n\n"))

# Create a "complete cases" version (no missing behavioral data)
merged_data_complete <- merged_data %>%
  filter(!is.na(accuracy) & !is.na(reaction_time) & !is.na(beta_value))

if (nrow(merged_data_complete) < nrow(merged_data)) {
  output_file_complete <- gsub("\\.csv$", "_complete.csv", output_file)
  output_path_complete <- file.path(output_dir, output_file_complete)
  write_csv(merged_data_complete, output_path_complete)
  
  cat(paste("✓ Complete cases dataset saved to:", output_file_complete, "\n"))
  cat(paste("  Observations:", nrow(merged_data_complete), "\n"))
  cat(paste("  Excluded:", nrow(merged_data) - nrow(merged_data_complete), 
            "observations with missing data\n\n"))
}


summary_file <- file.path(output_dir, "merge_summary.txt")
sink(summary_file)

cat("Generated:", as.character(Sys.time()), "\n")
cat(strrep("=", 80), "\n\n")

cat("INPUT FILES:\n")
cat("  Beta values:", beta_file, "\n")
cat("  Behavioral:", behavioral_file, "\n\n")

cat("OUTPUT FILES:\n")
cat(" ", output_file, "\n")
cat(" ", output_file_complete, "\n\n")

cat("  Total observations:", nrow(merged_data), "\n")
cat("  Complete cases:", nrow(merged_data_complete), "\n")
cat("  Participants:", n_distinct(merged_data$participant_id), "\n")
cat("  Conditions:", paste(unique(merged_data$condition), collapse = ", "), "\n")
cat("  Networks:", paste(unique(merged_data$network), collapse = ", "), "\n")
cat("  ROIs:", paste(unique(merged_data$roi), collapse = ", "), "\n\n")

cat("MISSING DATA:\n")
if (nrow(missing_summary) > 0) {
  print(missing_summary)
} else {
  cat("  No missing values\n")
}
cat("\n")

cat("DEMOGRAPHICS:\n")
print(demo_summary)
cat("\n")

cat("BETA VALUES BY NETWORK AND CONDITION:\n")
print(beta_summary)
cat("\n")

cat("BEHAVIORAL DATA:\n")
cat("Accuracy by condition:\n")
print(acc_summary)
cat("\nReaction time by condition:\n")
print(rt_summary)
cat("\n")

cat("ROI REPRESENTATION:\n")
print(roi_representation)
cat("\n")

cat("VOXEL SELECTION INFORMATION:\n")
print(voxel_summary)
cat("\n")


sink()
