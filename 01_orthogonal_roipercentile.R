################################################################################
# OPTION C: ORTHOGONAL FUNCTIONAL LOCALIZER APPROACH
# MODIFIED: Top 100% Percentile Voxel Selection
#
# This script implements voxel selection using an orthogonal contrast:
# - SELECTION: Based on main effect of task (all conditions > control)
# - TESTING: Condition-specific betas extracted from selected voxels
#
# KEY CHANGE: Selects top 100% of voxels within each ROI (not fixed N)
# This accounts for anatomical variability across ROIs and participants
#
# Contrast mapping:
# con_0001 = G_F (finiteness_error > control)
# con_0002 = G_G (correct_grammar > control)  
# con_0003 = G_P (plurality_error > control)
################################################################################

# Required packages
packages <- c("oro.nifti", "neurobase", "dplyr", "tidyr", "readr", "stringr", "ggplot2")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

################################################################################
# CONFIGURATION
################################################################################

# Directory containing your first-level contrast files
contrasts_dir <- "D:/Paper 1/Contrasts"

# ROI mask directory
roi_dir <- "C:/Users/anama/OneDrive/Documents/Diss1_DMPM_Paper/DMPM/ROIFiles"

# ROI mask filenames
roi_masks <- list(
  hippocampus_left = "resampled_Hippocampus_L.nii",
  hippocampus_right = "resampled_Hippocampus_R.nii",
  parahippocampal_left = "resampled_ParaHippocampal_L.nii",
  parahippocampal_right = "resampled_ParaHippocampal_R.nii",
  caudate_left = "resampled_Caudate_L.nii",
  caudate_right = "resampled_Caudate_R.nii",
  putamen_left = "resampled_Putamen_L.nii",
  putamen_right = "resampled_Putamen_R.nii"
)

# VOXEL SELECTION PARAMETERS - CHANGE BASED ON DESIRED THRESHOLD
percentile_cutoff <- 0.0  # Top 100% = 0th percentile and above, Top 25% = 0.75, Top 50% = 0.50, Top 75% = 0.24
# This means the top 100% of voxels by activation within each ROI

# Output directory
output_dir <- "C:/Users/anama/OneDrive/Documents/Diss1_DMPM_Paper/2025_Results"

# Output files
output_file <- "top100percentile_beta_values.csv"
selection_info_file <- "top100percentile_voxel_selection_info.csv"
qc_report_file <- "top100percentile_extraction_QC_report.csv"

################################################################################
# HELPER FUNCTIONS
################################################################################

# Function to compute main effect of task (localizer contrast)
# This averages activation across all three task conditions
compute_localizer_contrast <- function(con1_img, con2_img, con3_img) {
  # Convert to arrays
  con1_array <- as.vector(con1_img)
  con2_array <- as.vector(con2_img)
  con3_array <- as.vector(con3_img)
  
  # Compute mean across conditions (main effect)
  localizer_array <- (con1_array + con2_array + con3_array) / 3
  
  # Return as nifti with same structure
  localizer_img <- con1_img
  localizer_img@.Data <- array(localizer_array, dim = dim(con1_img))
  
  return(localizer_img)
}

# Function to select TOP PERCENTILE voxels from localizer contrast within ROI
# Modified to use percentile cutoff instead of fixed N voxels
select_top_percentile_voxels <- function(localizer_img, mask_img, percentile = 0.5) {
  
  # Get mask voxels
  mask_vals <- as.vector(mask_img)
  mask_indices <- which(mask_vals > 0)
  
  # Get localizer values within mask
  localizer_vals <- as.vector(localizer_img)
  roi_localizer <- localizer_vals[mask_indices]
  
  # Remove non-finite values
  finite_mask <- is.finite(roi_localizer)
  roi_localizer_finite <- roi_localizer[finite_mask]
  mask_indices_finite <- mask_indices[finite_mask]
  
  # Check if we have any voxels
  n_available <- length(roi_localizer_finite)
  
  if (n_available == 0) {
    return(list(
      voxel_indices = integer(0),
      n_selected = 0,
      n_available = 0,
      selection_threshold = NA,
      mean_selected_value = NA,
      percentile_used = NA
    ))
  }
  
  # Calculate percentile threshold
  threshold <- quantile(roi_localizer_finite, probs = percentile, type = 7)
  
  # Select voxels above percentile threshold
  above_threshold <- roi_localizer_finite >= threshold
  top_indices <- which(above_threshold)
  selected_voxel_indices <- mask_indices_finite[top_indices]
  
  # Summary statistics
  n_selected <- length(selected_voxel_indices)
  mean_value <- mean(roi_localizer_finite[top_indices])
  pct_selected <- 100 * n_selected / n_available
  
  return(list(
    voxel_indices = selected_voxel_indices,
    n_selected = n_selected,
    n_available = n_available,
    selection_threshold = threshold,
    mean_selected_value = mean_value,
    percentile_used = percentile,
    pct_of_roi = pct_selected
  ))
}

# Function to extract mean beta from specific voxels
extract_mean_from_voxels <- function(beta_img, voxel_indices) {
  if (length(voxel_indices) == 0) {
    return(NA)
  }
  
  beta_vals <- as.vector(beta_img)
  selected_betas <- beta_vals[voxel_indices]
  
  # Remove non-finite values
  selected_betas_finite <- selected_betas[is.finite(selected_betas)]
  
  if (length(selected_betas_finite) == 0) {
    return(NA)
  }
  
  return(mean(selected_betas_finite))
}

################################################################################
# STEP 1: LOAD ROI MASKS
################################################################################

cat("Configuration:\n")
cat("Contrast directory:", contrasts_dir, "\n")
cat("ROI mask directory:", roi_dir, "\n")
cat("Voxel percentile cutoff:", percentile_cutoff, "(= top 100%)\n") # Change based on chosen threshold
cat("Output directory:", output_dir, "\n\n")

# Load ROI masks
cat("Loading ROI masks...\n")
masks <- list()

for (roi_name in names(roi_masks)) {
  mask_path <- file.path(roi_dir, roi_masks[[roi_name]])
  if (file.exists(mask_path)) {
    masks[[roi_name]] <- readNIfTI(mask_path, reorient = FALSE)
    n_voxels <- sum(as.vector(masks[[roi_name]]) > 0)
    cat(sprintf("Loaded: %-100s (%d voxels total)\n", roi_name, n_voxels))
  } else {
    warning(paste("NOT FOUND:", roi_name, "at", mask_path))
  }
}

if (length(masks) == 0) {
  stop("\nERROR: No ROI masks loaded!\n")
}

cat(paste("\nSuccessfully loaded", length(masks), "ROI masks\n\n"))

################################################################################
# STEP 2: GET CONTRAST FILES AND ORGANIZE BY PARTICIPANT
################################################################################

cat("Scanning for contrast files...\n")

# Get all contrast files
all_files <- list.files(contrasts_dir, pattern = "G[FGP]C.*con_000[1-3]\\.nii", full.names = TRUE)

if (length(all_files) == 0) {
  stop("\nERROR: No contrast files found!\n")
}

cat(paste("Found", length(all_files), "contrast files\n\n"))

# Parse filenames
file_info <- data.frame(
  filepath = all_files,
  filename = basename(all_files),
  stringsAsFactors = FALSE
) %>%
  mutate(
    condition_prefix = str_extract(filename, "G[FGP]C"),  # Extract G_F, G_G, or G_P
    contrast_num = str_extract(filename, "con_000[1-3]"),
    subject_id = str_extract(filename, "[0-9]{4}"),
    condition = case_when(
      condition_prefix == "GFC" & contrast_num == "con_0001" ~ "finiteness_error",
      condition_prefix == "GGC" & contrast_num == "con_0002" ~ "correct_grammar",
      condition_prefix == "GPC" & contrast_num == "con_0003" ~ "plurality_error",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(condition), !is.na(subject_id))

# Get unique participants
participants <- unique(file_info$subject_id)
n_participants <- length(participants)

cat("Unique participants found:", n_participants, "\n")
cat("Files per participant should be: 3 (one per condition)\n\n")

# Check data completeness
files_per_subj <- file_info %>%
  group_by(subject_id) %>%
  summarise(n_files = n(), .groups = "drop")

cat("Data completeness check:\n")
expected_files <- 3 * length(masks)
cat(sprintf("Expected files per participant: %d\n", expected_files))
print(table(files_per_subj$n_files))

################################################################################
# STEP 3: EXTRACT VOXELS AND BETAS
################################################################################


cat("EXTRACTING VOXELS AND BETA VALUES\n")
cat(strrep("=", 80), "\n\n")

results_list <- list()
selection_info_list <- list()
counter <- 1
pb <- txtProgressBar(min = 0, max = n_participants, style = 3)

for (i in seq_along(participants)) {
  subj_id <- participants[i]
  
  tryCatch({
    # Get this participant's files
    subj_files <- file_info %>%
      filter(subject_id == subj_id) %>%
      arrange(condition)
    
    # Load participant's contrast files and create localizer
    con1_path <- subj_files$filepath[subj_files$condition == "finiteness_error"]
    con2_path <- subj_files$filepath[subj_files$condition == "correct_grammar"]
    con3_path <- subj_files$filepath[subj_files$condition == "plurality_error"]
    
    if (length(con1_path) == 0 || length(con2_path) == 0 || length(con3_path) == 0) {
      warning(paste("Participant", subj_id, "missing one or more contrast files"))
      setTxtProgressBar(pb, i)
      return()
    }
    
    con1_img <- readNIfTI(con1_path, reorient = FALSE)
    con2_img <- readNIfTI(con2_path, reorient = FALSE)
    con3_img <- readNIfTI(con3_path, reorient = FALSE)
    
    # Compute localizer (main effect)
    localizer_img <- compute_localizer_contrast(con1_img, con2_img, con3_img)
    
    # Process each ROI
    for (roi_name in names(masks)) {
      roi_mask <- masks[[roi_name]]
      
      # Select TOP 100% voxels using percentile method
      selection_result <- select_top_percentile_voxels(
        localizer_img,
        roi_mask,
        percentile = percentile_cutoff
      )
      
      # Store selection info
      selection_info_list[[length(selection_info_list) + 1]] <- data.frame(
        participant_id = subj_id,
        roi = roi_name,
        n_selected = selection_result$n_selected,
        n_available = selection_result$n_available,
        pct_of_roi = selection_result$pct_of_roi,
        selection_threshold = selection_result$selection_threshold,
        mean_selected_value = selection_result$mean_selected_value,
        percentile_cutoff = selection_result$percentile_used,
        stringsAsFactors = FALSE
      )
      
      # Extract betas for each condition from selected voxels
      beta_finiteness <- extract_mean_from_voxels(con1_img, selection_result$voxel_indices)
      beta_correct <- extract_mean_from_voxels(con2_img, selection_result$voxel_indices)
      beta_plurality <- extract_mean_from_voxels(con3_img, selection_result$voxel_indices)
      
      # Store results
      # Finiteness error
      results_list[[counter]] <- data.frame(
        participant_id = subj_id,
        condition = "finiteness_error",
        roi = roi_name,
        beta_value = beta_finiteness,
        n_voxels_used = selection_result$n_selected,
        pct_of_roi_used = round(selection_result$pct_of_roi, 1),
        stringsAsFactors = FALSE
      )
      counter <- counter + 1
      
      # Correct grammar
      results_list[[counter]] <- data.frame(
        participant_id = subj_id,
        condition = "correct_grammar",
        roi = roi_name,
        beta_value = beta_correct,
        n_voxels_used = selection_result$n_selected,
        pct_of_roi_used = round(selection_result$pct_of_roi, 1),
        stringsAsFactors = FALSE
      )
      counter <- counter + 1
      
      # Plurality error
      results_list[[counter]] <- data.frame(
        participant_id = subj_id,
        condition = "plurality_error",
        roi = roi_name,
        beta_value = beta_plurality,
        n_voxels_used = selection_result$n_selected,
        pct_of_roi_used = round(selection_result$pct_of_roi, 1),
        stringsAsFactors = FALSE
      )
      counter <- counter + 1
    }
    
  }, error = function(e) {
    warning(paste("Error processing participant", subj_id, ":", e$message))
  })
  
  setTxtProgressBar(pb, i)
}

close(pb)

################################################################################
# STEP 4: COMPILE RESULTS
################################################################################

cat("\n\nCompiling results...\n")

# Combine results
beta_data <- bind_rows(results_list) %>%
  arrange(participant_id, roi, condition)

selection_info <- bind_rows(selection_info_list) %>%
  arrange(participant_id, roi)

cat(sprintf("Total observations: %d\n", nrow(beta_data)))
cat(sprintf("Participants: %d\n", n_distinct(beta_data$participant_id)))
cat(sprintf("ROIs: %d\n", n_distinct(beta_data$roi)))
cat(sprintf("Conditions: %d\n", n_distinct(beta_data$condition)))

################################################################################
# STEP 5: QUALITY CONTROL CHECKS
################################################################################


cat("QUALITY CONTROL\n")
cat(strrep("=", 80), "\n\n")

# Check for missing values
n_missing <- sum(is.na(beta_data$beta_value))
pct_missing <- round(100 * n_missing / nrow(beta_data), 1)

cat(sprintf("Missing beta values: %d / %d (%0.1f%%)\n", 
            n_missing, nrow(beta_data), pct_missing))

# Voxel selection summary (now using percentile-based counts)
cat("\nVoxel Selection Summary (Top 100% Percentile):\n") # Change based on chosen threshold
selection_summary <- selection_info %>%
  summarise(
    mean_selected = round(mean(n_selected, na.rm = TRUE), 1),
    sd_selected = round(sd(n_selected, na.rm = TRUE), 1),
    min_selected = min(n_selected, na.rm = TRUE),
    max_selected = max(n_selected, na.rm = TRUE),
    mean_pct_of_roi = round(mean(pct_of_roi, na.rm = TRUE), 1),
    mean_available = round(mean(n_available, na.rm = TRUE), 1)
  )
print(selection_summary)

# Voxel selection by ROI
cat("\nVoxels Selected by ROI:\n")
selection_by_roi <- selection_info %>%
  group_by(roi) %>%
  summarise(
    mean_selected = round(mean(n_selected, na.rm = TRUE), 1),
    sd_selected = round(sd(n_selected, na.rm = TRUE), 1),
    min_selected = min(n_selected, na.rm = TRUE),
    max_selected = max(n_selected, na.rm = TRUE),
    mean_pct_of_roi = round(mean(pct_of_roi, na.rm = TRUE), 1),
    mean_threshold = round(mean(selection_threshold, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(mean_selected)

print(selection_by_roi)

# Verify that we're getting approximately 100% across the board
cat("\nPercentile Selection Verification (should be ~100%):\n")
pct_check <- selection_info %>%
  summarise(
    min_pct = min(pct_of_roi, na.rm = TRUE),
    max_pct = max(pct_of_roi, na.rm = TRUE),
    mean_pct = round(mean(pct_of_roi, na.rm = TRUE), 1),
    sd_pct = round(sd(pct_of_roi, na.rm = TRUE), 1)
  )
print(pct_check)

################################################################################
# STEP 6: SAVE OUTPUTS
################################################################################

cat("SAVING RESULTS\n")


# Save main beta values
write_csv(beta_data, file.path(output_dir, output_file))
cat(paste("Beta values saved to:", output_file, "\n"))

# Save voxel selection info
write_csv(selection_info, file.path(output_dir, selection_info_file))
cat(paste("Voxel selection info saved to:", selection_info_file, "\n"))

# Create QC report
qc_report <- list(
  extraction_date = as.character(Sys.time()),
  method = "Top 100% Percentile Voxel Selection",
  percentile_cutoff = percentile_cutoff,
  n_participants = n_distinct(beta_data$participant_id),
  n_rois = n_distinct(beta_data$roi),
  n_conditions = n_distinct(beta_data$condition),
  mean_voxels_selected = selection_summary$mean_selected,
  mean_pct_of_roi_selected = selection_summary$mean_pct_of_roi,
  missing_values_n = n_missing,
  missing_values_pct = pct_missing
)

write_csv(as.data.frame(qc_report), file.path(output_dir, qc_report_file))
cat(paste("QC report saved to:", qc_report_file, "\n"))

################################################################################
# STEP 7: DIAGNOSTIC VISUALIZATIONS
################################################################################

cat("\nGenerating diagnostic plots...\n")

# Plot 1: Voxel selection distribution
p1 <- ggplot(selection_info, aes(x = n_selected)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = mean(selection_info$n_selected), linetype = "dashed", 
             color = "red", size = 1) +
  labs(
    title = "Distribution of Voxels Selected (Top 100%)",
    subtitle = paste("Mean:", round(mean(selection_info$n_selected), 1), "voxels"),
    x = "Number of Voxels Selected",
    y = "Count"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "top100pct_voxel_selection_distribution.png"), 
       p1, width = 8, height = 6, dpi = 300)

# Plot 2: Voxel selection by ROI
p2 <- ggplot(selection_by_roi, aes(x = reorder(roi, mean_selected), 
                                   y = mean_selected)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_selected - sd_selected, 
                    ymax = mean_selected + sd_selected),
                width = 0.3) +
  coord_flip() +
  labs(
    title = "Mean Voxels Selected by ROI (Top 100%)",
    subtitle = "Error bars = ±1 SD",
    x = "ROI",
    y = "Mean Number of Voxels"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "top100pct_voxels_by_roi.png"), 
       p2, width = 8, height = 6, dpi = 300)

# Plot 3: Percentage of ROI selected (should cluster around 100%)
p2b <- ggplot(selection_info, aes(x = pct_of_roi)) +
  geom_histogram(bins = 20, fill = "coral", alpha = 0.7) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "red", size = 1) +
  labs(
    title = "Percentage of ROI Selected",
    subtitle = "Should cluster around 100%",
    x = "Percent of ROI (%)",
    y = "Count"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "top100pct_pct_of_roi.png"), 
       p2b, width = 8, height = 6, dpi = 300)

# Plot 4: Beta value distribution by condition
beta_data_clean <- beta_data %>%
  filter(!is.na(beta_value))

p3 <- ggplot(beta_data_clean, aes(x = condition, y = beta_value, 
                                  fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Beta Value Distribution by Condition",
    subtitle = "All ROIs combined, top 100% voxels",
    x = "Condition",
    y = "Beta Value"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "top100pct_beta_distribution.png"), 
       p3, width = 8, height = 6, dpi = 300)

cat("Diagnostic plots saved\n")


cat("SUMMARY:\n")
cat(sprintf("Participants processed: %d\n", n_distinct(beta_data$participant_id)))
cat(sprintf("Total observations: %d\n", nrow(beta_data)))
cat(sprintf("Missing values: %d (%0.1f%%)\n", n_missing, pct_missing))
cat(sprintf("Mean voxels selected per ROI: %0.1f\n", 
            selection_summary$mean_selected))
cat(sprintf("Mean %% of ROI selected: %0.1f%%\n", 
            selection_summary$mean_pct_of_roi))
cat("Files saved to:", output_dir, "\n\n")