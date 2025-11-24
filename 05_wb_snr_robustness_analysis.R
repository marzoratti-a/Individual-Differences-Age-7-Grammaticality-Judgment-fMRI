# Loading libraries for analyses
library(lme4) # This is the package that actually fits mixed effects models
library(lmerTest) # This package allows for significance tests
library(performance) # This package allows for computation of the intraclass correlation coefficient
library(summarytools) # This allows for a summary of the dataframe
library(psych) # More summary statistics
library(sjPlot) # For model diagnostics & other useful plots
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(HLMdiag) # For model diagnostics
library(lmtest) # To condut LRT 
library(devtools)
library(imputeTS)
library(dplyr)
library(broom.mixed)

# libraries for plotting
library(lubridate)
library(ggplot2)
library(scales)
library(gridExtra)
library(ggthemes)
library(gridExtra)
library("Hmisc")


######################## WB MODELS + BOOTSTRAPS ###############################

#-------------------------------------------------------------------
# Data import & coding
#-------------------------------------------------------------------
wb_dat <- read.csv("C:/Users/anama/OneDrive/Documents/Diss1_DMPM_Paper/2025_Results/wb_merged_data_for_analysis.csv")

wb_dat <- wb_dat %>%
  mutate(
    condition   = factor(condition,
                         levels = c("correct_grammar", "plurality_error", "finiteness_error"),
                         labels = c("No Error", "Easy Error", "Hard Error")),
    hemisphere  = factor(hemisphere, levels = c("right", "left")),
    participant_id = factor(participant_id),
    sex         = factor(sex, levels = c("0", "1"), labels = c("Male","Female")),
    handedness  = factor(handedness, levels = c("3","4","5")),
    beta_value_cz = scale(beta_value, center = TRUE, scale = TRUE)[,1],
    accuracy_c    = scale(accuracy,   center = TRUE, scale = TRUE)[,1]
  )

wb_dat$handedness <- relevel(wb_dat$handedness, ref = "3")
contrasts(wb_dat$handedness) <- contr.treatment(levels(wb_dat$handedness))
wb_dat$hemisphere <- relevel(wb_dat$hemisphere, ref = "right")
contrasts(wb_dat$hemisphere) <- contr.treatment(levels(wb_dat$hemisphere))
wb_dat$sex <- relevel(wb_dat$sex, ref = "Male")
contrasts(wb_dat$sex) <- contr.treatment(levels(wb_dat$sex))

wb_dat$condition  <- relevel(wb_dat$condition, ref = "No Error")
contrasts(wb_dat$condition) <- contr.treatment(levels(wb_dat$condition))

#-------------------------------------------------------------------
# ICC helpers (unchanged)
#-------------------------------------------------------------------
icc.boot <- function(data, x) {
  irr::icc(data[x, ], model = "twoway", type = "agreement", unit = "single")[[7]]
}

calc.icc <- function(y) {
  icc_results <- performance::icc(y)
  icc_results$ICC_adjusted
}

extract_random_slope <- function(model) {
  as.numeric(VarCorr(model)$participant_id[2, 2])
}

extract_random_intercept <- function(model) {
  as.numeric(VarCorr(model)$participant_id[1, 1])
}

#-------------------------------------------------------------------
# Fit WB models
#-------------------------------------------------------------------

## RQ1: mean_beta unconditional
wbmod1 <- lmer(mean_beta ~ 1 + (1 | participant_id),
               REML = FALSE, data = wb_dat)
summary(wbmod1)

tab_model(wbmod1)

wbicc <- quantile(bootMer(wbmod1, calc.icc, nsim = 100)$t,
                  c(0.025, 0.975), na.rm = TRUE)
cat("\nModel 1 WB ICC (mean_beta):", wbicc, "\n")

## RQ1: mean_beta with condition etc.
wbmod2mean.s <- lmer(
  mean_beta ~ condition + hemisphere + sex + handedness +
    (1 + condition | participant_id),
  data = wb_dat,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
summary(wbmod2mean.s)
tab_model(wbmod2mean.s)

## RQ3: accuracy unconditional
wbmod1acc <- lmer(accuracy ~ 1 + (1 | participant_id),
                  REML = FALSE, data = wb_dat)
summary(wbmod1acc)

wbicc_acc <- quantile(bootMer(wbmod1acc, calc.icc, nsim = 100)$t,
                      c(0.025, 0.975), na.rm = TRUE)
cat("\nModel 1 WB ICC (accuracy):", wbicc_acc, "\n")

## RQ3: accuracy_c ~ mean_beta * condition
wbmod2acc.s <- lmer(
  accuracy_c ~ mean_beta * condition + hemisphere +
    sex + handedness +
    (1 + beta_value_cz | participant_id),
  data = wb_dat,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
summary(wbmod2acc.s)

tab_model(wbmod2acc.s)

#-------------------------------------------------------------------
# OUTPUT DIR
#-------------------------------------------------------------------
output_dir <- "C:/Users/anama/OneDrive/Documents/Diss1_DMPM_Paper/2025_Results"

#-------------------------------------------------------------------
# 1) EXPORT FIXED-EFFECT COEFFICIENTS
#-------------------------------------------------------------------
extract_wb_coefs <- function(model, model_label) {
  broom.mixed::tidy(
    model,
    effects  = "fixed",
    conf.int = TRUE
  ) %>%
    dplyr::select(
      term,
      estimate,
      std.error,
      df,
      statistic,
      p.value,
      conf.low,
      conf.high
    ) %>%
    dplyr::mutate(model = model_label, .before = 1)
}

wb_coefs <- dplyr::bind_rows(
  extract_wb_coefs(wbmod2mean.s, "wbmod2mean.s"),
  extract_wb_coefs(wbmod2acc.s,  "wbmod2acc.s")
)

wb_coef_file <- file.path(output_dir, "wb_models_fixed_effects_coefficients.csv")
write.csv(wb_coefs, wb_coef_file, row.names = FALSE)
cat("Exported WB fixed-effect coefficients to:", wb_coef_file, "\n\n")

#-------------------------------------------------------------------
# 2) EXPORT RANDOM-EFFECT BLUPs (per participant)
#-------------------------------------------------------------------
extract_wb_ranefs <- function(model, model_label) {
  broom.mixed::tidy(
    model,
    effects = "ran_vals"
  ) %>%
    dplyr::select(
      group,   # e.g., "participant_id"
      level,   # cluster ID
      term,    # (Intercept), conditionEasy Error, etc.
      estimate
    ) %>%
    dplyr::mutate(model = model_label, .before = 1)
}

wb_ranefs <- dplyr::bind_rows(
  extract_wb_ranefs(wbmod2mean.s, "wbmod2mean.s"),
  extract_wb_ranefs(wbmod2acc.s,  "wbmod2acc.s")
)

wb_ranefs_file <- file.path(output_dir, "wb_models_random_effects_BLUPs.csv")
write.csv(wb_ranefs, wb_ranefs_file, row.names = FALSE)
cat("Exported WB random-effect BLUPs to:", wb_ranefs_file, "\n\n")

#-------------------------------------------------------------------
# 3) BOOTSTRAP CIs FOR FIXED EFFECTS (N = 1000)
#-------------------------------------------------------------------
set.seed(20251114)  # reproducibility

boot_fixef_ci <- function(mod, model_label, nsim = 1000) {
  bb <- bootMer(
    mod,
    FUN = function(fit) fixef(fit),
    nsim = nsim,
    use.u = TRUE,
    type  = "parametric",
    parallel = "no"
  )
  
  orig     <- fixef(mod)
  boot_mat <- bb$t
  
  ci_mat <- t(apply(
    boot_mat, 2,
    quantile, probs = c(0.025, 0.975), na.rm = TRUE
  ))
  
  data.frame(
    model     = model_label,
    term      = names(orig),
    estimate  = as.numeric(orig),
    boot_mean = apply(boot_mat, 2, mean, na.rm = TRUE),
    boot_se   = apply(boot_mat, 2, sd,   na.rm = TRUE),
    ci_lower  = ci_mat[, 1],
    ci_upper  = ci_mat[, 2],
    row.names = NULL
  )
}

cat("================================================================================\n")
cat("BOOTSTRAP CIs FOR FIXED EFFECTS (WB MODELS, N = 1000)\n")
cat("================================================================================\n\n")

wb_boot_fixef <- dplyr::bind_rows(
  boot_fixef_ci(wbmod2mean.s, "wbmod2mean.s"),
  boot_fixef_ci(wbmod2acc.s,  "wbmod2acc.s")
)

wb_boot_fixef_file <- file.path(output_dir, "wb_bootstrap_fixed_effects_1000.csv")
write.csv(wb_boot_fixef, wb_boot_fixef_file, row.names = FALSE)
cat("Exported WB bootstrap fixed-effect CIs to:", wb_boot_fixef_file, "\n\n")

#-------------------------------------------------------------------
# 4) BOOTSTRAP CIs FOR RANDOM-EFFECT SDs (N = 1000)
#-------------------------------------------------------------------
boot_ranef_var_ci <- function(mod, model_label, nsim = 1000,
                              return_sd = TRUE) {
  # variance components from original model (exclude covariances)
  vc_orig <- as.data.frame(VarCorr(mod))
  vc_vars <- vc_orig[vc_orig$var2 == "" | is.na(vc_orig$var2), ]
  
  # labels like "participant_id:(Intercept)", "participant_id:conditionEasy Error"
  comp_names <- paste0(vc_vars$grp, ":", vc_vars$var1)
  
  cat("Components to extract from", model_label, ":\n")
  print(comp_names)
  cat("\n")
  
  get_vc <- function(fit) {
    vc_curr <- as.data.frame(VarCorr(fit))
    vc_curr <- vc_curr[vc_curr$var2 == "" | is.na(vc_curr$var2), ]
    curr_names <- paste0(vc_curr$grp, ":", vc_curr$var1)
    curr_vals  <- structure(vc_curr$vcov, names = curr_names)
    out <- curr_vals[comp_names]
    names(out) <- comp_names
    out
  }
  
  orig_vals <- vc_vars$vcov
  if (return_sd) orig_vals <- sqrt(pmax(orig_vals, 0))
  names(orig_vals) <- comp_names
  
  cat("Running parametric bootstrap for random effects in", model_label,
      "(N =", nsim, ")\n")
  
  set.seed(20251114)
  bb <- suppressWarnings(
    bootMer(
      mod,
      FUN   = get_vc,
      nsim  = nsim,
      use.u = TRUE,
      type  = "parametric",
      parallel = "no"
    )
  )
  
  boot_mat <- bb$t
  if (return_sd) {
    boot_mat <- sqrt(pmax(boot_mat, 0))
  }
  
  ci_mat <- t(apply(
    boot_mat, 2,
    quantile, probs = c(0.025, 0.975), na.rm = TRUE
  ))
  
  data.frame(
    model      = model_label,
    component  = comp_names,
    estimate   = as.numeric(orig_vals),
    boot_mean  = apply(boot_mat, 2, mean, na.rm = TRUE),
    boot_se    = apply(boot_mat, 2, sd,   na.rm = TRUE),
    ci_lower   = ci_mat[, 1],
    ci_upper   = ci_mat[, 2],
    scale      = if (return_sd) "SD" else "Variance",
    row.names  = NULL
  )
}

wb_boot_re_sd <- dplyr::bind_rows(
  boot_ranef_var_ci(wbmod2mean.s, "wbmod2mean.s", nsim = 1000, return_sd = TRUE),
  boot_ranef_var_ci(wbmod2acc.s,  "wbmod2acc.s",  nsim = 1000, return_sd = TRUE)
)

wb_boot_re_sd_file <- file.path(output_dir, "wb_bootstrap_random_effects_sd_1000.csv")
write.csv(wb_boot_re_sd, wb_boot_re_sd_file, row.names = FALSE)
cat("Exported WB bootstrap random-effect SD CIs to:", wb_boot_re_sd_file, "\n\n")







######################## MRI QC Correlations #####################################################

# Load necessary packages
library("ggpubr")
library("ggplot2")
library(sjPlot)
library(car)
library(corrplot)
library(psych)
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(broom.mixed)

# Install openxlsx if not already installed
if (!require("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
  library("openxlsx")
}


# Importing data for RQ1 (DV= dm/pm, L1= condition, L2= participant_id)
qc_dat = read.csv("C:/Users/anama/OneDrive/Documents/Diss1_DMPM_Paper/qc_Gram.csv")


qc_dat <- qc_dat %>%
  mutate(
    cond_numeric = dplyr::recode(condition,
                                 "correct_grammar" = 0,
                                 "plurality_error" = 1,
                                 "finiteness_error" = 2),
    participant_id_numeric = as.numeric(participant_id)
  )

# Factor coding with explicit levels
qc_dat <- qc_dat %>%
  mutate(
    condition   = factor(condition,
                         levels = c("correct_grammar", "plurality_error", "finiteness_error"),
                         labels = c("No Error", "Easy Error", "Hard Error")),
    hemisphere  = factor(hemisphere, levels = c("right", "left")),
    network     = factor(network,    levels = c("PM", "DM")),
    participant_id = factor(participant_id),
    sex         = factor(sex, levels = c("0", "1"), labels = c("Male","Female")),  
    handedness  = factor(handedness, levels = c("3","4","5")),
    n_voxels_c    = scale(n_voxels_used, center = TRUE, scale = TRUE)[,1],
    
    # First coerce to numeric (if they were character/factor), then scale
    beta_value_num = as.numeric(beta_value),
    acc_num        = as.numeric(accuracy),
    
    beta_value_cz  = scale(beta_value_num, center = TRUE, scale = TRUE)[, 1],
    accuracy_c          = scale(acc_num,       center = TRUE, scale = TRUE)[, 1]
  )


qc_dat$handedness <- relevel(qc_dat$handedness, ref = "3")
contrasts(qc_dat$handedness) <- contr.treatment(levels(qc_dat$handedness))
qc_dat$hemisphere <- relevel(qc_dat$hemisphere, ref = "right")
contrasts(qc_dat$hemisphere) <- contr.treatment(levels(qc_dat$hemisphere))
qc_dat$sex <- relevel(qc_dat$sex, ref = "Male")
contrasts(qc_dat$sex) <- contr.treatment(levels(qc_dat$sex))

#enforce treatment coding for condition (No Error = reference)
qc_dat$condition  <- relevel(qc_dat$condition, ref = "No Error")
contrasts(qc_dat$condition) <- contr.treatment(levels(qc_dat$condition))


# List of variables to standardize
vars_to_standardize <- c( "num_repaired","efc","gcor","aor","aqi","snr","tsnr","fber")

# Loop to standardize variables
for (var in vars_to_standardize) {
  qc_dat[[var]] <- as.numeric(scale(qc_dat[[var]]))
}

# Define the ICC function
icc.boot <- function(data, x) {
  irr::icc(data[x, ], model = "twoway", type = "agreement", unit = "single")[[7]]
}

#Create function to calculate ICC from fitted model
calc.icc <- function(y) {
  icc_results<- performance::icc(y)
  icc_results$ICC_adjusted
}

# Define a function to extract the random slope term variance component
extract_random_slope <- function(model) {
  as.numeric(VarCorr(model)$participant_id[2, 2])
}

# Define a function to extract the random intercept term variance component
extract_random_intercept <- function(model) {
  as.numeric(VarCorr(model)$participant_id[1, 1])
}


# Set the output directory and file name
output_dir  <- "C:/Users/anama/OneDrive/Documents/Diss1_DMPM_Paper/2025_Results"
output_file <- file.path(output_dir, "111425_QC_correlations.xlsx")

# Function to perform Pearson correlation test and extract results
perform_correlation <- function(var1, var2, data) {
  x <- data[[var1]]
  y <- data[[var2]]
  
  # Check that both variables are numeric
  if (!is.numeric(x) | !is.numeric(y)) {
    return(list(r = NA, p_value = NA, se = NA, valid = FALSE))
  }
  
  test <- cor.test(x, y, method = "pearson")
  
  # Calculate standard error: SE = sqrt((1 - r^2) / (n - 2))
  r <- test$estimate
  n <- length(x[!is.na(x) & !is.na(y)])
  se <- sqrt((1 - r^2) / (n - 2))
  
  list(
    r = r,
    p_value = test$p.value,
    se = se,
    valid = TRUE
  )
}

data_pm_qc <- droplevels(dplyr::filter(qc_dat, network == "PM"))
cat("PM data: N =", nrow(data_pm_qc), "observations\n\n")

data_dm_qc <- droplevels(dplyr::filter(qc_dat, network == "DM"))
cat("DM data: N =", nrow(data_dm_qc), "observations\n\n")

# Function to run a set of correlations on a given data frame
run_correlation_block <- function(data, tests) {
  res_list <- lapply(tests, function(tt) {
    out <- perform_correlation(tt$var1, tt$var2, data)
    data.frame(
      Variable_1 = tt$var1,
      Variable_2 = tt$var2,
      r          = out$r,
      p_value    = out$p_value,
      SE         = out$se,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res_list)
}

# 1) PM-only correlations: choose an actual DV (e.g., beta_value_cz) vs QC metrics
tests_pm <- list(
  list(var1 = "beta_value_cz", var2 = "num_repaired"),
  list(var1 = "beta_value_cz", var2 = "aor"),
  list(var1 = "beta_value_cz", var2 = "snr"),
  list(var1 = "beta_value_cz", var2 = "fber")
)

# 2) DM-only correlations
tests_dm <- list(
  list(var1 = "beta_value_cz", var2 = "num_repaired"),
  list(var1 = "beta_value_cz", var2 = "aor"),
  list(var1 = "beta_value_cz", var2 = "snr"),
  list(var1 = "beta_value_cz", var2 = "fber")
)

# 3) Global correlations for cond_numeric / participant_id_numeric vs QC
tests_global <- list(
  list(var1 = "cond_numeric",           var2 = "num_repaired"),
  list(var1 = "cond_numeric",           var2 = "aor"),
  list(var1 = "cond_numeric",           var2 = "snr"),
  list(var1 = "cond_numeric",           var2 = "fber"),
  list(var1 = "participant_id_numeric", var2 = "num_repaired"),
  list(var1 = "participant_id_numeric", var2 = "aor"),
  list(var1 = "participant_id_numeric", var2 = "snr"),
  list(var1 = "participant_id_numeric", var2 = "fber")
)

# Run blocks
results_pm     <- run_correlation_block(data_pm_qc, tests_pm)
results_dm     <- run_correlation_block(data_dm_qc, tests_dm)
results_global <- run_correlation_block(qc_dat,    tests_global)

results_pm$Block     <- "PM"
results_dm$Block     <- "DM"
results_global$Block <- "Global"

results_all <- rbind(results_pm, results_dm, results_global)

names(results_all) <- c("Variable 1","Variable 2","Correlation (r)",
                        "P-value","Standard Error","Block")

wb <- createWorkbook()
addWorksheet(wb, "Correlations")
writeData(wb, "Correlations", results_all, startRow = 1, startCol = 1)
saveWorkbook(wb, output_file, overwrite = TRUE)


### SNR ##########################
cat("\nModel 2 snr\n")
cat("\nPM\n")
mod2pm_snr <- lmerTest::lmer(
  beta_value_cz ~ condition + snr + n_voxels_c + sex + handedness +
    (1 + condition | participant_id) + (1 | roi),
  data = data_pm_qc,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
print(summary(mod2pm_snr))
tab_model(mod2pm_snr)

cat("\nDM\n")
mod2dm_snr  <- lmerTest::lmer(
  beta_value_cz ~ condition + snr + n_voxels_c + hemisphere + sex + handedness +
    (1 + condition | participant_id) + (1 | roi),
  data = data_dm_qc,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
print(summary(mod2dm_snr))
tab_model(mod2dm_snr)


cat("\nModel 4 snr\n")
cat("\nPM\n")
mod4pm_snr  <- lmerTest::lmer(
  accuracy_c ~ beta_value_cz * condition + snr + hemisphere + n_voxels_c +
    sex + handedness +
    (1 + beta_value_cz | participant_id),
  data = data_pm_qc,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

print(summary(mod4pm_snr))
tab_model(mod4pm_snr)

cat("\nDM\n")
mod4dm_snr  <- lmerTest::lmer(
  accuracy_c ~ beta_value_cz * condition + snr + hemisphere + n_voxels_c +
    sex + handedness +
    (1 + beta_value_cz | participant_id),
  data = data_dm_qc,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

print(summary(mod4dm_snr))
tab_model(mod4dm_snr)


### fber ##########################
cat("\nModel 2 fber\n")
cat("\nPM\n")
mod2pm_fber <- lmerTest::lmer(
  beta_value_cz ~ condition + fber + n_voxels_c + sex + handedness +
    (1 + condition | participant_id) + (1 | roi),
  data = data_pm_qc,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
print(summary(mod2pm_fber))
tab_model(mod2pm_fber)

cat("\nDM\n")
mod2dm_fber  <- lmerTest::lmer(
  beta_value_cz ~ condition + fber + n_voxels_c + hemisphere + sex + handedness +
    (1 + condition | participant_id) + (1 | roi),
  data = data_dm_qc,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
print(summary(mod2dm_fber))
tab_model(mod2dm_fber)


cat("\nModel 4 fber\n")
cat("\nPM\n")
mod4pm_fber  <- lmerTest::lmer(
  accuracy_c ~ beta_value_cz * condition + fber + hemisphere + n_voxels_c +
    sex + handedness +
    (1 + beta_value_cz | participant_id),
  data = data_pm_qc,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

print(summary(mod4pm_fber))
tab_model(mod4pm_fber)

cat("\nDM\n")
mod4dm_fber  <- lmerTest::lmer(
  accuracy_c ~ beta_value_cz * condition + fber + hemisphere + n_voxels_c +
    sex + handedness +
    (1 + beta_value_cz | participant_id),
  data = data_dm_qc,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

print(summary(mod4dm_fber))
tab_model(mod4dm_fber)


extract_qc_varcorr <- function(model, model_label) {
  vc <- as.data.frame(VarCorr(model))
  vc %>%
    dplyr::mutate(
      model = model_label,
      .before = 1
    )
}


qc_varcorr <- dplyr::bind_rows(
  extract_qc_varcorr(mod2pm_snr,  "mod2pm_snr"),
  extract_qc_varcorr(mod2dm_snr,  "mod2dm_snr"),
  extract_qc_varcorr(mod4pm_snr,  "mod4pm_snr"),
  extract_qc_varcorr(mod4dm_snr,  "mod4dm_snr"),
  extract_qc_varcorr(mod2pm_fber, "mod2pm_fber"),
  extract_qc_varcorr(mod2dm_fber, "mod2dm_fber"),
  extract_qc_varcorr(mod4pm_fber, "mod4pm_fber"),
  extract_qc_varcorr(mod4dm_fber, "mod4dm_fber")
)

qc_varcorr_file <- file.path(output_dir, "qc_models_random_effects_varcorr.csv")
write.csv(qc_varcorr, qc_varcorr_file, row.names = FALSE)
cat("Exported QC random-effect variance components to:", qc_varcorr_file, "\n\n")


################################################################################
# EXPORT FIXED-EFFECT COEFFICIENTS FOR QC MODELS
################################################################################

library(broom.mixed)
library(dplyr)

extract_qc_coefs <- function(model, model_label) {
  broom.mixed::tidy(
    model,
    effects  = "fixed",
    conf.int = TRUE
  ) %>%
    select(
      term,
      estimate,
      std.error,
      df,
      statistic,
      p.value,
      conf.low,
      conf.high
    ) %>%
    mutate(model = model_label, .before = 1)
}

qc_coefs <- bind_rows(
  extract_qc_coefs(mod2pm_snr,  "mod2pm_snr"),
  extract_qc_coefs(mod2dm_snr,  "mod2dm_snr"),
  extract_qc_coefs(mod4pm_snr,  "mod4pm_snr"),
  extract_qc_coefs(mod4dm_snr,  "mod4dm_snr"),
  extract_qc_coefs(mod2pm_fber, "mod2pm_fber"),
  extract_qc_coefs(mod2dm_fber, "mod2dm_fber"),
  extract_qc_coefs(mod4pm_fber, "mod4pm_fber"),
  extract_qc_coefs(mod4dm_fber, "mod4dm_fber")
)

qc_coef_file <- file.path(output_dir, "qc_models_fixed_effects_coefficients.csv")
write.csv(qc_coefs, qc_coef_file, row.names = FALSE)
cat("Exported QC model coefficients to:", qc_coef_file, "\n\n")

################################################################################
# BOOTSTRAP CIs FOR FIXED EFFECTS (QC MODELS, N = 1000)
################################################################################

set.seed(20251114)  # for reproducibility

boot_fixef_ci <- function(mod, model_label, nsim = 1000) {
  bb <- bootMer(
    mod,
    FUN = function(fit) fixef(fit),
    nsim = nsim,
    use.u = TRUE,
    type  = "parametric",
    parallel = "no"
  )
  
  orig     <- fixef(mod)
  boot_mat <- bb$t
  
  ci_mat <- t(apply(
    boot_mat, 2,
    quantile, probs = c(0.025, 0.975), na.rm = TRUE
  ))
  
  data.frame(
    model     = model_label,
    term      = names(orig),
    estimate  = as.numeric(orig),
    boot_mean = apply(boot_mat, 2, mean, na.rm = TRUE),
    boot_se   = apply(boot_mat, 2, sd,   na.rm = TRUE),
    ci_lower  = ci_mat[, 1],
    ci_upper  = ci_mat[, 2],
    row.names = NULL
  )
}

cat("================================================================================\n")
cat("BOOTSTRAP CIs FOR FIXED EFFECTS (QC MODELS, N = 1000)\n")
cat("================================================================================\n\n")

boot_qc_fixef <- dplyr::bind_rows(
  boot_fixef_ci(mod2pm_snr,  "mod2pm_snr"),
  boot_fixef_ci(mod2dm_snr,  "mod2dm_snr"),
  boot_fixef_ci(mod4pm_snr,  "mod4pm_snr"),
  boot_fixef_ci(mod4dm_snr,  "mod4dm_snr"),
  boot_fixef_ci(mod2pm_fber, "mod2pm_fber"),
  boot_fixef_ci(mod2dm_fber, "mod2dm_fber"),
  boot_fixef_ci(mod4pm_fber, "mod4pm_fber"),
  boot_fixef_ci(mod4dm_fber, "mod4dm_fber")
)

qc_boot_fixef_file <- file.path(output_dir, "qc_bootstrap_fixed_effects_1000.csv")
write.csv(boot_qc_fixef, qc_boot_fixef_file, row.names = FALSE)
cat("Exported QC bootstrap fixed-effect CIs to:", qc_boot_fixef_file, "\n\n")


################################################################################
# BOOTSTRAP CIs FOR RANDOM EFFECTS (QC MODELS, N = 1000; SD SCALE)
################################################################################

boot_ranef_var_ci <- function(mod, model_label, nsim = 1000,
                              return_sd = FALSE) {
  # Get variance components from original model (exclude covariances)
  vc_orig <- as.data.frame(VarCorr(mod))
  vc_vars <- vc_orig[vc_orig$var2 == "" | is.na(vc_orig$var2), ]
  
  # Create component labels: e.g. "participant_id:(Intercept)", "participant_id:conditionEasy Error"
  comp_names <- paste0(vc_vars$grp, ":", vc_vars$var1)
  
  cat("Components to extract from", model_label, ":\n")
  print(comp_names)
  cat("\n")
  
  # Function to extract variance components from a single model fit
  get_vc <- function(fit) {
    vc_curr <- as.data.frame(VarCorr(fit))
    vc_curr <- vc_curr[vc_curr$var2 == "" | is.na(vc_curr$var2), ]
    
    # Create a named vector for matching
    curr_names <- paste0(vc_curr$grp, ":", vc_curr$var1)
    curr_vals <- structure(
      vc_curr$vcov,
      names = curr_names
    )
    
    # Extract in the same order as the template, filling missing with NA
    out <- curr_vals[comp_names]
    names(out) <- comp_names
    out
  }
  
  # Get original variance/SD values
  orig_vals <- vc_vars$vcov
  if (return_sd) orig_vals <- sqrt(pmax(orig_vals, 0))
  names(orig_vals) <- comp_names
  
  # Parametric bootstrap with seed for reproducibility
  cat("Running parametric bootstrap for random effects (N =", nsim, ")...\n")
  set.seed(20251114)
  
  bb <- suppressWarnings(
    bootMer(
      mod,
      FUN   = get_vc,
      nsim  = nsim,
      use.u = TRUE,
      type  = "parametric",
      parallel = "no"
    )
  )
  
  boot_mat <- bb$t  # nsim x n_components matrix
  
  # Convert to SD if requested
  if (return_sd) {
    boot_mat <- sqrt(pmax(boot_mat, 0))
  }
  
  # Calculate 95% CIs
  ci_mat <- t(apply(
    boot_mat, 2,
    quantile, probs = c(0.025, 0.975), na.rm = TRUE
  ))
  
  # Build output data frame
  out <- data.frame(
    model      = model_label,
    component  = comp_names,
    estimate   = as.numeric(orig_vals),
    boot_mean  = apply(boot_mat, 2, mean, na.rm = TRUE),
    boot_se    = apply(boot_mat, 2, sd,   na.rm = TRUE),
    ci_lower   = ci_mat[, 1],
    ci_upper   = ci_mat[, 2],
    scale      = if (return_sd) "SD" else "Variance",
    row.names  = NULL
  )
  
  cat("Bootstrap completed for", model_label, "\n\n")
  out
}


boot_qc_re_sd <- dplyr::bind_rows(
  boot_ranef_var_ci(mod2pm_snr,  "mod2pm_snr",  nsim = 1000, return_sd = TRUE),
  boot_ranef_var_ci(mod2dm_snr,  "mod2dm_snr",  nsim = 1000, return_sd = TRUE),
  boot_ranef_var_ci(mod4pm_snr,  "mod4pm_snr",  nsim = 1000, return_sd = TRUE),
  boot_ranef_var_ci(mod4dm_snr,  "mod4dm_snr",  nsim = 1000, return_sd = TRUE),
  boot_ranef_var_ci(mod2pm_fber, "mod2pm_fber", nsim = 1000, return_sd = TRUE),
  boot_ranef_var_ci(mod2dm_fber, "mod2dm_fber", nsim = 1000, return_sd = TRUE),
  boot_ranef_var_ci(mod4pm_fber, "mod4pm_fber", nsim = 1000, return_sd = TRUE),
  boot_ranef_var_ci(mod4dm_fber, "mod4dm_fber", nsim = 1000, return_sd = TRUE)
)

qc_boot_re_sd_file <- file.path(output_dir, "qc_bootstrap_random_effects_sd_1000.csv")
write.csv(boot_qc_re_sd, qc_boot_re_sd_file, row.names = FALSE)
cat("Exported QC bootstrap random-effect SD CIs to:", qc_boot_re_sd_file, "\n\n")


################################################################################
# EXPORT MODEL DIAGNOSTICS: Random Effects, ICC, and REML Criterion
################################################################################

library(lme4)
library(lmerTest)
library(performance)
library(dplyr)
library(broom.mixed)
library(openxlsx)

# Set output directory
output_dir <- "C:/Users/anama/OneDrive/Documents/Diss1_DMPM_Paper/2025_Results"
output_file <- file.path(output_dir, "Model_Diagnostics_RandomEffects_ICC_REML.xlsx")

################################################################################
# FUNCTION: Extract Random Effects, ICC, and REML Criterion
################################################################################

extract_model_diagnostics <- function(model, model_label) {
  # Get variance components
  vc_df <- as.data.frame(VarCorr(model))
  vc_vars <- vc_df[vc_df$var2 == "" | is.na(vc_df$var2), ]
  
  # Extract components
  components_list <- list()
  for (i in seq_len(nrow(vc_vars))) {
    comp_name <- paste0(vc_vars$grp[i], ":", vc_vars$var1[i])
    components_list[[comp_name]] <- vc_vars$vcov[i]
  }
  
  # Calculate ICC
  icc_val <- tryCatch(
    performance::icc(model)$ICC_adjusted,
    error = function(e) NA_real_
  )
  
  # Extract REML criterion (-2 * logLik)
  reml_criterion <- -2 * as.numeric(logLik(model))
  
  # Extract residual variance
  residual_var <- sigma(model)^2
  
  # Build output: one row per component
  output_rows <- list()
  
  if (length(components_list) > 0) {
    for (j in seq_len(length(components_list))) {
      comp_name <- names(components_list)[j]
      comp_var <- components_list[[j]]
      comp_sd <- sqrt(pmax(comp_var, 0))
      
      output_rows[[j]] <- data.frame(
        Model = model_label,
        Component = comp_name,
        Variance = comp_var,
        SD = comp_sd,
        ICC = icc_val,
        REML_Criterion = reml_criterion,
        Residual_Variance = residual_var,
        LogLik = as.numeric(logLik(model)),
        AIC = AIC(model),
        BIC = BIC(model),
        stringsAsFactors = FALSE
      )
    }
  } else {
    # No random components (shouldn't happen, but handle edge case)
    output_rows[[1]] <- data.frame(
      Model = model_label,
      Component = NA_character_,
      Variance = NA_real_,
      SD = NA_real_,
      ICC = icc_val,
      REML_Criterion = reml_criterion,
      Residual_Variance = residual_var,
      LogLik = as.numeric(logLik(model)),
      AIC = AIC(model),
      BIC = BIC(model),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, output_rows)
}

cat("Extracting diagnostics from WB models...\n")

wb_diagnostics <- bind_rows(
  extract_model_diagnostics(wbmod2mean.s, "wbmod2mean.s"),
  extract_model_diagnostics(wbmod2acc.s,  "wbmod2acc.s")
)

cat("Extracting diagnostics from QC models...\n")

qc_diagnostics <- bind_rows(
  extract_model_diagnostics(mod2pm_snr,  "mod2pm_snr"),
  extract_model_diagnostics(mod2dm_snr,  "mod2dm_snr"),
  extract_model_diagnostics(mod4pm_snr,  "mod4pm_snr"),
  extract_model_diagnostics(mod4dm_snr,  "mod4dm_snr"),
  extract_model_diagnostics(mod2pm_fber, "mod2pm_fber"),
  extract_model_diagnostics(mod2dm_fber, "mod2dm_fber"),
  extract_model_diagnostics(mod4pm_fber, "mod4pm_fber"),
  extract_model_diagnostics(mod4dm_fber, "mod4dm_fber")
)

################################################################################
# WRITE TO EXCEL WORKBOOK
################################################################################

wb <- createWorkbook()

# Sheet 1: Summary
addWorksheet(wb, "Summary")
summary_text <- tibble(
  Note = "Model Diagnostics: Random Effects Variance/SD, ICC, and REML Criterion",
  Generated = as.character(Sys.time()),
  Description = "Random effects components are listed per model. ICC and REML criterion are repeated for reference."
)
writeData(wb, "Summary", summary_text)

# Sheet 2: WB Models
addWorksheet(wb, "WB Models")
writeData(wb, "WB Models", wb_diagnostics, startRow = 1, colNames = TRUE)

# Sheet 3: QC Models
addWorksheet(wb, "QC Models")
writeData(wb, "QC Models", qc_diagnostics, startRow = 1, colNames = TRUE)

# Save workbook
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("Output saved to:", output_file, "\n\n")

cat("WB Models Summary:\n")
print(wb_diagnostics)

cat("\n\nQC Models Summary:\n")
print(qc_diagnostics)
