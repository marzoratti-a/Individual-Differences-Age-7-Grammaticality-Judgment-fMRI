################################################################################
# Grammar × Memory Networks Analysis - ORTHOGONAL LOCALIZER

################################################################################

library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(broom.mixed)

# NOTE: We intentionally DO NOT set global sum-to-zero contrasts; we use
# treatment coding for planned contrasts and locally switch to sum-coding
# only when printing supplemental Type-III ANOVAs.

cat("================================================================================\n")
cat("LOADING DATA & PREPARING ANALYSIS\n")
cat("================================================================================\n\n")

################################################################################
# LOAD & PREPARE DATA
################################################################################

output_dir <- "C:/Users/anama/OneDrive/Documents/Diss1_DMPM_Paper/2025_Results/1107_Results/"
merged_file <- "C:/Users/anama/OneDrive/Documents/Diss1_DMPM_Paper/2025_Results/merged_data_for_analysis.csv"

if (!file.exists(merged_file)) {
  stop("ERROR: Merged data file not found!\n")
}

data <- read.csv(merged_file)

cat("Loaded dataset\n")
cat("  Participants:", dplyr::n_distinct(data$participant_id), "\n")
cat("  Observations:", nrow(data), "\n\n")

# Voxel filtering
voxel_minimum <- 20
data_before <- nrow(data)
data <- data %>%
  filter(n_voxels_used >= voxel_minimum)
data_after <- nrow(data)
data <- droplevels(data)  #### CHANGED: drop unused levels after filtering

cat("Voxel filtering (minimum", voxel_minimum, "):\n")
cat("  Retained:", data_after, "/", data_before, "observations\n\n")

# Factor coding with explicit levels
data <- data %>%
  mutate(
    condition   = factor(condition,
                         levels = c("correct_grammar", "plurality_error", "finiteness_error"),
                         labels = c("No Error", "Easy Error", "Hard Error")),
    network     = factor(network,    levels = c("PM", "DM")),
    hemisphere  = factor(hemisphere, levels = c("right", "left")),
    participant_id = factor(participant_id),
    roi         = factor(roi),
    sex         = factor(sex, levels = c("0", "1"), labels = c("Male","Female")),  
    handedness  = factor(handedness, levels = c("3","4","5")),
    # Standardize continuous
    beta_value_cz = scale(beta_value,    center = TRUE, scale = TRUE)[,1],
    n_voxels_c    = scale(n_voxels_used, center = TRUE, scale = TRUE)[,1],
    accuracy_c    = scale(accuracy,      center = TRUE, scale = TRUE)[,1]
  )

data$handedness <- relevel(data$handedness, ref = "3")
contrasts(data$handedness) <- contr.treatment(levels(data$handedness))
data$hemisphere <- relevel(data$hemisphere, ref = "right")
contrasts(data$hemisphere) <- contr.treatment(levels(data$hemisphere))
data$network <- relevel(data$network, ref = "PM")
contrasts(data$network) <- contr.treatment(levels(data$network))
data$sex <- relevel(data$sex, ref = "Male")
contrasts(data$sex) <- contr.treatment(levels(data$sex))

#enforce treatment coding for condition (No Error = reference)
data$condition  <- relevel(data$condition, ref = "No Error")
contrasts(data$condition) <- contr.treatment(levels(data$condition))

stopifnot(all(levels(data$condition) == c("No Error","Easy Error","Hard Error")))
print(contrasts(data$condition))
colnames(model.matrix(~ condition * network + hemisphere + n_voxels_c + sex + handedness, data))[1:20]

cat("Data preparation complete\n")

################################################################################
# HELPER FUNCTIONS
################################################################################
refit_REML_with_p <- function(mod_ml, data = NULL,
                              control = lmerControl(optimizer = "bobyqa",
                                                    optCtrl = list(maxfun = 2e5))) {
  if (is.null(data)) data <- model.frame(mod_ml)  # pull the original frame
  lmerTest::lmer(formula(mod_ml), data = data, REML = TRUE, control = control)
}

# For estimating ICCs

get_vc <- function(mod) {
  vc <- as.data.frame(VarCorr(mod))
  sigma2 <- sigma(mod)^2
  list(vc = vc, sigma2 = sigma2)
}

icc_crossed_random_intercepts <- function(mod, g_part = "participant_id", g_roi = "roi") {
  x <- get_vc(mod); vc <- x$vc; sig2 <- x$sigma2
  v_p   <- vc$vcov[vc$grp == g_part & vc$var1 == "(Intercept)"][1]
  v_roi <- vc$vcov[vc$grp == g_roi   & vc$var1 == "(Intercept)"][1]
  den   <- v_p + v_roi + sig2
  data.frame(
    Group = c("participant", "roi", "total_cluster"),
    ICC   = c(v_p/den, v_roi/den, (v_p + v_roi)/den)
  )
}

# Helper: fit a temporary sum-coded version and print Type-III ANOVA
anova_type3_sum <- function(formula, data, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))) {
  dsum <- data
  # apply sum coding per-column so we don't rely on global options
  if (is.factor(dsum$condition))   contrasts(dsum$condition)   <- contr.sum(length(levels(dsum$condition)))
  if (is.factor(dsum$network))     contrasts(dsum$network)     <- contr.sum(length(levels(dsum$network)))
  if (is.factor(dsum$hemisphere))  contrasts(dsum$hemisphere)  <- contr.sum(length(levels(dsum$hemisphere)))
  if (is.factor(dsum$sex))         contrasts(dsum$sex)         <- contr.sum(length(levels(dsum$sex)))
  if (is.factor(dsum$handedness))  contrasts(dsum$handedness)  <- contr.sum(length(levels(dsum$handedness)))
  m_sum <- lmerTest::lmer(formula, data = dsum, REML = FALSE, control = control)
  print(lmerTest::anova(m_sum, type = "III"))
}

cat("\nContrasts check (condition):\n"); print(contrasts(data$condition))

# Supplemental Type-III ANOVA under sum coding
anova_type3_sum_df <- function(formula, data,
                               fixed_factors = NULL,
                               control = lmerControl(optimizer = "bobyqa",
                                                     optCtrl = list(maxfun = 2e5))) {
  dsum <- droplevels(data)
  
  # If caller doesn't specify, discover candidate fixed factors from the formula
  if (is.null(fixed_factors)) {
    # variables that appear symbolically in the fixed-effects part
    vars_in_form <- all.vars(delete.response(terms(formula)))
    fixed_factors <- intersect(vars_in_form, names(dsum))
  }
  
  # Apply contr.sum ONLY to factors used in the model AND having >= 2 levels
  for (v in fixed_factors) {
    if (!is.null(dsum[[v]]) && is.factor(dsum[[v]])) {
      k <- nlevels(dsum[[v]])
      if (k >= 2) {
        contrasts(dsum[[v]]) <- contr.sum(k)
      }
    }
  }
  # Fit and return Type-III ANOVA table
  m_sum <- lmerTest::lmer(formula, data = dsum, REML = FALSE, control = control)
  out <- as.data.frame(anova(m_sum, type = 3))
  out$Effect <- rownames(out); rownames(out) <- NULL
  out <- out[, c("Effect", setdiff(names(out), "Effect"))]
  out
}

summ_file <- paste0(output_dir, "model_summaries.txt")
sink(summ_file)  # START capturing

cat("Analysis started\n")
cat("", file = summ_file)  # truncate/create

append_summary <- function(obj, title, anova_type3 = FALSE) {
  cat("\n============================================================\n", file = summ_file, append = TRUE)
  cat(title, "\n", file = summ_file, append = TRUE)
  cat("============================================================\n", file = summ_file, append = TRUE)
  capture.output(summary(obj), file = summ_file, append = TRUE)
  if (anova_type3) {
    capture.output({
      local({
        op <- options(contrasts = c("contr.sum","contr.poly"))
        on.exit(options(op), add = TRUE)
        print(anova(obj, type = "III"))
      })
    }, file = summ_file, append = TRUE)
  }
}


cat("================================================================================\n")
cat("RQ1 & H1: NEURAL ACTIVITY ACROSS CONDITIONS\n")
cat("================================================================================\n\n")

model_rq1_lrt <- lmer(
  beta_value_cz ~ condition * network + hemisphere + n_voxels_c +
    sex + handedness +
    (1 | participant_id) +
    (1 | roi),
  data = data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

print(summary(model_rq1_lrt))
cat("Singularity (ML):", isSingular(model_rq1_lrt), "\n\n")

model_rq1 <- refit_REML_with_p(model_rq1_lrt)
cat("Singularity (REML):", isSingular(model_rq1), "\n")
r2_rq1 <- r2(model_rq1)
cat("R² - Marginal:", round(r2_rq1$R2_marginal, 4),
    "| Conditional:", round(r2_rq1$R2_conditional, 4), "\n\n")

icc_rq1 <- icc_crossed_random_intercepts(model_rq1)

cat("Supplemental Type-III ANOVA for RQ1 (sum-coded refit)\n")
rq1_type3 <- anova_type3_sum_df(
  beta_value_cz ~ condition * network + hemisphere + n_voxels_c + sex + handedness +
    (1 | participant_id) + (1 | roi),
  data
)
write.csv(rq1_type3, paste0(output_dir, "rq1_type3_anova.csv"), row.names = FALSE)

# append to text log:
capture.output({
  cat("\n--- RQ1 supplemental Type-III ANOVA ---\n")
  print(rq1_type3)
}, file = summ_file, append = TRUE)

# --- RQ1: PM-only and DM-only models ---

# Subsets
data_pm_rq1 <- droplevels(dplyr::filter(data, network == "PM"))
data_dm_rq1 <- droplevels(dplyr::filter(data, network == "DM"))

# PM-only
model_rq1_pm_lrt <- lmer(
  beta_value_cz ~ condition + n_voxels_c + sex + handedness +
    (1 | participant_id) + (1 | roi),
  data = data_pm_rq1,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
model_rq1_pm <- refit_REML_with_p(model_rq1_pm_lrt, data = data_pm_rq1)
r2_rq1_pm <- r2(model_rq1_pm)

# DM-only
model_rq1_dm_lrt <- lmer(
  beta_value_cz ~ condition + hemisphere + n_voxels_c + sex + handedness +
    (1 | participant_id) + (1 | roi),
  data = data_dm_rq1,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
model_rq1_dm <- refit_REML_with_p(model_rq1_dm_lrt, data = data_dm_rq1)
r2_rq1_dm <- r2(model_rq1_dm)


# EMMs & planned contrasts for PM/DM (No Error = control)
emm_rq1_pm <- emmeans(model_rq1_pm, ~ condition)
emm_rq1_dm <- emmeans(model_rq1_dm, ~ condition)

con_rq1_pm <- contrast(emm_rq1_pm, method = "trt.vs.ctrl", ref = "No Error", adjust = "bonferroni")
con_rq1_dm <- contrast(emm_rq1_dm, method = "trt.vs.ctrl", ref = "No Error", adjust = "bonferroni")

write.csv(as.data.frame(con_rq1_pm), paste0(output_dir, "rq1_pm_condition_contrasts.csv"), row.names = FALSE)
write.csv(as.data.frame(con_rq1_dm), paste0(output_dir, "rq1_dm_condition_contrasts.csv"), row.names = FALSE)

# Append summaries to log
append_summary(model_rq1_pm, "RQ1 PM-only (REML) — summary")
append_summary(model_rq1_dm, "RQ1 DM-only (REML) — summary")

cat("================================================================================\n")
cat("RQ2 & H2: INDIVIDUAL DIFFERENCES IN CONDITION EFFECTS\n")
cat("================================================================================\n\n")

# PART A: Overall individual differences (pooled across networks)
cat("PART A: Testing overall individual differences in condition effects\n")
cat("(Pooled across both PM and DM networks)\n\n")

model2_baseline_lrt <- lmer(
  beta_value_cz ~ condition * network + hemisphere + n_voxels_c +
    sex + handedness +
    (1 | participant_id) +
    (1 | roi),
  data = data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# Use uncorrelated random slopes (double-bar) for stability
model2_slopes_lrt <- lmer(
  beta_value_cz ~ condition * network + hemisphere + n_voxels_c +
    sex + handedness +
    (1 + condition | participant_id) +
    (1 | roi),
  data = data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

comparison_rq2 <- anova(model2_baseline_lrt, model2_slopes_lrt)
cat("Likelihood ratio test for random slopes:\n")
print(comparison_rq2)
cat("\n")

if (comparison_rq2$`Pr(>Chisq)`[2] < 0.05) {
  model2 <- refit_REML_with_p(model2_slopes_lrt)
  cat("Using RQ2 slopes model (uncorrelated) after LRT.\n\n")
} else {
  model2 <- refit_REML_with_p(model2_baseline_lrt)
  cat("Using RQ2 intercept-only model (no slope improvement).\n\n")
}

cat("Singularity (REML, chosen RQ2 model):", isSingular(model2), "\n")
r2_rq2 <- r2(model2)
cat("R² - Marginal:", round(r2_rq2$R2_marginal, 4),
    "| Conditional:", round(r2_rq2$R2_conditional, 4), "\n\n")

############## Network-specific individual differences

# ---- PM Network ----
cat("--- PM NETWORK ---\n")
cat("Fitting PM-only models with and without random slopes\n\n")

data_pm_rq2 <- droplevels(filter(data, network == "PM"))
cat("PM data: N =", nrow(data_pm_rq2), "observations\n\n")

# Baseline: random intercepts only
model2_pm_baseline <- lmer(
  beta_value_cz ~ condition + n_voxels_c + sex + handedness +
    (1 | participant_id) + (1 | roi),
  data = data_pm_rq2,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# With random slopes for condition
model2_pm_slopes <- lmerTest::lmer(
  beta_value_cz ~ condition + n_voxels_c + sex + handedness +
    (1 + condition | participant_id) + (1 | roi),
  data = data_pm_rq2,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(model2_pm_slopes)

comparison_pm <- anova(model2_pm_baseline, model2_pm_slopes)
cat("PM random slopes likelihood ratio test:\n")
print(comparison_pm)
cat("Singularity (PM slopes):", isSingular(model2_pm_slopes), "\n\n")

if (comparison_pm$`Pr(>Chisq)`[2] < 0.05) {
  model2_pm_final <- refit_REML_with_p(model2_pm_slopes, data = data_pm_rq2)
  cat("PM: Random slopes significantly improve fit (ΔAIC =",
      round(comparison_pm$AIC[1] - comparison_pm$AIC[2], 1), ")\n\n")
} else {
  model2_pm_final <- refit_REML_with_p(model2_pm_baseline, data = data_pm_rq2)
  cat("PM: Random slopes do not improve fit; using intercept-only model\n\n")
}

# ---- DM Network ----
cat("--- DM NETWORK ---\n")
cat("Fitting DM-only models with and without random slopes\n\n")

data_dm_rq2 <- droplevels(filter(data, network == "DM"))
cat("DM data: N =", nrow(data_dm_rq2), "observations\n\n")

# Baseline: random intercepts only
model2_dm_baseline <- lmer(
  beta_value_cz ~ condition + hemisphere + n_voxels_c + sex + handedness +
    (1 | participant_id) + (1 | roi),
  data = data_dm_rq2,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# With random slopes for condition
model2_dm_slopes <- lmer(
  beta_value_cz ~ condition + hemisphere + n_voxels_c + sex + handedness +
    (1 + condition | participant_id) + (1 | roi),
  data = data_dm_rq2,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(model2_dm_slopes)

comparison_dm <- anova(model2_dm_baseline, model2_dm_slopes)
cat("DM random slopes likelihood ratio test:\n")
print(comparison_dm)
cat("Singularity (DM slopes):", isSingular(model2_dm_slopes), "\n\n")

if (comparison_dm$`Pr(>Chisq)`[2] < 0.05) {
  model2_dm_final <- refit_REML_with_p(model2_dm_slopes, data = data_dm_rq2)
  cat("DM: Random slopes significantly improve fit (ΔAIC =",
      round(comparison_dm$AIC[1] - comparison_dm$AIC[2], 1), ")\n\n")
} else {
  model2_dm_final <- refit_REML_with_p(model2_dm_baseline, data = data_dm_rq2)
  cat("DM: Random slopes do not improve fit; using intercept-only model\n\n")
}

# Model fit comparison (R²)
r2_rq2_pm <- r2(model2_pm_final)
r2_rq2_dm <- r2(model2_dm_final)
cat("Model Fit (R²):\n")
cat("  PM - Marginal:", round(r2_rq2_pm$R2_marginal, 4),
    "| Conditional:", round(r2_rq2_pm$R2_conditional, 4), "\n")
cat("  DM - Marginal:", round(r2_rq2_dm$R2_marginal, 4),
    "| Conditional:", round(r2_rq2_dm$R2_conditional, 4), "\n\n")

# AIC comparison (note: different datasets, so only within-network use)
cat("Model Comparison (AIC - lower is better):\n")
cat("  PM AIC:", round(AIC(model2_pm_final), 1), "\n")
cat("  DM AIC:", round(AIC(model2_dm_final), 1), "\n")
cat("  (Note: AICs are not directly comparable across PM vs DM because they fit different subsets.)\n\n")




model2_baseline_lrt <- lmer(
  beta_value_cz ~ condition * network + hemisphere + n_voxels_c +
    sex + handedness +
    (1 | participant_id) +
    (1 | roi),
  data = data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

model2_slopes_lrt <- lmer(
  beta_value_cz ~ condition * network + hemisphere + n_voxels_c +
    sex + handedness +
    (1 + condition | participant_id) +
    (1 | roi),
  data = data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

comparison_rq2 <- anova(model2_baseline_lrt, model2_slopes_lrt)
print(comparison_rq2)
cat("\n")

if (comparison_rq2$`Pr(>Chisq)`[2] < 0.05) {
  model2 <- refit_REML_with_p(model2_slopes_lrt)
  cat("Using RQ2 slopes model (uncorrelated) after LRT.\n\n")
} else {
  model2 <- refit_REML_with_p(model2_baseline_lrt)
  cat("Using RQ2 intercept-only model (no slope improvement).\n\n")
}

cat("Singularity (REML, chosen RQ2 model):", isSingular(model2), "\n")
r2_rq2 <- r2(model2)
cat("R² - Marginal:", round(r2_rq2$R2_marginal, 4),
    "| Conditional:", round(r2_rq2$R2_conditional, 4), "\n\n")

r2_rq2 <- r2(model2)
cat("R² - Marginal:", round(r2_rq2$R2_marginal, 4),
    "| Conditional:", round(r2_rq2$R2_conditional, 4), "\n\n")

r2_rq2 <- r2(model2)
cat("R² - Marginal:", round(r2_rq2$R2_marginal, 4),
    "| Conditional:", round(r2_rq2$R2_conditional, 4), "\n\n")


cat("================================================================================\n")
cat("TYPE III ANOVA FOR NETWORK-SPECIFIC MODELS\n")
cat("================================================================================\n\n")

# ---- PM Network Type III ANOVA ----
cat("--- PM Network: Omnibus Test of Condition Effect ---\n")
rq2_pm_type3 <- anova_type3_sum_df(
  beta_value_cz ~ condition + n_voxels_c + sex + handedness +
    (1 + condition | participant_id) + (1 | roi),
  data_pm_rq2
)
print(rq2_pm_type3)
write.csv(rq2_pm_type3,
          paste0(output_dir, "rq1_pm_type3_anova.csv"),
          row.names = FALSE)
cat("Exported: rq1_pm_type3_anova.csv\n\n")

# ---- DM Network Type III ANOVA ----
cat("--- DM Network: Omnibus Test of Condition Effect ---\n")
rq2_dm_type3 <- anova_type3_sum_df(
  beta_value_cz ~ condition + hemisphere + n_voxels_c + sex + handedness +
    (1 + condition | participant_id) + (1 | roi),
  data_dm_rq2
)
print(rq2_dm_type3)
write.csv(rq2_dm_type3,
          paste0(output_dir, "rq1_dm_type3_anova.csv"),
          row.names = FALSE)
cat("Exported: rq2_dm_type3_anova.csv\n\n")

# Append to summary file
capture.output({
  cat("\n--- RQ2 PM supplemental Type-III ANOVA ---\n")
  print(rq2_pm_type3)
  cat("\n--- RQ2 DM supplemental Type-III ANOVA ---\n")
  print(rq2_dm_type3)
}, file = summ_file, append = TRUE)

cat("--- Interpretation of Type III ANOVAs ---\n")
pm_cond_p <- rq2_pm_type3$`Pr(>F)`[rq2_pm_type3$Effect == "condition"]
dm_cond_p <- rq2_dm_type3$`Pr(>F)`[rq2_dm_type3$Effect == "condition"]

cat("PM Network - Condition effect: F =", 
    round(rq2_pm_type3$`F value`[rq2_pm_type3$Effect == "condition"], 2),
    ", p =", round(pm_cond_p, 4), 
    ifelse(pm_cond_p < .05, "(SIGNIFICANT)", "(not significant)"), "\n")
cat("DM Network - Condition effect: F =", 
    round(rq2_dm_type3$`F value`[rq2_dm_type3$Effect == "condition"], 2),
    ", p =", round(dm_cond_p, 4),
    ifelse(dm_cond_p < .05, "(SIGNIFICANT)", "(not significant)"), "\n\n")



cat("================================================================================\n")
cat("RQ3 & H3: BRAIN-BEHAVIOR RELATIONSHIPS\n")
cat("SEPARATE models for PM and DM to preserve network-specific interactions\n")
cat("================================================================================\n\n")

cat("Accuracy distribution:\n")
print(quantile(data$accuracy, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE))

# RQ3: PM NETWORK
cat("RQ3 PM: Does PM neural activity predict accuracy?\n\n")

data_pm <- droplevels(filter(data, network == "PM"))   

model3_pm_lrt <- lmer(
  accuracy_c ~ beta_value_cz * condition + hemisphere + n_voxels_c +
    sex + handedness +
    (1 | participant_id),
  data = data_pm,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

print(summary(model3_pm_lrt))

cat("Type-III ANOVA (supplemental; local sum-to-zero contrasts):\n")
local({
  op <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(op), add = TRUE)
  print(anova(model3_pm_lrt, type = "III"))
})

cat("Singularity (ML):", isSingular(model3_pm_lrt), "\n\n")

model3_pm <- refit_REML_with_p(model3_pm_lrt)
cat("Singularity (REML):", isSingular(model3_pm), "\n")
r2_rq3_pm <- r2(model3_pm)
cat("R² - Marginal:", round(r2_rq3_pm$R2_marginal, 4),
    "| Conditional:", ifelse(is.na(r2_rq3_pm$R2_conditional), "NA", round(r2_rq3_pm$R2_conditional, 4)), "\n\n")

# RQ3: DM NETWORK
cat("RQ3 DM: Does DM neural activity predict accuracy?\n\n")

data_dm <- droplevels(filter(data, network == "DM"))  

model3_dm_lrt <- lmer(
  accuracy_c ~ beta_value_cz * condition + hemisphere + n_voxels_c +
    sex + handedness +
    (1 | participant_id),
  data = data_dm,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

print(summary(model3_dm_lrt))

cat("Singularity (ML):", isSingular(model3_dm_lrt), "\n\n")

model3_dm <- refit_REML_with_p(model3_dm_lrt)
cat("Singularity (REML):", isSingular(model3_dm), "\n")
r2_rq3_dm <- r2(model3_dm)
cat("R² - Marginal:", round(r2_rq3_dm$R2_marginal, 4),
    "| Conditional:", ifelse(is.na(r2_rq3_dm$R2_conditional), "NA",
                             round(r2_rq3_dm$R2_conditional, 4)), "\n\n")

icc_rq1 <- icc_crossed_random_intercepts(model_rq1)

cat("Supplemental Type-III ANOVA for RQ3 PM (sum-coded refit)\n")


rq3_pm_type3 <- anova_type3_sum_df(
  accuracy_c ~ beta_value_cz * condition + hemisphere + n_voxels_c + sex + handedness +
    (1 + beta_value_cz | participant_id),
  data = data_pm, 
  fixed_factors = c("condition","hemisphere","sex","handedness") 
)

write.csv(rq3_pm_type3, paste0(output_dir, "rq3_pm_type3_anova.csv"), row.names = FALSE)

capture.output({
  cat("\n--- RQ3 PM supplemental Type-III ANOVA ---\n")
  print(rq3_pm_type3)
}, file = summ_file, append = TRUE)

cat("Supplemental Type-III ANOVA for RQ3 DM (sum-coded refit)\n")
rq3_dm_type3 <- anova_type3_sum_df(
  accuracy_c ~ beta_value_cz * condition + hemisphere + n_voxels_c + sex + handedness +
    (1 + beta_value_cz | participant_id),
  data_dm
)
write.csv(rq3_dm_type3, paste0(output_dir, "rq3_dm_type3_anova.csv"), row.names = FALSE)

capture.output({
  cat("\n--- RQ3 DM supplemental Type-III ANOVA ---\n")
  print(rq3_dm_type3)
}, file = summ_file, append = TRUE)


model_rq3_pooled_lrt <- lmer(
  accuracy_c ~ beta_value_cz * condition * network +
    hemisphere + n_voxels_c + sex + handedness +
    (1 | participant_id),
  data  = data,         # full dataset with both networks
  REML  = FALSE,
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl   = list(maxfun = 2e5))
)

summary(model_rq3_pooled_lrt)
isSingular(model_rq3_pooled_lrt)

rq3_pooled_type3 <- anova_type3_sum_df(
  accuracy_c ~ beta_value_cz * condition +  beta_value_cz * network +
    hemisphere + n_voxels_c + sex + handedness +
    (1 + beta_value_cz | participant_id),
  data,
  fixed_factors = c("condition","network","hemisphere","sex","handedness")
)

print(rq3_pooled_type3)
write.csv(rq3_pooled_type3,
          file = paste0(output_dir, "rq3_pooled_type3_anova.csv"),
          row.names = FALSE)


cat("================================================================================\n")
cat("SIMPLE SLOPES: Neural activity → accuracy per condition (unconstrained)\n")
cat("================================================================================\n\n")

# Use emmeans to get slopes per condition
emtrends_pm_sum <- emtrends(model3_pm_sum, ~ condition, var = "beta_value_cz")
emtrends_dm_sum <- emtrends(model3_dm_sum, ~ condition, var = "beta_value_cz")

cat("PM Network - Simple slopes (neural activity → accuracy) by condition:\n")
print(emtrends_pm_sum)
cat("\n")

cat("DM Network - Simple slopes (neural activity → accuracy) by condition:\n")
print(emtrends_dm_sum)
cat("\n")

# Convert emmeans objects to data frames (S4 → data frame)
emtrends_pm_sum_df <- as.data.frame(emtrends_pm_sum)
emtrends_dm_sum_df <- as.data.frame(emtrends_dm_sum)

# Build table with correct column names
simple_slopes_sum <- rbind(
  data.frame(
    network = "PM",
    condition = emtrends_pm_sum_df$condition,
    slope = emtrends_pm_sum_df$beta_value_cz.trend,
    SE = emtrends_pm_sum_df$SE,
    lower_CI = emtrends_pm_sum_df$lower.CL,
    upper_CI = emtrends_pm_sum_df$upper.CL
  ),
  data.frame(
    network = "DM",
    condition = emtrends_dm_sum_df$condition,
    slope = emtrends_dm_sum_df$beta_value_cz.trend,
    SE = emtrends_dm_sum_df$SE,
    lower_CI = emtrends_dm_sum_df$lower.CL,
    upper_CI = emtrends_dm_sum_df$upper.CL
  )
)

cat("\n\nSimple Slopes Summary (Sum-to-Zero Contrasts):\n")
print(simple_slopes_sum)

write.csv(simple_slopes_sum, 
          paste0(output_dir, "rq3_simple_slopes_by_condition_sum_coded.csv"), 
          row.names = FALSE)


cat("Exported: rq3_simple_slopes_by_condition_sum_coded.csv\n\n")


cat("================================================================================\n")
cat("RQ4 & H4: INDIVIDUAL DIFFERENCES IN BRAIN-BEHAVIOR\n")
cat("================================================================================\n\n")


model4_pm_baseline_lrt <- lmer(
  accuracy_c ~ beta_value_cz * condition + hemisphere + n_voxels_c +
    sex + handedness +
    (1 | participant_id),
  data = data_pm,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

model4_pm_slopes_lrt <- lmer(
  accuracy_c ~ beta_value_cz * condition + hemisphere + n_voxels_c +
    sex + handedness +
    (1 + beta_value_cz | participant_id),
  data = data_pm,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(model4_pm_slopes_lrt)

comparison_rq4_pm <- anova(model4_pm_baseline_lrt, model4_pm_slopes_lrt)
print(comparison_rq4_pm)
cat("Singularity (PM slopes, ML):", isSingular(model4_pm_slopes_lrt), "\n\n")

if (comparison_rq4_pm$`Pr(>Chisq)`[2] < 0.05) {
  model4_pm <- refit_REML_with_p(model4_pm_slopes_lrt)
  cat("Random slopes improve PM fit → Using slopes model\n\n")
} else {
  model4_pm <- refit_REML_with_p(model4_pm_baseline_lrt)
  cat("Random slopes do NOT improve PM fit → Using intercept-only model\n\n")
}

cat("Singularity (REML, chosen RQ4 PM):", isSingular(model4_pm), "\n")
r2_rq4_pm <- r2(model4_pm)
cat("R² - Marginal:", round(r2_rq4_pm$R2_marginal, 4),
    "| Conditional:", ifelse(is.na(r2_rq4_pm$R2_conditional), "NA", round(r2_rq4_pm$R2_conditional, 4)), "\n\n")

# RQ4: DM NETWORK
cat("RQ4 DM: Does DM brain-behavior relationship vary across individuals?\n\n")

model4_dm_baseline_lrt <- lmer(
  accuracy_c ~ beta_value_cz * condition + hemisphere + n_voxels_c +
    sex + handedness +
    (1 | participant_id),
  data = data_dm,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

model4_dm_slopes_lrt <- lmer(
  accuracy_c ~ beta_value_cz * condition + hemisphere + n_voxels_c +
    sex + handedness +
    (1 + beta_value_cz | participant_id),
  data = data_dm,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(model4_pm_slopes_lrt)

comparison_rq4_dm <- anova(model4_dm_baseline_lrt, model4_dm_slopes_lrt)
print(comparison_rq4_dm)
cat("Singularity (DM slopes, ML):", isSingular(model4_dm_slopes_lrt), "\n\n")

if (comparison_rq4_dm$`Pr(>Chisq)`[2] < 0.05) {
  model4_dm <- refit_REML_with_p(model4_dm_slopes_lrt)
  cat("Random slopes improve DM fit → Using slopes model\n\n")
} else {
  model4_dm <- refit_REML_with_p(model4_dm_baseline_lrt)
  cat("Random slopes do NOT improve DM fit → Using intercept-only model\n\n")
}

cat("Singularity (REML, chosen RQ4 DM):", isSingular(model4_dm), "\n")
r2_rq4_dm <- r2(model4_dm)
cat("R² - Marginal:", round(r2_rq4_dm$R2_marginal, 4),
    "| Conditional:", ifelse(is.na(r2_rq4_dm$R2_conditional), "NA", round(r2_rq4_dm$R2_conditional, 4)), "\n\n")

sink()


################################################################################
# SENSITIVITY ANALYSES
################################################################################

cat("================================================================================\n")
cat("SENSITIVITY ANALYSES\n")
cat("================================================================================\n\n")

# Sensitivity 1: Right hemisphere only
cat("Sensitivity 1: Right hemisphere only\n")

data_right <- droplevels(data %>% filter(hemisphere == "right"))  # ### CHANGED

sens1 <- lmer(
  beta_value_cz ~ condition * network + n_voxels_c +
    sex + handedness +
    (1 | participant_id) + (1 | roi),
  data = data_right,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

print(summary(sens1))
cat("\n")

# Sensitivity 2: Top ROIs only
cat("Sensitivity 2: Most common ROIs only\n")
roi_counts <- data %>% count(roi) %>% arrange(desc(n))
top_rois <- roi_counts$roi[1:5]
cat("Top 5 ROIs:", paste(top_rois, collapse = ", "), "\n\n")

sens2_dat <- droplevels(filter(data, roi %in% top_rois))  # ### CHANGED
sens2 <- lmer(
  beta_value_cz ~ condition * network + hemisphere + n_voxels_c +
    sex + handedness +
    (1 | participant_id) + (1 | roi),
  data = sens2_dat,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

print(summary(sens2))
cat("\n")

################################################################################
# GENERATE FIGURES
################################################################################

cat("================================================================================\n")
cat("GENERATING FIGURES\n")
cat("================================================================================\n\n")

# Figure 1: Main Results (RQ1) using per-network EMMs
png(paste0(output_dir, "main_results_plot.png"), width = 1200, height = 700)

emm_df_pm <- as.data.frame(emm_rq1_pm) %>% dplyr::mutate(network = "PM")
emm_df_dm <- as.data.frame(emm_rq1_dm) %>% dplyr::mutate(network = "DM")
emm_df    <- dplyr::bind_rows(emm_df_pm, emm_df_dm)

p1 <- ggplot(emm_df, aes(x = condition, y = emmean, fill = network)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8, color = "black", linewidth = 1) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.2, position = position_dodge(0.9), linewidth = 1) +
  scale_fill_manual(values = c("PM" = "#80B1D3", "DM" = "#FB8072")) +
  labs(x = "Grammatical Condition",
       y = "Neural Activity (β, z-scored)",
       title = "RQ1: Condition Effects by Memory Network (per-network models)",
       subtitle = "Estimated marginal means with 95% CIs; pooled interaction model retained separately") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray60"),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    panel.grid.major.y = element_line(color = "gray90")
  )

print(p1)
dev.off()
cat("main_results_plot.png\n\n")


# Figure 2: Brain-Behavior (separate by network)

png(paste0(output_dir, "brain_behavior_plot.png"), width = 1200, height = 700)

p2 <- ggplot(data, aes(x = beta_value, y = accuracy, color = condition, shape = condition)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, aes(fill = condition)) +
  facet_wrap(~ network, scales = "free") +
  scale_color_manual(values = c("No Error" = "#1B9E77", "Easy Error" = "#D95F02", "Hard Error" = "#7570B3")) +
  scale_fill_manual(values = c("No Error" = "#1B9E77", "Easy Error" = "#D95F02", "Hard Error" = "#7570B3")) +
  labs(x = "Neural Activity (β)",
       y = "Accuracy (0–1 proportion)",
       title = "RQ3: Brain–Behavior Relationships: Network-Specific Effects",
       subtitle = "Does neural activity predict accuracy differently in PM vs DM? (visualization)") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray60"),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90")
  )

print(p2)
dev.off()
cat("brain_behavior_plot.png\n\n")

################################################################################
# EXPORT TABLES (Coefficients, Variance Components, R², Planned Contrasts)
################################################################################

cat("================================================================================\n")
cat("EXPORTING TABLES\n")
cat("================================================================================\n\n")

# Table 1: Coefficients - FIXED VERSION (no df column)
cat("Exporting model coefficients\n")

extract_coef <- function(model, model_name) {
  cf <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  # broom.mixed will include p.value if model is lmerModLmerTest
  cf %>%
    dplyr::select(dplyr::any_of(c("term","estimate","std.error","conf.low","conf.high","statistic","df","p.value"))) %>%
    dplyr::mutate(model = model_name) %>%
    dplyr::select(model, dplyr::everything())
}

all_coef <- bind_rows(
  extract_coef(model_rq1,    "RQ1 (Pooled): Neural Activity"),
  extract_coef(model_rq1_pm, "RQ1 (PM-only): Neural Activity"),
  extract_coef(model_rq1_dm, "RQ1 (DM-only): Neural Activity"),
  extract_coef(model2_pm_final, "RQ2: PM Final Slopes Model"),
  extract_coef(model2_dm_final, "RQ2: DM Final Slopes Model"),
  extract_coef(model3_pm,    "RQ3: Brain–Behavior (PM)"),
  extract_coef(model3_dm,    "RQ3: Brain–Behavior (DM)"),
  extract_coef(model4_pm,    "RQ4: Individual Diff (PM)"),
  extract_coef(model4_dm,    "RQ4: Individual Diff (DM)"),
  extract_coef(sens1,        "Sensitivity: Right Hemisphere"),
  extract_coef(sens2,        "Sensitivity: Top ROIs")
)

write.csv(all_coef, paste0(output_dir, "model_coefficients.csv"), row.names = FALSE)
cat("model_coefficients.csv\n\n")

# Table 1b: Planned contrasts (treatment-vs-control) ### CHANGED: new block
cat("Exporting planned contrasts (Easy–No, Hard–No)\n")

# RQ1: contrasts within each network
con_rq1_pm <- contrast(emm_rq1_pm, method = "trt.vs.ctrl", ref = "No Error", adjust = "bonferroni")
con_rq1_dm <- contrast(emm_rq1_dm, method = "trt.vs.ctrl", ref = "No Error", adjust = "bonferroni")
write.csv(as.data.frame(con_rq1_pm), paste0(output_dir, "rq1_pm_simple_slope_contrasts.csv"), row.names = FALSE)
write.csv(as.data.frame(con_rq1_dm), paste0(output_dir, "rq1_dm_simple_slope_contrasts.csv"), row.names = FALSE)
cat("rq1_condition_contrasts_by_network.csv\n")

# RQ3: simple slopes of beta_value_cz on accuracy within each condition, by network
emm_tr_pm <- emtrends(model3_pm, ~ condition, var = "beta_value_cz")
emm_tr_dm <- emtrends(model3_dm, ~ condition, var = "beta_value_cz")
con_tr_pm  <- contrast(emm_tr_pm, "trt.vs.ctrl", ref = "No Error", adjust = "bonferroni")
con_tr_dm  <- contrast(emm_tr_dm, "trt.vs.ctrl", ref = "No Error", adjust = "bonferroni")
write.csv(as.data.frame(con_tr_pm), paste0(output_dir, "rq3_pm_simple_slope_contrasts.csv"), row.names = FALSE)
write.csv(as.data.frame(con_tr_dm), paste0(output_dir, "rq3_dm_simple_slope_contrasts.csv"), row.names = FALSE)
cat("rq3_pm_simple_slope_contrasts.csv\n")
cat("rq3_dm_simple_slope_contrasts.csv\n\n")

# Table 2: Variance Components
# helper already defined above: re_variances()

vc_rq2    <- re_variances(model2,    group = "participant_id", slope_prefix = "condition")
vc_rq4_pm <- re_variances(model4_pm, group = "participant_id", slope_prefix = "beta_value_cz")
vc_rq4_dm <- re_variances(model4_dm, group = "participant_id", slope_prefix = "beta_value_cz")

variance_summary <- data.frame(
  Analysis = c("RQ2: Condition Effects", "RQ4: Brain–Behavior (PM)", "RQ4: Brain–Behavior (DM)"),
  Random_Intercept_Var = c(vc_rq2$ri, vc_rq4_pm$ri, vc_rq4_dm$ri),
  Random_Slope_Var     = c(vc_rq2$rs, vc_rq4_pm$rs, vc_rq4_dm$rs)
)

#write.csv(variance_summary, paste0(output_dir, "variance_components.csv"), row.names = FALSE)
cat("variance_components.csv\n\n")



cat("Exporting effect sizes (R²)\n")

r2_summary <- data.frame(
  Model = c("RQ1 (Pooled): Neural Activity",
            "RQ1 (PM-only): Neural Activity",
            "RQ1 (DM-only): Neural Activity",
            "RQ2 (PM-only): Individual Differences",
            "RQ2 (DM-only): Individual Differences",
            "RQ3: Brain–Behavior (PM)", 
            "RQ3: Brain–Behavior (DM)",
            "RQ4: Brain–Behavior Variation (PM)", 
            "RQ4: Brain–Behavior Variation (DM)"),
  Marginal_R2   = c(round(r2_rq1$R2_marginal, 4),
                    round(r2_rq1_pm$R2_marginal, 4),
                    round(r2_rq1_dm$R2_marginal, 4),
                    round(r2_rq2_pm$R2_marginal, 4),          
                    round(r2_rq2_dm$R2_marginal, 4),          
                    round(r2_rq3_pm$R2_marginal, 4), 
                    round(r2_rq3_dm$R2_marginal, 4),
                    round(r2_rq4_pm$R2_marginal, 4), 
                    round(r2_rq4_dm$R2_marginal, 4)),
  Conditional_R2 = c(round(r2_rq1$R2_conditional, 4),
                     round(r2_rq1_pm$R2_conditional, 4),
                     round(r2_rq1_dm$R2_conditional, 4),
                     round(r2_rq2_pm$R2_conditional, 4),      
                     round(r2_rq2_dm$R2_conditional, 4),      
                     round(r2_rq3_pm$R2_conditional, 4), 
                     round(r2_rq3_dm$R2_conditional, 4),
                     ifelse(is.na(r2_rq4_pm$R2_conditional), NA, round(r2_rq4_pm$R2_conditional, 4)),
                     ifelse(is.na(r2_rq4_dm$R2_conditional), NA, round(r2_rq4_dm$R2_conditional, 4)))
)

write.csv(r2_summary, paste0(output_dir, "effect_sizes_r2.csv"), row.names = FALSE)
cat("effect_sizes_r2.csv\n\n")

# Print for verification
cat("\nEffect Sizes Summary:\n")
print(r2_summary)
append_summary(model_rq1,  "RQ1 final (REML) — summary")
append_summary(model2,     "RQ2 final (REML) — summary")
append_summary(model3_pm,  "RQ3 PM final (REML) — summary & Type-III ANOVA", anova_type3 = TRUE)
append_summary(model3_dm,  "RQ3 DM final (REML) — summary & Type-III ANOVA", anova_type3 = TRUE)
append_summary(model4_pm,  "RQ4 PM final (REML) — summary")
append_summary(model4_dm,  "RQ4 DM final (REML) — summary")

cat("model_summaries.txt\n\n")

################################################################################
# SAVE MODELS
################################################################################

save(
  model_rq1, model2,
  model3_pm, model3_dm,
  model4_pm, model4_dm,
  sens1, sens2,
  file = paste0(output_dir, "all_models_final.RData")
)



cat("================================================================================\n")
cat("BOOTSTRAP CIs FOR FIXED EFFECTS (N = 1000)\n")
cat("================================================================================\n\n")


set.seed(20251114)  # for reproducibility

boot_fixef_ci <- function(mod, model_label, nsim = 1000) {
  # Parametric bootstrap of fixed effects
  bb <- bootMer(
    mod,
    FUN = function(fit) fixef(fit),
    nsim = nsim,
    use.u = TRUE,         # keep random effects structure
    type  = "parametric",
    parallel = "no"       # change to "multicore" or "snow" if you want later
  )
  
  orig     <- fixef(mod)
  boot_mat <- bb$t        # nsim x p matrix of bootstrapped coefficients
  
  # Percentile 95% CIs
  ci_mat <- t(apply(
    boot_mat, 2,
    quantile, probs = c(0.025, 0.975), na.rm = TRUE
  ))
  
  out <- data.frame(
    model     = model_label,
    term      = names(orig),
    estimate  = as.numeric(orig),
    boot_mean = apply(boot_mat, 2, mean, na.rm = TRUE),
    boot_se   = apply(boot_mat, 2, sd,   na.rm = TRUE),
    ci_lower  = ci_mat[, 1],
    ci_upper  = ci_mat[, 2],
    row.names = NULL
  )
  
  out
}

boot_results <- dplyr::bind_rows(
  boot_fixef_ci(model2_pm_slopes,      "model2_pm_slopes"),
  boot_fixef_ci(model2_dm_slopes,      "model2_dm_slopes"),
  boot_fixef_ci(model4_pm_slopes_lrt,  "model4_pm_slopes_lrt"),
  boot_fixef_ci(model4_dm_slopes_lrt,  "model4_dm_slopes_lrt")
)

print(boot_results)

boot_file <- paste0(output_dir, "bootstrap_fixed_effects_1000.csv")
write.csv(boot_results, boot_file, row.names = FALSE)
cat("Exported bootstrap CIs to:", boot_file, "\n\n")



cat("================================================================================\n")
cat("BOOTSTRAP CIs FOR RANDOM EFFECTS (N = 1000)\n")
cat("================================================================================\n\n")


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

### USING THE FUNCTIONS

boot_re_sd_results <- dplyr::bind_rows(
  boot_ranef_var_ci(model2_pm_slopes,     "model2_pm_slopes",     nsim = 1000, return_sd = TRUE),
  boot_ranef_var_ci(model2_dm_slopes,     "model2_dm_slopes",     nsim = 1000, return_sd = TRUE),
  boot_ranef_var_ci(model4_pm_slopes_lrt, "model4_pm_slopes_lrt", nsim = 1000, return_sd = TRUE),
  boot_ranef_var_ci(model4_dm_slopes_lrt, "model4_dm_slopes_lrt", nsim = 1000, return_sd = TRUE)
)

boot_re_sd_file <- paste0(output_dir, "bootstrap_random_effects_sd_1000.csv")
write.csv(boot_re_sd_results, boot_re_sd_file, row.names = FALSE)