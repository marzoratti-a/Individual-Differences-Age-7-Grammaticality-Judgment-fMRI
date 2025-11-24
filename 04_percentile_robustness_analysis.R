################################################################################
# DMPM Paper - Robustness (RS-only) with per-config tabs for RQ1 and RQ3
# Tabs created for each config:
#   - Top25 (Base)
#   - Top25 + CELF (adds language_score if present)
#   - Top25 + KBIT (adds IQ_score if present)
#   - Top50
#   - Top75
#   - Top100
#
# RQ1 (β outcome): PM no hemisphere; DM with hemisphere; RS on condition
# RQ3 (ACC outcome): hemisphere kept; RS on beta
# All fits reported with REML=TRUE 
################################################################################

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(performance)
  library(broom.mixed)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(openxlsx)
  library(rlang)
})

## ---------------------------------------------------------------------------
## Paths
## ---------------------------------------------------------------------------
results_dir <- "C:/Users/anama/OneDrive/Documents/Diss1_DMPM_Paper/2025_Results/"
merged_file <- file.path(results_dir, "merged_data_for_analysis.csv")

beta_files <- list(
  Top25  = file.path(results_dir, "Robustness_Analysis/top25percentile_beta_values.csv"),
  Top50  = file.path(results_dir, "Robustness_Analysis/top50percentile_beta_values.csv"),
  Top75  = file.path(results_dir, "Robustness_Analysis/top75percentile_beta_values.csv"),
  Top100 = file.path(results_dir, "Robustness_Analysis/top100percentile_beta_values.csv")
)

wb_path <- file.path(results_dir, "DMPM_Robustness_111625.xlsx")

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------
recode_condition_labels <- function(x) {
  dplyr::recode(
    x,
    "correct_grammar"  = "No Error",
    "plurality_error"  = "Easy Error",
    "finiteness_error" = "Hard Error",
    "No Error"   = "No Error",
    "Easy Error" = "Easy Error",
    "Hard Error" = "Hard Error",
    .default = x
  )
}

standardize_beta_df <- function(df_raw) {
  df <- df_raw
  if (!"participant_id" %in% names(df)) {
    alt <- intersect(c("participant","sub","Subject","id"), names(df))
    if (length(alt)) df <- dplyr::rename(df, participant_id = !!rlang::sym(alt[1]))
  }
  if (!"n_voxels_used" %in% names(df)) {
    alt <- intersect(c("n_voxels_u","n_voxels","nVoxels"), names(df))
    if (length(alt)) df <- dplyr::rename(df, n_voxels_used = !!rlang::sym(alt[1]))
  }
  if (!"beta_value" %in% names(df)) {
    alt <- intersect(c("beta","beta_val","beta_z","betaValue"), names(df))
    if (length(alt)) df <- dplyr::rename(df, beta_value = !!rlang::sym(alt[1]))
  }
  needed <- c("participant_id","condition","roi","beta_value","n_voxels_used")
  missing <- setdiff(needed, names(df))
  if (length(missing)) stop("Robustness file is missing columns: ", paste(missing, collapse=", "))
  df %>%
    dplyr::mutate(
      participant_id = as.character(participant_id),
      roi            = as.character(roi),
      condition      = recode_condition_labels(condition)
    )
}

tidy_fixed <- function(mod, lbl) {
  broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE) %>%
    dplyr::select(dplyr::any_of(c("term","estimate","std.error",
                                  "conf.low","conf.high","statistic","df","p.value"))) %>%
    dplyr::mutate(Model = lbl, .before = 1)
}

varcomp_row <- function(mod, lbl, slope_label) {
  vc <- as.data.frame(VarCorr(mod))
  sig2 <- sigma(mod)^2
  v_part <- vc$vcov[vc$grp == "participant_id" & vc$var1 == "(Intercept)"][1]
  v_roi  <- vc$vcov[vc$grp == "roi"            & vc$var1 == "(Intercept)"][1]
  
  slope_components <- vc[vc$grp == "participant_id" & vc$var1 != "(Intercept)", ]
  slope_components <- slope_components %>%
    dplyr::distinct(.data$grp, .data$var1, .keep_all = TRUE)
  
  if (nrow(slope_components) > 0) {
    output <- tibble::tibble(
      Model = lbl,
      Component = slope_components$var1,
      Random_Intercept_Var_participant = v_part,
      Random_Intercept_Var_roi = v_roi,
      Random_Slope_Var = slope_components$vcov,
      Slope_Type = slope_label,
      Residual_Var = sig2,
      ICC = (v_part + v_roi) / (v_part + v_roi + sig2),
      REML_Criterion = -2 * as.numeric(logLik(mod)),  # ADD HERE
      AIC = AIC(mod), 
      BIC = BIC(mod), 
      logLik = as.numeric(logLik(mod)),
      isSingular = isSingular(mod)
    )
  } else {
    output <- tibble::tibble(
      Model = lbl,
      Component = NA_character_,
      Random_Intercept_Var_participant = v_part,
      Random_Intercept_Var_roi = v_roi,
      Random_Slope_Var = NA_real_,
      Slope_Type = slope_label,
      Residual_Var = sig2,
      ICC = (v_part + v_roi) / (v_part + v_roi + sig2),
      REML_Criterion = -2 * as.numeric(logLik(mod)),  # AND ADD HERE TOO
      AIC = AIC(mod), 
      BIC = BIC(mod), 
      logLik = as.numeric(logLik(mod)),
      isSingular = isSingular(mod)
    )
  }
  output
}

append_covars <- function(df, add_covars) {
  if (is.null(add_covars) | !length(add_covars)) return(list(df=df, used_covars=character(0)))
  used <- character(0)
  for (v in add_covars) {
    if (v %in% names(df)) {
      if (is.numeric(df[[v]])) {
        df[[paste0(v, "_c")]] <- scale(df[[v]], center = TRUE, scale = TRUE)[,1]
        used <- c(used, paste0(v, "_c"))
      } else {
        used <- c(used, v)
      }
    }
  }
  list(df=df, used_covars=used)
}

write_block <- function(wb, sheet, coefs, contrasts, varcomp, r2) {
  addWorksheet(wb, sheet)
  start <- 1
  writeData(wb, sheet, paste0("Block: ", sheet), startRow = start); start <- start + 2
  if (nrow(coefs))     { writeData(wb, sheet, "Fixed Effects (REML)", startRow = start); start <- start + 1
  writeData(wb, sheet, coefs, startRow = start); start <- start + nrow(coefs) + 2 }
  if (nrow(contrasts)){ writeData(wb, sheet, "Contrasts", startRow = start); start <- start + 1
  writeData(wb, sheet, contrasts, startRow = start); start <- start + nrow(contrasts) + 2 }
  if (nrow(varcomp))  { writeData(wb, sheet, "Random Effects, ICC, Fit", startRow = start); start <- start + 1
  writeData(wb, sheet, varcomp, startRow = start); start <- start + nrow(varcomp) + 2 }
  if (nrow(r2))       { writeData(wb, sheet, "R2 (performance::r2)", startRow = start); start <- start + 1
  writeData(wb, sheet, r2, startRow = start); start <- start + nrow(r2) + 2 }
}

## ---------------------------------------------------------------------------
## Load + prepare merged data
## ---------------------------------------------------------------------------
if (!file.exists(merged_file)) stop("Merged data not found at: ", merged_file)
data <- read.csv(merged_file)

voxel_minimum <- 20
data <- data %>% filter(n_voxels_used >= voxel_minimum) %>% droplevels()

data <- data %>%
  mutate(
    condition   = factor(condition,
                         levels = c("correct_grammar","plurality_error","finiteness_error"),
                         labels = c("No Error","Easy Error","Hard Error")),
    network     = factor(network,    levels = c("PM","DM")),
    hemisphere  = factor(hemisphere, levels = c("right","left")),
    participant_id = factor(participant_id),
    roi         = factor(roi),
    sex         = factor(sex, levels = c("0","1"), labels = c("Male","Female")),
    handedness  = factor(handedness, levels = c("3","4","5")),
    beta_value_cz = scale(beta_value,    center = TRUE, scale = TRUE)[,1],
    n_voxels_c    = scale(n_voxels_used, center = TRUE, scale = TRUE)[,1],
    accuracy_c    = scale(accuracy,      center = TRUE, scale = TRUE)[,1]
  )

data$condition  <- stats::relevel(data$condition, ref = "No Error")
contrasts(data$condition)  <- contr.treatment(levels(data$condition))
data$network    <- stats::relevel(data$network, ref = "PM")
data$hemisphere <- stats::relevel(data$hemisphere, ref = "right")
data$sex        <- stats::relevel(data$sex, ref = "Male")
data$handedness <- stats::relevel(data$handedness, ref = "3")

behavioral_keys <- c("participant_id","roi","condition")
behavioral_keep <- c(behavioral_keys, "network","hemisphere","accuracy","sex","handedness","n_voxels_used")
behavioral_slice <- data %>%
  dplyr::select(dplyr::any_of(behavioral_keep)) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    participant_id = as.character(participant_id),
    roi            = as.character(roi),
    condition      = as.character(condition)
  )

## ---------------------------------------------------------------------------
## Model fitters (RS-only)
## ---------------------------------------------------------------------------
fit_rq3_accuracy_rs <- function(df, label = "RQ3 ACC RS", add_covars = NULL) {
  ac <- append_covars(df, add_covars); df <- ac$df; used_covars <- ac$used_covars
  covar_terms <- if (length(used_covars)) paste("+", paste(used_covars, collapse = " + ")) else ""
  
  data_pm <- droplevels(filter(df, network == "PM"))
  data_dm <- droplevels(filter(df, network == "DM"))
  
  fm_txt <- paste0(
    "accuracy_c ~ beta_value_cz * condition + hemisphere + n_voxels_c + sex + handedness ",
    covar_terms, " + (1 + beta_value_cz | participant_id)"
  )
  fm <- as.formula(fm_txt)
  
  m_pm <- lmerTest::lmer(fm, data = data_pm, REML = TRUE,
                         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  m_dm <- lmerTest::lmer(fm, data = data_dm, REML = TRUE,
                         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  
  coefs <- bind_rows(
    tidy_fixed(m_pm, paste0(label," | PM")),
    tidy_fixed(m_dm, paste0(label," | DM"))
  )
  
  emm_tr_pm <- emtrends(m_pm, ~ condition, var = "beta_value_cz")
  emm_tr_dm <- emtrends(m_dm, ~ condition, var = "beta_value_cz")
  con_pm <- as.data.frame(contrast(emm_tr_pm, "trt.vs.ctrl", ref = "No Error", adjust = "bonferroni"))
  con_dm <- as.data.frame(contrast(emm_tr_dm, "trt.vs.ctrl", ref = "No Error", adjust = "bonferroni"))
  con_pm$Model <- paste0(label," | PM")
  con_dm$Model <- paste0(label," | DM")
  contrasts_tbl <- bind_rows(con_pm, con_dm) %>% relocate(Model)
  
  varcomp <- bind_rows(
    varcomp_row(m_pm, paste0(label," | PM"), "beta"),
    varcomp_row(m_dm, paste0(label," | DM"), "beta")
  )
  
  r2_tbl <- tibble(
    Model = c(paste0(label," | PM"), paste0(label," | DM")),
    Marginal_R2   = c(performance::r2(m_pm)$R2_marginal, performance::r2(m_dm)$R2_marginal),
    Conditional_R2= c(performance::r2(m_pm)$R2_conditional, performance::r2(m_dm)$R2_conditional)
  )
  
  list(coefs = coefs, contrasts = contrasts_tbl, varcomp = varcomp, r2 = r2_tbl)
}

fit_rq1_beta_rs <- function(df, label = "RQ1 β RS", add_covars = NULL) {
  ac <- append_covars(df, add_covars); df <- ac$df; used_covars <- ac$used_covars
  covar_terms <- if (length(used_covars)) paste("+", paste(used_covars, collapse = " + ")) else ""
  
  data_pm <- droplevels(filter(df, network == "PM"))
  data_dm <- droplevels(filter(df, network == "DM"))
  
  fm_pm_txt <- paste0(
    "beta_value_cz ~ condition + n_voxels_c + sex + handedness ",
    covar_terms, " + (1 + condition | participant_id) + (1 | roi)"
  )
  fm_dm_txt <- paste0(
    "beta_value_cz ~ condition + hemisphere + n_voxels_c + sex + handedness ",
    covar_terms, " + (1 + condition | participant_id) + (1 | roi)"
  )
  
  m_pm <- lmerTest::lmer(as.formula(fm_pm_txt), data = data_pm, REML = TRUE,
                         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  m_dm <- lmerTest::lmer(as.formula(fm_dm_txt), data = data_dm, REML = TRUE,
                         control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  
  coefs <- bind_rows(
    tidy_fixed(m_pm, paste0(label," | PM")),
    tidy_fixed(m_dm, paste0(label," | DM"))
  )
  
  cond_contrasts <- function(mod, mdl_lab) {
    emm <- emmeans(mod, ~ condition)
    out <- as.data.frame(contrast(emm, "trt.vs.ctrl", ref = "No Error", adjust = "bonferroni"))
    out$Model <- mdl_lab
    relocate(out, Model)
  }
  contrasts_tbl <- bind_rows(
    cond_contrasts(m_pm, paste0(label," | PM")),
    cond_contrasts(m_dm, paste0(label," | DM"))
  )
  
  varcomp <- bind_rows(
    varcomp_row(m_pm, paste0(label," | PM"), "condition"),
    varcomp_row(m_dm, paste0(label," | DM"), "condition")
  )
  
  r2_tbl <- tibble(
    Model = c(paste0(label," | PM"), paste0(label," | DM")),
    Marginal_R2   = c(performance::r2(m_pm)$R2_marginal, performance::r2(m_dm)$R2_marginal),
    Conditional_R2= c(performance::r2(m_pm)$R2_conditional, performance::r2(m_dm)$R2_conditional)
  )
  
  list(coefs = coefs, contrasts = contrasts_tbl, varcomp = varcomp, r2 = r2_tbl)
}


## ---------------------------------------------------------------------------
## BOOTSTRAP FUNCTIONS
## ---------------------------------------------------------------------------

boot_fixef_ci <- function(mod, model_label, nsim = 1000, data = NULL) {
  if (is.null(data)) {
    data <- model.frame(mod)
  }
  
  bb <- suppressWarnings(
    bootMer(
      mod,
      FUN = function(fit) fixef(fit),
      nsim = nsim,
      use.u = TRUE,
      type = "parametric",
      parallel = "no"
    )
  )
  
  orig <- fixef(mod)
  boot_mat <- bb$t
  ci_mat <- t(apply(boot_mat, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  
  data.frame(
    model = model_label,
    term = names(orig),
    estimate = as.numeric(orig),
    boot_mean = apply(boot_mat, 2, mean, na.rm = TRUE),
    boot_se = apply(boot_mat, 2, sd, na.rm = TRUE),
    ci_lower = ci_mat[, 1],
    ci_upper = ci_mat[, 2],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

boot_ranef_var_ci <- function(mod, model_label, nsim = 1000, return_sd = TRUE, data = NULL) {
  if (is.null(data)) {
    data <- model.frame(mod)
  }
  
  vc_orig <- as.data.frame(VarCorr(mod))
  vc_vars <- vc_orig[vc_orig$var2 == "" | is.na(vc_orig$var2), ]
  comp_names <- paste0(vc_vars$grp, ":", vc_vars$var1)
  
  get_vc <- function(fit) {
    vc_curr <- as.data.frame(VarCorr(fit))
    vc_curr <- vc_curr[vc_curr$var2 == "" | is.na(vc_curr$var2), ]
    curr_names <- paste0(vc_curr$grp, ":", vc_curr$var1)
    curr_vals <- structure(vc_curr$vcov, names = curr_names)
    out <- curr_vals[comp_names]
    names(out) <- comp_names
    out
  }
  
  orig_vals <- vc_vars$vcov
  if (return_sd) orig_vals <- sqrt(pmax(orig_vals, 0))
  names(orig_vals) <- comp_names
  
  set.seed(20251114)
  bb <- suppressWarnings(
    bootMer(
      mod, 
      FUN = get_vc, 
      nsim = nsim, 
      use.u = TRUE, 
      type = "parametric", 
      parallel = "no"
    )
  )
  
  boot_mat <- bb$t
  if (return_sd) boot_mat <- sqrt(pmax(boot_mat, 0))
  ci_mat <- t(apply(boot_mat, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  
  data.frame(
    model = model_label,
    component = comp_names,
    estimate = as.numeric(orig_vals),
    boot_mean = apply(boot_mat, 2, mean, na.rm = TRUE),
    boot_se = apply(boot_mat, 2, sd, na.rm = TRUE),
    ci_lower = ci_mat[, 1],
    ci_upper = ci_mat[, 2],
    scale = if (return_sd) "SD" else "Variance",
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

## ---------------------------------------------------------------------------
## Config runner with bootstrap
## ---------------------------------------------------------------------------

run_config_with_bootstrap <- function(config_key,
                                      base_label,
                                      beta_path,
                                      add_covars = NULL,
                                      celf_expected = "language_score",
                                      kbit_expected = "IQ_score",
                                      nsim_boot = 1000) {
  
  if (!file.exists(beta_path)) {
    msg <- paste0("Missing ", config_key, " beta file: ", beta_path)
    run_log <<- add_row(run_log, Timestamp = as.character(Sys.time()), Note = msg)
    return(list(boot_fixef = data.frame(), boot_ranef = data.frame()))
  }
  
  # -----------------------------
  # Load / merge / prep data
  # -----------------------------
  bet <- read.csv(beta_path) %>% standardize_beta_df()
  
  merged_rb <- inner_join(
    behavioral_slice, bet,
    by = c("participant_id", "roi", "condition"),
    suffix = c(".behav", ".beta")
  )
  
  if ("n_voxels_used.behav" %in% names(merged_rb) | "n_voxels_used.beta" %in% names(merged_rb)) {
    merged_rb <- merged_rb %>%
      mutate(n_voxels_used = coalesce(.data$`n_voxels_used.behav`, .data$`n_voxels_used.beta`)) %>%
      select(-any_of(c("n_voxels_used.behav", "n_voxels_used.beta")))
  }
  
  merged_rb <- merged_rb %>%
    mutate(
      condition      = factor(condition, levels = c("No Error", "Easy Error", "Hard Error")),
      network        = factor(network,    levels = c("PM", "DM")),
      hemisphere     = factor(hemisphere, levels = c("right", "left")),
      participant_id = factor(participant_id),
      roi            = factor(roi),
      n_voxels_c     = scale(n_voxels_used, center = TRUE, scale = TRUE)[, 1],
      accuracy_c     = scale(accuracy,     center = TRUE, scale = TRUE)[, 1],
      beta_value_cz  = scale(beta_value,   center = TRUE, scale = TRUE)[, 1]
    )
  
  # carry covars from original data if needed
  for (cv in c(celf_expected, kbit_expected)) {
    if (cv %in% names(data) && !(cv %in% names(merged_rb))) {
      merged_rb <- merged_rb %>%
        left_join(distinct(data, participant_id, !!rlang::sym(cv)), by = "participant_id")
    }
  }
  
  # Center covariates once here
  ac <- append_covars(merged_rb, add_covars)
  merged_rb  <- ac$df
  used_covars <- ac$used_covars
  covar_terms <- if (length(used_covars)) paste("+", paste(used_covars, collapse = " + ")) else ""
  
  data_pm <- droplevels(filter(merged_rb, network == "PM"))
  data_dm <- droplevels(filter(merged_rb, network == "DM"))
  
  rq1_label <- paste0("RQ1 | ", base_label)
  rq3_label <- paste0("RQ3 | ", base_label)
  
  # -----------------------------
  # Main RS fits (for Excel sheets)
  # -----------------------------
  cat("\n--- Fitting RQ3 (", base_label, ") ---\n")
  rq3 <- fit_rq3_accuracy_rs(
    merged_rb,
    label = paste0("RQ3 ACC RS (", base_label, ")"),
    add_covars = add_covars
  )
  if (rq3_label %in% names(wb)) removeWorksheet(wb, rq3_label)
  write_block(wb, rq3_label, rq3$coefs, rq3$contrasts, rq3$varcomp, rq3$r2)

  cat("\n--- Fitting RQ1 (", base_label, ") ---\n")
  rq1 <- fit_rq1_beta_rs(
    merged_rb,
    label = paste0("RQ1 β RS (", base_label, ")"),
    add_covars = add_covars
  )
  if (rq1_label %in% names(wb)) removeWorksheet(wb, rq1_label)
  write_block(wb, rq1_label, rq1$coefs, rq1$contrasts, rq1$varcomp, rq1$r2)

  # -----------------------------
  # Models for bootstrap
  # -----------------------------
  fm_rq3 <- as.formula(paste0(
    "accuracy_c ~ beta_value_cz * condition + hemisphere + n_voxels_c + sex + handedness ",
    covar_terms, " + (1 + beta_value_cz | participant_id)"
  ))
  
  fm_pm <- as.formula(paste0(
    "beta_value_cz ~ condition + n_voxels_c + sex + handedness ",
    covar_terms, " + (1 + condition | participant_id) + (1 | roi)"
  ))
  
  fm_dm <- as.formula(paste0(
    "beta_value_cz ~ condition + hemisphere + n_voxels_c + sex + handedness ",
    covar_terms, " + (1 + condition | participant_id) + (1 | roi)"
  ))
  
  cat("\n--- Fitting models for bootstrap (", base_label, ") ---\n")

  boot_fixef_rq3_pm <- data.frame()
  boot_fixef_rq3_dm <- data.frame()
  boot_fixef_rq1_pm <- data.frame()
  boot_fixef_rq1_dm <- data.frame()
  boot_ranef_rq3_pm <- data.frame()
  boot_ranef_rq3_dm <- data.frame()
  boot_ranef_rq1_pm <- data.frame()
  boot_ranef_rq1_dm <- data.frame()

  tryCatch({
    m_rq3_pm <- lmerTest::lmer(fm_rq3, data = data_pm, REML = TRUE,
                               control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
    m_rq3_dm <- lmerTest::lmer(fm_rq3, data = data_dm, REML = TRUE,
                               control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
    
    m_rq1_pm <- lmerTest::lmer(fm_pm, data = data_pm, REML = TRUE,
                               control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
    m_rq1_dm <- lmerTest::lmer(fm_dm, data = data_dm, REML = TRUE,
                               control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

    cat("\n--- Bootstrapping Fixed Effects (", base_label, ") ---\n")
    boot_fixef_rq3_pm <- suppressWarnings(boot_fixef_ci(m_rq3_pm, paste0("RQ3 PM | ", base_label), nsim = nsim_boot))
    boot_fixef_rq3_dm <- suppressWarnings(boot_fixef_ci(m_rq3_dm, paste0("RQ3 DM | ", base_label), nsim = nsim_boot))
    boot_fixef_rq1_pm <- suppressWarnings(boot_fixef_ci(m_rq1_pm, paste0("RQ1 PM | ", base_label), nsim = nsim_boot))
    boot_fixef_rq1_dm <- suppressWarnings(boot_fixef_ci(m_rq1_dm, paste0("RQ1 DM | ", base_label), nsim = nsim_boot))

    cat("\n--- Bootstrapping Random Effects (", base_label, ") ---\n")
    boot_ranef_rq3_pm <- suppressWarnings(boot_ranef_var_ci(m_rq3_pm, paste0("RQ3 PM | ", base_label), nsim = nsim_boot, return_sd = TRUE))
    boot_ranef_rq3_dm <- suppressWarnings(boot_ranef_var_ci(m_rq3_dm, paste0("RQ3 DM | ", base_label), nsim = nsim_boot, return_sd = TRUE))
    boot_ranef_rq1_pm <- suppressWarnings(boot_ranef_var_ci(m_rq1_pm, paste0("RQ1 PM | ", base_label), nsim = nsim_boot, return_sd = TRUE))
    boot_ranef_rq1_dm <- suppressWarnings(boot_ranef_var_ci(m_rq1_dm, paste0("RQ1 DM | ", base_label), nsim = nsim_boot, return_sd = TRUE))

  }, error = function(e) {
    cat("\nWARNING - Error in bootstrap for", base_label, ":\n")
    cat(e$message, "\n")
  })
  
  run_log <<- add_row(run_log, Timestamp = as.character(Sys.time()),
                      Note = paste("Completed config:", base_label))

  list(
    boot_fixef = bind_rows(boot_fixef_rq3_pm, boot_fixef_rq3_dm,
                           boot_fixef_rq1_pm, boot_fixef_rq1_dm),
    boot_ranef = bind_rows(boot_ranef_rq3_pm, boot_ranef_rq3_dm,
                           boot_ranef_rq1_pm, boot_ranef_rq1_dm)
  )
}
  

## ---------------------------------------------------------------------------
## Workbook + Run Log
## ---------------------------------------------------------------------------
wb <- createWorkbook()
addWorksheet(wb, "Run Log")
run_log <- tibble(Timestamp = as.character(Sys.time()),
                  Note = "Started workbook assembly (RS-only, per-config tabs).")
writeData(wb, "Run Log", run_log)

## ---------------------------------------------------------------------------
## EXECUTE CONFIGURATIONS
## ---------------------------------------------------------------------------

cat("================================================================================\n")
cat("ROBUSTNESS ANALYSIS WITH BOOTSTRAP CIs\n")
cat("================================================================================\n\n")

all_boot_fixef <- list()
all_boot_ranef <- list()

# Top25 (Base)
result_top25 <- run_config_with_bootstrap("Top25", "Top25 (Base)", beta_files$Top25, add_covars = NULL, nsim_boot = 1000)
all_boot_fixef$Top25 <- result_top25$boot_fixef
all_boot_ranef$Top25 <- result_top25$boot_ranef

# Top25 + CELF
celf_cov <- if ("language_score" %in% names(data)) "language_score" else NULL
if (is.null(celf_cov)) {
  run_log <- add_row(run_log, Timestamp = as.character(Sys.time()),
                     Note = "Top25 + CELF requested, but 'language_score' not found; running without CELF.")
}
result_celf <- run_config_with_bootstrap("Top25", "Top25 + CELF", beta_files$Top25, add_covars = celf_cov, nsim_boot = 1000)
all_boot_fixef$CELF <- result_celf$boot_fixef
all_boot_ranef$CELF <- result_celf$boot_ranef

# Top25 + KBIT
kbit_cov <- if ("IQ_score" %in% names(data)) "IQ_score" else NULL
if (is.null(kbit_cov)) {
  run_log <- add_row(run_log, Timestamp = as.character(Sys.time()),
                     Note = "Top25 + KBIT requested, but 'IQ_score' not found; running without KBIT.")
}
result_kbit <- run_config_with_bootstrap("Top25", "Top25 + KBIT", beta_files$Top25, add_covars = kbit_cov, nsim_boot = 1000)
all_boot_fixef$KBIT <- result_kbit$boot_fixef
all_boot_ranef$KBIT <- result_kbit$boot_ranef

# Top50
result_top50 <- run_config_with_bootstrap("Top50", "Top50", beta_files$Top50, add_covars = NULL, nsim_boot = 1000)
all_boot_fixef$Top50 <- result_top50$boot_fixef
all_boot_ranef$Top50 <- result_top50$boot_ranef

# Top75
result_top75 <- run_config_with_bootstrap("Top75", "Top75", beta_files$Top75, add_covars = NULL, nsim_boot = 1000)
all_boot_fixef$Top75 <- result_top75$boot_fixef
all_boot_ranef$Top75 <- result_top75$boot_ranef

# Top100
result_top100 <- run_config_with_bootstrap("Top100", "Top100", beta_files$Top100, add_covars = NULL, nsim_boot = 1000)
all_boot_fixef$Top100 <- result_top100$boot_fixef
all_boot_ranef$Top100 <- result_top100$boot_ranef

## ---------------------------------------------------------------------------
## Save main workbook
## ---------------------------------------------------------------------------
addWorksheet(wb, "Run Log (final)")
writeData(wb, "Run Log (final)", run_log, startRow = 1, colNames = TRUE)
saveWorkbook(wb, wb_path, overwrite = TRUE)
message("Main analysis outputs written to: ", wb_path)

## ---------------------------------------------------------------------------
## CREATE BOOTSTRAP WORKBOOKS
## ---------------------------------------------------------------------------

# Fixed Effects Workbook
wb_fixef <- createWorkbook()
addWorksheet(wb_fixef, "Summary")
writeData(wb_fixef, "Summary", tibble(
  Note = "Bootstrap Fixed Effects CIs (N=1000) for RQ1 and RQ3 across all robustness configurations",
  Generated = as.character(Sys.time())
))

for (config_name in names(all_boot_fixef)) {
  if (nrow(all_boot_fixef[[config_name]]) > 0) {
    sheet_name <- substr(config_name, 1, 31)
    addWorksheet(wb_fixef, sheet_name)
    writeData(wb_fixef, sheet_name, all_boot_fixef[[config_name]])
  }
}

fixef_path <- file.path(results_dir, "Bootstrap_FixedEffects_RQ1RQ3.xlsx")
saveWorkbook(wb_fixef, fixef_path, overwrite = TRUE)
cat("Fixed Effects Bootstrap workbook saved to:", fixef_path, "\n")

# Random Effects Workbook
wb_ranef <- createWorkbook()
addWorksheet(wb_ranef, "Summary")
writeData(wb_ranef, "Summary", tibble(
  Note = "Bootstrap Random Effects SDs (N=1000) for RQ1 and RQ3 across all robustness configurations",
  Generated = as.character(Sys.time())
))

for (config_name in names(all_boot_ranef)) {
  if (nrow(all_boot_ranef[[config_name]]) > 0) {
    sheet_name <- substr(config_name, 1, 31)
    addWorksheet(wb_ranef, sheet_name)
    writeData(wb_ranef, sheet_name, all_boot_ranef[[config_name]])
  }
}

ranef_path <- file.path(results_dir, "Bootstrap_RandomEffects_RQ1RQ3.xlsx")
saveWorkbook(wb_ranef, ranef_path, overwrite = TRUE)
cat("Random Effects Bootstrap workbook saved to:", ranef_path, "\n\n")
