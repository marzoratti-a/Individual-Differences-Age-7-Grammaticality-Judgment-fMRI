# Loading libraries for analyses
# update.packages(ask = FALSE)
library(lme4) # This is the package that actually fits mixed effects models
library(lmerTest) # This package allows for significance tests
library(performance) # This package allows for computation of the intraclass correlation coefficient
library(summarytools) # This allows for a summary of the dataframe
library(psych) # More summary statistics
library(sjPlot) # For model diagnostics & other useful plots
library(sjmisc)
library(sjlabelled)
library(HLMdiag) # For model diagnostics
library(lmtest) # To condut LRT 
library(devtools)
library(imputeTS)
library(irr)
library(boot)
library(dplyr)
library(officer)
library(flextable)

# libraries for plotting
library(lubridate)
library(ggplot2)
library(scales)
library(gridExtra)
library(ggthemes)
library(gridExtra)
library("Hmisc")

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
  as.numeric(VarCorr(model)$sub[2, 2])
}

# Define a function to extract the random intercept term variance component
extract_random_intercept <- function(model) {
  as.numeric(VarCorr(model)$sub[1, 1])
}

# Importing data for RQ1 (DV= dm/pm, L1= cond, L2= sub)
dat = read.csv("C:/Users/anama/OneDrive/Documents/DMPM_Paper/Gram_ROISep.csv")

# Transforming variables to numeric/continuous
dat$handedness<- as.numeric(as.character(dat$handedness))
dat$kbit2<- as.numeric(as.character(dat$kbit2))

# Interpolating missing values linearly (only 2 per variable)
dat$handedness <- na_interpolation(dat,option = 'linear')$handedness
dat$kbit2 <- na_interpolation(dat,option = 'linear')$kbit2
dat$rt <- na_interpolation(dat,option = 'linear')$rt


# Find average RT per condition
aggregate(dat$rt, list(dat$cond), FUN=mean) 

# Centering reaction time by cond
dat$rtc <- ifelse(dat$cond == 0, dat$rt - 2.32683, 
                  ifelse(dat$cond == 1, dat$rt - 2.256074, 
                         ifelse(dat$cond == 2, dat$rt - 3.638021, dat$rt)))

# Turn from character to numeric
# Loop through each column
for (col in names(dat)) {
  dat[[col]] <- as.numeric(dat[[col]])
}


library(data.table)
library(mice)
# Here is where you would input your "save data set" command of choice. I wasn't 
# sure if .txt or .csv would be moreinit = mice(dat[2:nrow(dat)], maxit = 0) 
init = mice(dat[8:13], maxit = 0) 
meth = init$method
predM = init$predictorMatrix

# Setting variables to not be considered as a predictors during imputation
# In this case, we don't want timestamps to be factored into the calculation for imputing EEG data
predM[1:7] <- 0

# Setting method for imputing (using "norm" for normal imputation)
cols_to_impute <- names(dat)[8:13]

# Use lapply to set the imputation method for each column
meth[cols_to_impute] <- lapply(meth[cols_to_impute], function(x) "norm")

set.seed(103)
imputed = mice(dat, method = meth, predictorMatrix = predM, m = 5)

dat_imp <- complete(imputed)
dat_imp$cond_b <- dat_imp$cond
dat_imp$cond<-as.factor(dat_imp$cond)
dat_imp$cond <- factor(dat_imp$cond, levels = c("0", "1", "2"), labels = c("No Error", "Easy Error", "Hard Error"))

dat_imp$acc_og <- dat_imp$acc
dat_imp$acc <- asin(sqrt(dat_imp$acc))


########################################################################################
#####

## RQ1 ###########################################

# Define the output directory and file name
output_dir <- "C:\\Users\\anama\\OneDrive\\Documents\\DMPM_Paper"
output_file <- file.path(output_dir, "RQ1_Models_output_sl.txt")

# Redirect output to a text file
sink(output_file)

# Model 1 DM: Unconditional Model
cat("\nModel 1 DM: Unconditional Model\n")
mod1dmL <- lmer(dmL ~ 1 + (1 | sub), REML = FALSE, data = dat_imp)
print(summary(mod1dmL))

dmLicc <- quantile(bootMer(mod1dmL, calc.icc, nsim = 10000)$t, c(0.025, 0.975))
cat("\nModel 1 DM - dmL ICC:", dmLicc, "\n")

# Perform the bootstrapping
boot_intercept <- bootMer(mod1dmL, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("DML Random Intercept Term CI:", intercept_ci, "\n")


mod1dmR <- lmer(dmR ~ 1 + (1 | sub), REML = FALSE, data = dat_imp)
print(summary(mod1dmR))

dmRicc <- quantile(bootMer(mod1dmR, calc.icc, nsim = 10000)$t, c(0.025, 0.975))
cat("\nModel 1 DM - dmR ICC:", dmRicc, "\n")


# Perform the bootstrapping
boot_intercept <- bootMer(mod1dmR, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("DMR Random Intercept Term CI:", intercept_ci, "\n")




# Model 2 DM: Level-1 Model
cat("\nModel 2 DM: Level-1 Model\n")
mod2dmL <- lmer(dmL ~ cond + (1  | sub), REML = FALSE, data = dat_imp)
print(summary(mod2dmL))

# Perform the bootstrapping
boot_intercept <- bootMer(mod2dmL, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 dmL '1 | sub':", intercept_ci, "\n")



mod2dmR <- lmer(dmR ~ cond + (1 | sub), REML = FALSE, data = dat_imp)
print(summary(mod2dmR))

# Perform the bootstrapping
boot_intercept <- bootMer(mod2dmR, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 dmR '1 | sub':", intercept_ci, "\n")



# Likelihood Ratio Test for DM Models
cat("Likelihood Ratio Test for DM Models\n")
lrt_dmL <- lrtest(mod1dmL, mod2dmL)
print(lrt_dmL)

lrt_dmR <- lrtest(mod1dmR, mod2dmR)
print(lrt_dmR)



# Model 3 DM: Level-1 and Level-2 Model
cat("\nModel 3 DM: Level-1 and Level-2 Model\n")
mod3dmL <- lmer(dmL ~ cond + handedness + female + (1 | sub), REML = FALSE, data = dat_imp)
print(summary(mod3dmL))

# Perform the bootstrapping
boot_intercept <- bootMer(mod3dmL, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI:", intercept_ci, "\n")




mod3dmR <- lmer(dmR ~ cond + handedness + female + (1  | sub), REML = FALSE, data = dat_imp)
print(summary(mod3dmR))

# Perform the bootstrapping
boot_intercept <- bootMer(mod3dmR, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI:", intercept_ci, "\n")



# Likelihood Ratio Test for DM Models
cat("Likelihood Ratio Test for DM Models\n")
lrt_dmL_3 <- lrtest(mod2dmL, mod3dmL)
print(lrt_dmL_3)

lrt_dmR_3 <- lrtest(mod2dmR, mod3dmR)
print(lrt_dmR_3)


### Checking Change in Model fit with Random Slopes for Best-Fitting DM Models ####

# Model 2 DM: Level-1 Model
cat("\nModel 2 DM: Level-1 Model\n")
mod2dmL_s <- lmer(dmL ~ cond + (1 + cond| sub), REML = FALSE, data = dat_imp)
print(summary(mod2dmL_s))

# Perform the bootstrapping
boot_intercept <- bootMer(mod2dmL_s, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 dmL '1 | sub':", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(mod2dmL_s, extract_random_slope, nsim = 10000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope Model 2 dmL 'condition | sub':", slope_ci, "\n")

cat("Likelihood Ratio Test for DM Models\n")
lrt_dmL_s <- lrtest(mod2dmL, mod2dmL_s)
print(lrt_dmL_s)




# Model 1 PM: Unconditional Model
cat("\nModel 1 PM: Unconditional Model\n")
mod1pmL <- lmer(pmL ~ 1 + (1 | sub), REML = FALSE, data = dat_imp)
print(summary(mod1pmL))

pmLicc <- quantile(bootMer(mod1pmL, calc.icc, nsim = 10000)$t, c(0.025, 0.975))
cat("\nModel 1 PM - pmL ICC:", pmLicc, "\n")
calc.icc(mod1pmL)

# Perform the bootstrapping
boot_intercept <- bootMer(mod1pmL, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI:", intercept_ci, "\n")


mod1pmR <- lmer(pmR ~ 1 + (1 | sub), REML = FALSE, data = dat_imp)
print(summary(mod1pmR))

pmRicc <- quantile(bootMer(mod1pmR, calc.icc, nsim = 10000)$t, c(0.025, 0.975))
cat("\nModel 1 PM - pmR ICC:", pmRicc, "\n")
calc.icc(mod1pmR)

# Perform the bootstrapping
boot_intercept <- bootMer(mod1pmR, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI:", intercept_ci, "\n")





# Model 2 PM: Level-1 Model
cat("\nModel 2 PM: Level-1 Model\n")
mod2pmL <- lmer(pmL ~ cond + (1 | sub), REML = FALSE, data = dat_imp)
print(summary(mod2pmL))

# Perform the bootstrapping
boot_intercept <- bootMer(mod2pmL, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 pmL '1 | sub':", intercept_ci, "\n")


mod2pmR <- lmer(pmR ~ cond + (1 | sub), REML = FALSE, data = dat_imp)
print(summary(mod2pmR))

# Perform the bootstrapping
boot_intercept <- bootMer(mod2pmR, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 pmR '1 | sub':", intercept_ci, "\n")



# Likelihood Ratio Test for PM Models
cat("Likelihood Ratio Test for PM Models\n")
lrt_pmL <- lrtest(mod1pmL, mod2pmL)
print(lrt_pmL)

lrt_pmR <- lrtest(mod1pmR, mod2pmR)
print(lrt_pmR)



# Model 3 PM: Level-1 and Level-2 Model
cat("\nModel 3 PM: Level-1 and Level-2 Model\n")
mod3pmL <- lmer(pmL ~ cond + handedness + female + (1  | sub), REML = FALSE, data = dat_imp)
print(summary(mod3pmL))

# Perform the bootstrapping
boot_intercept <- bootMer(mod3pmL, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI:", intercept_ci, "\n")



mod3pmR <- lmer(pmR ~ cond + handedness + female + (1  | sub), REML = FALSE, data = dat_imp)
print(summary(mod3pmR))

# Perform the bootstrapping
boot_intercept <- bootMer(mod3pmR, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI:", intercept_ci, "\n")



# Likelihood Ratio Test for PM Models
cat("Likelihood Ratio Test for PM Models\n")
lrt_pmL <- lrtest(mod2pmL, mod3pmL)
print(lrt_pmL)

lrt_pmR <- lrtest(mod2pmR, mod3pmR)
print(lrt_pmR)


### Checking Change in Model fit with Random Slopes for Applicable Best-Fitting PM Models ####
### Right PM Best fitting was unconditional 

# Model 2 PM: Level-1 Model
cat("\nModel 2 PM: Level-1 Model\n")
mod2pmL_s <- lmer(pmL ~ cond + (1 + cond| sub), REML = FALSE, data = dat_imp)
print(summary(mod2pmL_s))

# Perform the bootstrapping
boot_intercept <- bootMer(mod2pmL_s, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 pmL '1 | sub':", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(mod2pmL_s, extract_random_slope, nsim = 10000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope Model 2 pmL 'condition | sub':", slope_ci, "\n")

cat("Likelihood Ratio Test for PM Models\n")
lrt_pmL_s <- lrtest(mod2pmL, mod2pmL_s)
print(lrt_pmL_s)




### Checking Benefits of Random Effects
## Estimating fixed (f) models
mod2dmL_f <- lm(dmL ~ cond, data = dat_imp)
mod2dmL_f_summary <- summary(mod2dmL_f)

cat("Fixed Effect - dmL\n")
print(mod1dmL_f_summary)

mod1dmR_f <- lm(dmR ~ 1, data = dat_imp)
mod1dmR_f_summary <- summary(mod1dmR_f)

cat("Fixed Effect  -  dmR\n")
print(mod1dmR_f_summary)

mod2pmL_f <- lm(pmL ~ cond, data = dat_imp)
mod2pmL_f_summary <- summary(mod2pmL_f)

cat("Fixed Effect -  pmL\n")
print(mod2pmL_f_summary)

mod2pmR_f <- lm(pmR~ cond, data = dat_imp)
mod2pmR_f_summary <- summary(mod1pmR_f)

cat("Fixed Effect - pmR\n")
print(mod2pmR_f_summary)


lrt_pmL_f <- lrtest(mod2pmL_s, mod2pmL_f)
print(lrt_pmL_f)

lrt_pmR_f <- lrtest(mod1pmR, mod1pmR_f)
print(lrt_pmR_f)

lrt_dmL_f <- lrtest(mod2dmL_s, mod2dmL_f)
print(lrt_dmL_f)

lrt_dmR_f <- lrtest(mod1dmR, mod1dmR_f)
print(lrt_dmR_f)



# AIC/BIC comparisons of  Fixed vs Random Effect Models

### DML ###
# Calculate AIC and BIC for both models
aic_lm <- AIC(mod2dmL_f)
aic_lmer <- AIC(mod2dmL_s)
bic_lm <- BIC(mod2dmL_f)
bic_lmer <- BIC(mod2dmL_s)

# Print model name
cat("\n########## Model: Left DM ##########\n")

# Print results
cat("AIC for Fixed model 1 dmL", aic_lm, "\n")
cat("AIC for Mixed model 1 dmL", aic_lmer, "\n")
cat("BIC for Fixed model 1 dmL", bic_lm, "\n")
cat("BIC for Mixed model 1 dmL", bic_lmer, "\n")

# Calculate the difference in AIC
aic_diff <- abs(aic_lm - aic_lmer)

# Categorize and print aic_diff
if (aic_diff > 2) {
  cat("The difference in AIC is", aic_diff, "indicating a significant improvement in model fit for the model with the lower AIC.\n")
} else {
  cat("The difference in AIC is", aic_diff, "indicating no significant improvement in model fit.\n")
}

# Calculate the difference in BIC
bic_diff <- bic_lm - bic_lmer

# Categorize the strength of evidence and print bic_diff
if (bic_diff > 10) {
  cat("BIC diff =", bic_diff, "- very strong evidence for the more complex model\n")
} else if (bic_diff > 6) {
  cat("BIC diff =", bic_diff, "- strong evidence for the more complex model\n")
} else if (bic_diff > 2) {
  cat("BIC diff =", bic_diff, "- positive evidence for the more complex model\n")
} else if (bic_diff > 0) {
  cat("BIC diff =", bic_diff, "- weak evidence for the more complex model\n")
} else {
  cat("BIC diff =", bic_diff, "- no evidence for the more complex model\n")
}


### DMR ###
# Calculate AIC and BIC for both models
aic_lm <- AIC(mod1dmR_f)
aic_lmer <- AIC(mod1dmR)
bic_lm <- BIC(mod1dmR_f)
bic_lmer <- BIC(mod1dmR)

# Print model name
cat("\n########## Model: Right DM ##########\n")

# Print results
cat("AIC for Fixed model 1 dmR", aic_lm, "\n")
cat("AIC for Mixed model 1 dmR", aic_lmer, "\n")
cat("BIC for Fixed model 1 dmR", bic_lm, "\n")
cat("BIC for Mixed model 1 dmR", bic_lmer, "\n")

# Calculate the difference in AIC
aic_diff <- abs(aic_lm - aic_lmer)

# Categorize and print aic_diff
if (aic_diff > 2) {
  cat("The difference in AIC is", aic_diff, "indicating a significant improvement in model fit for the model with the lower AIC.\n")
} else {
  cat("The difference in AIC is", aic_diff, "indicating no significant improvement in model fit.\n")
}

# Calculate the difference in BIC
bic_diff <- bic_lm - bic_lmer

# Categorize the strength of evidence and print bic_diff
if (bic_diff > 10) {
  cat("BIC diff =", bic_diff, "- very strong evidence for the more complex model\n")
} else if (bic_diff > 6) {
  cat("BIC diff =", bic_diff, "- strong evidence for the more complex model\n")
} else if (bic_diff > 2) {
  cat("BIC diff =", bic_diff, "- positive evidence for the more complex model\n")
} else if (bic_diff > 0) {
  cat("BIC diff =", bic_diff, "- weak evidence for the more complex model\n")
} else {
  cat("BIC diff =", bic_diff, "- no evidence for the more complex model\n")
}


### PML ###
# Calculate AIC and BIC for both models
aic_lm <- AIC(mod2pmL_f)
aic_lmer <- AIC(mod2pmL_s)
bic_lm <- BIC(mod2pmL_f)
bic_lmer <- BIC(mod2pmL_s)

# Print model name
cat("\n########## Model: Left PM ##########\n")

# Print results
cat("AIC for Fixed model 2 pmL:", aic_lm, "\n")
cat("AIC for Mixed model 2 pmL:", aic_lmer, "\n")
cat("BIC for Fixed model 2 pmL:", bic_lm, "\n")
cat("BIC for Mixed model 2 pmL:", bic_lmer, "\n")

# Calculate the difference in AIC
aic_diff <- abs(aic_lm - aic_lmer)

# Categorize and print aic_diff
if (aic_diff > 2) {
  cat("The difference in AIC is", aic_diff, "indicating a significant improvement in model fit for the model with the lower AIC.\n")
} else {
  cat("The difference in AIC is", aic_diff, "indicating no significant improvement in model fit.\n")
}

# Calculate the difference in BIC
bic_diff <- bic_lm - bic_lmer

# Categorize the strength of evidence and print bic_diff
if (bic_diff > 10) {
  cat("BIC diff =", bic_diff, "- very strong evidence for the more complex model\n")
} else if (bic_diff > 6) {
  cat("BIC diff =", bic_diff, "- strong evidence for the more complex model\n")
} else if (bic_diff > 2) {
  cat("BIC diff =", bic_diff, "- positive evidence for the more complex model\n")
} else if (bic_diff > 0) {
  cat("BIC diff =", bic_diff, "- weak evidence for the more complex model\n")
} else {
  cat("BIC diff =", bic_diff, "- no evidence for the more complex model\n")
}





### PMR ###
# Calculate AIC and BIC for both models
aic_lm <- AIC(mod2pmR_f)
aic_lmer <- AIC(mod2pmR)
bic_lm <- BIC(mod2pmR_f)
bic_lmer <- BIC(mod2pmR)

# Print model name
cat("\n########## Model: Right PM ##########\n")

# Print results
cat("AIC for Fixed model 1 pmR:", aic_lm, "\n")
cat("AIC for Mixed model 1 pmR:", aic_lmer, "\n")
cat("BIC for Fixed model 1 pmR:", bic_lm, "\n")
cat("BIC for Mixed model 1 pmR:", bic_lmer, "\n")

# Calculate the difference in AIC
aic_diff <- abs(aic_lm - aic_lmer)

# Categorize and print aic_diff
if (aic_diff > 2) {
  cat("The difference in AIC is", aic_diff, "indicating a significant improvement in model fit for the model with the lower AIC.\n")
} else {
  cat("The difference in AIC is", aic_diff, "indicating no significant improvement in model fit.\n")
}

# Calculate the difference in BIC
bic_diff <- bic_lm - bic_lmer

# Categorize the strength of evidence and print bic_diff
if (bic_diff > 10) {
  cat("BIC diff =", bic_diff, "- very strong evidence for the more complex model\n")
} else if (bic_diff > 6) {
  cat("BIC diff =", bic_diff, "- strong evidence for the more complex model\n")
} else if (bic_diff > 2) {
  cat("BIC diff =", bic_diff, "- positive evidence for the more complex model\n")
} else if (bic_diff > 0) {
  cat("BIC diff =", bic_diff, "- weak evidence for the more complex model\n")
} else {
  cat("BIC diff =", bic_diff, "- no evidence for the more complex model\n")
}




# Stop redirecting output to the text file
sink()





### RQ2 ###########################################
#### Distribution of Level-1 Across conds ####
# Function to calculate mean and standard error
mean_se <- function(x) {
  mean_val <- mean(x, na.rm = TRUE)
  se_val <- sd(x, na.rm = TRUE) / sqrt(length(x))
  return(c(y = mean_val, ymin = mean_val - se_val, ymax = mean_val + se_val))
}

# accplot
accplot <- ggplot( dat_imp, aes(cond, acc)) +
  stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.5) +
  ylab("Accuracy (%)") +
  geom_smooth(method = lm, se = TRUE, fill = "gray") +  # Add shaded confidence interval
  xlab("Condition") +
  theme_minimal() +  # Use the default theme
  theme(panel.grid.major.x = element_blank(),
       panel.grid.minor.x = element_blank())  # Remove vertical gridlines# Remove vertical gridlines

# rtplot
rtplot <- ggplot( dat_imp, aes(cond, rtc)) +
  stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.5) +
  geom_smooth(method = lm, se = TRUE, fill = "gray") +  # Add shaded confidence interval
  xlab("Condition") +
  theme_minimal()+ # Use the default theme
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())  # Remove vertical gridlines

library(stringr)
y <- c("Reaction Time (Difference from Condition Mean [s])")
wrapped_title <- str_wrap(y, width = 20)  # Wrap at the desired width

# Add the wrapped title to the plot
rtplot <- rtplot + ylab(wrapped_title)

# Use grid.arrange from gridExtra to arrange the plots
grid.arrange(accplot, rtplot)

### START MODELS



#### RQ2 ###########################################

# Set the output directory and file name
output_dir <- "C:/Users/anama/OneDrive/Documents/DMPM_Paper"
output_file <- file.path(output_dir, "RQ2_Models_output_sl.txt")

# Redirect output to a text file
sink(output_file)

#### Model 1: Unconditional Model ####
mod1acc <- lmer(acc ~ 1 + (1 | sub), REML = FALSE, data = dat_imp)
mod1acc_summary <- summary(mod1acc)
icc <- calc.icc(mod1acc)
boot.icc <- bootMer(mod1acc, calc.icc, nsim = 10000)

accicc <- quantile(boot.icc$t, c(0.025, 0.975), na.rm = TRUE)

# Perform the bootstrapping
boot_intercept <- bootMer(mod1acc, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI:", intercept_ci, "\n")



cat("\n Model 1: Unconditional Model\n")
print(mod1acc_summary)
cat("\n Model 1: ICC and CI\n")
cat("ICC: ", icc, "\n95% CI: ", accicc, "\n")
cat("\n Model 1: Table\n")
print(tab_model(mod1acc))



### Mod 2: Random Intercept Only ##############

mod2accdmL_i <- lmer(acc ~ dmL * cond + (1 | sub), REML = FALSE, data = dat_imp)
mod2accdmL_i_summary <- summary(mod2accdmL_i)

cat("Random Intercept Only - Accuracy dmL\n")
print(mod2accdmL_i_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accdmL_i, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 accdmL Intercept Only '1 | sub':", intercept_ci, "\n")




mod2accdmR_i <- lmer(acc ~ dmR * cond + (1 | sub), REML = FALSE, data = dat_imp)
mod2accdmR_i_summary <- summary(mod2accdmR_i)

cat("Random Intercept Only - Accuracy dmR\n")
print(mod2accdmR_i_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accdmR_i, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 accdmR Intercept Only '1 | sub':", intercept_ci, "\n")




mod2accpmL_i <- lmer(acc ~ pmL * cond + (1 | sub), REML = FALSE, data = dat_imp)
mod2accpmL_i_summary <- summary(mod2accpmL_i)

cat("Random Intercept Only - Accuracy pmL\n")
print(mod2accpmL_i_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accpmL_i, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 accpmL Intercept Only '1 | sub':", intercept_ci, "\n")




mod2accpmR_i <- lmer(acc ~ pmR * cond + (1 | sub), REML = FALSE, data = dat_imp)
mod2accpmR_i_summary <- summary(mod2accpmR_i)

cat("Random Intercept Only - Accuracy pmR\n")
print(mod2accpmR_i_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accpmR_i, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 accpmR Intercept Only '1 | sub':", intercept_ci, "\n")

# Likelihood ratio tests for Model 2
lrtest_dmL <- lrtest(mod1acc, mod2accdmL_i)
lrtest_dmR <- lrtest(mod1acc, mod2accdmR_i)
lrtest_pmL <- lrtest(mod1acc, mod2accpmL_i)
lrtest_pmR <- lrtest(mod1acc, mod2accpmR_i)

cat("LR Tests for Model 2\n")
print(lrtest_dmL)
print(lrtest_dmR)
print(lrtest_pmL)
print(lrtest_pmR)



#### Model 2: Level-1 Model Linear, Random Slopes + Intercepts ####
mod2accdmL <- lmer(acc_og ~ dmL * cond + (1 + dmL | sub), REML = FALSE, data = dat_imp)
mod2accdmL_summary <- summary(mod2accdmL)

cat("\n Model 2: Level-1 Model Linear - Accuracy dmL\n")
print(mod2accdmL_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accdmL, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 accdmL Accuracy '1 | sub':", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(mod2accdmL, extract_random_slope, nsim = 10000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope Model 2 accdmL Accuracy 'condition | sub':", slope_ci, "\n")



mod2accdmR <- lmer(acc ~ dmR * cond + (1 + dmR | sub), REML = FALSE, data = dat_imp)
mod2accdmR_summary <- summary(mod2accdmR)

cat("\n Model 2: Level-1 Model Linear - dmR\n")
print(mod2accdmR_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accdmR, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 accdmR '1 | sub':", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(mod2accdmR, extract_random_slope, nsim = 10000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope Model 2 accdmR 'condition | sub':", slope_ci, "\n")




mod2accpmL <- lmer(acc_og ~ pmL * cond + (1 + pmL | sub), REML = FALSE, data = dat_imp)
mod2accpmL_summary <- summary(mod2accpmL)

cat("\n Model 2: Level-1 Model Linear - pmL\n")
print(mod2accpmL_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accpmL, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 accpmL '1 | sub':", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(mod2accpmL, extract_random_slope, nsim = 10000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope Model 2 accpmL 'condition | sub':", slope_ci, "\n")




mod2accpmR <- lmer(acc ~ pmR * cond + (1 + pmR | sub), REML = FALSE, data = dat_imp)
mod2accpmR_summary <- summary(mod2accpmR)

cat("\n Model 2: Level-1 Model Linear - pmR\n")
print(mod2accpmR_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accpmR, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 accpmR '1 | sub':", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(mod2accpmR, extract_random_slope, nsim = 10000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope Model 2 accpmR 'condition | sub':", slope_ci, "\n")



# Likelihood ratio tests for Model 3
lrtest_dmL_i <- lrtest(mod2accdmL_i, mod2accdmL)
lrtest_dmR_i <- lrtest(mod2accdmR_i, mod2accdmR)
lrtest_pmL_i <- lrtest(mod2accpmL_i, mod2accpmL)
lrtest_pmR_i <- lrtest(mod2accpmR_i, mod2accpmR)

cat("LR Tests for Intercept-Only vs Intercept + Slope\n")
print(lrtest_dmL_i)
print(lrtest_dmR_i)
print(lrtest_pmL_i)
print(lrtest_pmR_i)


### Checking Benefits of Random Effects
mod2accdmL_f <- lm(acc ~ dmL * cond, data = dat_imp)
mod2accdmL_f_summary <- summary(mod2accdmL_f)

cat("Random Slopes - Accuracy dmL\n")
print(mod2accdmL_f_summary)

mod2accdmR_f <- lm(acc ~ dmR * cond, data = dat_imp)
mod2accdmR_f_summary <- summary(mod2accdmR_f)

cat("Random Slopes - Accuracy dmR\n")
print(mod2accdmR_f_summary)

mod2accpmL_f <- lm(acc ~ pmL * cond, data = dat_imp)
mod2accpmL_f_summary <- summary(mod2accpmL_f)

cat("Random Slopes - Accuracy pmL\n")
print(mod2accpmL_f_summary)

mod2accpmR_f <- lm(acc ~ pmR * cond, data = dat_imp)
mod2accpmR_f_summary <- summary(mod2accpmR_f)

cat("Random Slopes - Accuracy pmR\n")
print(mod2accpmR_f_summary)




lrt_pmL_f <- lrtest(mod2accpmL_i, mod2accpmL_f)
print(lrt_pmL_f)

lrt_pmR_f <- lrtest(mod2accpmR_i, mod2accpmR_f)
print(lrt_pmR_f)

lrt_dmL_f <- lrtest(mod2accdmL, mod2accdmL_f)
print(lrt_dmL_f)

lrt_dmR_f <- lrtest(mod2accdmR_i, mod2accdmR_f)
print(lrt_dmR_f)




#### Model 3: Level-1 and Level-2 Model Random Intercept Only ####
mod3accdmL <- lmer(acc ~ dmL * cond + handedness + female + (1 + dmL | sub), REML = FALSE, data = dat_imp)
mod3accdmL_summary <- summary(mod3accdmL)

cat("\n Model 3: Level-1 and Level-2 - dmL\n")
print(mod3accdmL_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod3accdmL, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI:", intercept_ci, "\n")


# Perform the bootstrapping
boot_slope <- bootMer(mod3accdmL, extract_random_slope, nsim = 10000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope:", slope_ci, "\n")




mod3accdmR <- lmer(acc ~ dmR * cond + handedness + female + (1| sub), REML = FALSE, data = dat_imp)
mod3accdmR_summary <- summary(mod3accdmR)

cat("\n Model 3: Level-1 and Level-2 - dmR\n")
print(mod3accdmR_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod3accdmR, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI:", intercept_ci, "\n")




mod3accpmL <- lmer(acc ~ pmL * cond + handedness + female + (1 + pmL| sub), REML = FALSE, data = dat_imp)
mod3accpmL_summary <- summary(mod3accpmL)

cat("\n Model 3: Level-1 and Level-2 - pmL\n")
print(mod3accpmL_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod3accpmL, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI:", intercept_ci, "\n")




mod3accpmR <- lmer(acc ~ pmR * cond + handedness + female + (1| sub), REML = FALSE, data = dat_imp)
mod3accpmR_summary <- summary(mod3accpmR)

cat("\n Model 3: Level-1 and Level-2 - pmR\n")
print(mod3accpmR_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod3accpmR, extract_random_intercept, nsim = 10000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI", intercept_ci, "\n")


# Likelihood ratio tests for Model 3
lrtest_dmL_3 <- lrtest(mod2accdmL, mod3accdmL)
lrtest_dmR_3 <- lrtest(mod2accdmR, mod3accdmR)
lrtest_pmL_3 <- lrtest(mod2accpmL, mod3accpmL)
lrtest_pmR_3 <- lrtest(mod2accpmR, mod3accpmR)

cat("LR Tests for Model 3\n")
print(lrtest_dmL_3)
print(lrtest_dmR_3)
print(lrtest_pmL_3)
print(lrtest_pmR_3)


# AIC/BIC comparisons of  Fixed vs Random Effect Models

### DML Acc ###
# Calculate AIC and BIC for both models
aic_lm <- AIC(mod2accdmL_f)
aic_lmer <- AIC(mod2accdmL)
bic_lm <- BIC(mod2accdmL_f)
bic_lmer <- BIC(mod2accdmL)

# Print model name
cat("\n########## Model: DML Acc ##########\n")

# Print results
cat("AIC for Fixed model:", aic_lm, "\n")
cat("AIC for Mixed model:", aic_lmer, "\n")
cat("BIC for Fixed model:", bic_lm, "\n")
cat("BIC for Mixed model:", bic_lmer, "\n")

# Calculate the difference in AIC
aic_diff <- abs(aic_lm - aic_lmer)

# Categorize and print aic_diff
if (aic_diff > 2) {
  cat("The difference in AIC is", aic_diff, "indicating a significant improvement in model fit for the model with the lower AIC.\n")
} else {
  cat("The difference in AIC is", aic_diff, "indicating no significant improvement in model fit.\n")
}

# Calculate the difference in BIC
bic_diff <- bic_lm - bic_lmer

# Categorize the strength of evidence and print bic_diff
if (bic_diff > 10) {
  cat("BIC diff =", bic_diff, "- very strong evidence for the more complex model\n")
} else if (bic_diff > 6) {
  cat("BIC diff =", bic_diff, "- strong evidence for the more complex model\n")
} else if (bic_diff > 2) {
  cat("BIC diff =", bic_diff, "- positive evidence for the more complex model\n")
} else if (bic_diff > 0) {
  cat("BIC diff =", bic_diff, "- weak evidence for the more complex model\n")
} else {
  cat("BIC diff =", bic_diff, "- no evidence for the more complex model\n")
}




### DMR Acc###
# Calculate AIC and BIC for both models
aic_lm <- AIC(mod2accdmR_f)
aic_lmer <- AIC(mod2accdmR_i)
bic_lm <- BIC(mod2accdmR_f)
bic_lmer <- BIC(mod2accdmR_i)

# Print model name
cat("\n########## Model: DMR Acc ##########\n")

# Print results
cat("AIC for Fixed model:", aic_lm, "\n")
cat("AIC for Mixed model:", aic_lmer, "\n")
cat("BIC for Fixed model:", bic_lm, "\n")
cat("BIC for Mixed model:", bic_lmer, "\n")

# Calculate the difference in AIC
aic_diff <- abs(aic_lm - aic_lmer)

# Categorize and print aic_diff
if (aic_diff > 2) {
  cat("The difference in AIC is", aic_diff, "indicating a significant improvement in model fit for the model with the lower AIC.\n")
} else {
  cat("The difference in AIC is", aic_diff, "indicating no significant improvement in model fit.\n")
}

# Calculate the difference in BIC
bic_diff <- bic_lm - bic_lmer

cat("BIC diff = ", bic_diff)

# Categorize the strength of evidence and print bic_diff
if (bic_diff > 10) {
  cat("BIC diff =", bic_diff, "- very strong evidence for the more complex model\n")
} else if (bic_diff > 6) {
  cat("BIC diff =", bic_diff, "- strong evidence for the more complex model\n")
} else if (bic_diff > 2) {
  cat("BIC diff =", bic_diff, "- positive evidence for the more complex model\n")
} else if (bic_diff > 0) {
  cat("BIC diff =", bic_diff, "- weak evidence for the more complex model\n")
} else {
  cat("BIC diff =", bic_diff, "- no evidence for the more complex model\n")
}





### PML ###
# Calculate AIC and BIC for both models
aic_lm <- AIC(mod2accpmL_f)
aic_lmer <- AIC(mod2accpmL_i)
bic_lm <- BIC(mod2accpmL_f)
bic_lmer <- BIC(mod2accpmL_i)

# Print model name
cat("\n########## Model: PML##########\n")

# Print results
cat("AIC for Fixed model:", aic_lm, "\n")
cat("AIC for Mixed model:", aic_lmer, "\n")
cat("BIC for Fixed model:", bic_lm, "\n")
cat("BIC for Mixed model:", bic_lmer, "\n")

# Calculate the difference in AIC
aic_diff <- abs(aic_lm - aic_lmer)

# Categorize and print aic_diff
if (aic_diff > 2) {
  cat("The difference in AIC is", aic_diff, "indicating a significant improvement in model fit for the model with the lower AIC.\n")
} else {
  cat("The difference in AIC is", aic_diff, "indicating no significant improvement in model fit.\n")
}

# Calculate the difference in BIC
bic_diff <- bic_lm - bic_lmer

# Categorize the strength of evidence and print bic_diff
if (bic_diff > 10) {
  cat("BIC diff =", bic_diff, "- very strong evidence for the more complex model\n")
} else if (bic_diff > 6) {
  cat("BIC diff =", bic_diff, "- strong evidence for the more complex model\n")
} else if (bic_diff > 2) {
  cat("BIC diff =", bic_diff, "- positive evidence for the more complex model\n")
} else if (bic_diff > 0) {
  cat("BIC diff =", bic_diff, "- weak evidence for the more complex model\n")
} else {
  cat("BIC diff =", bic_diff, "- no evidence for the more complex model\n")
}





### Acc PMR ###
# Calculate AIC and BIC for both models
aic_lm <- AIC(mod2accpmR_f)
aic_lmer <- AIC(mod2accpmR_i)
bic_lm <- BIC(mod2accpmR_f)
bic_lmer <- BIC(mod2accpmR_i)

# Print model name
cat("\n########## Model: PMR Accuracy ##########\n")

# Print results
cat("AIC for Fixed model:", aic_lm, "\n")
cat("AIC for Mixed model:", aic_lmer, "\n")
cat("BIC for Fixed model:", bic_lm, "\n")
cat("BIC for Mixed model:", bic_lmer, "\n")

# Calculate the difference in AIC
aic_diff <- abs(aic_lm - aic_lmer)

# Categorize and print aic_diff
if (aic_diff > 2) {
  cat("The difference in AIC is", aic_diff, "indicating a significant improvement in model fit for the model with the lower AIC.\n")
} else {
  cat("The difference in AIC is", aic_diff, "indicating no significant improvement in model fit.\n")
}

# Calculate the difference in BIC
bic_diff <- bic_lm - bic_lmer

# Categorize the strength of evidence and print bic_diff
if (bic_diff > 10) {
  cat("BIC diff =", bic_diff, "- very strong evidence for the more complex model\n")
} else if (bic_diff > 6) {
  cat("BIC diff =", bic_diff, "- strong evidence for the more complex model\n")
} else if (bic_diff > 2) {
  cat("BIC diff =", bic_diff, "- positive evidence for the more complex model\n")
} else if (bic_diff > 0) {
  cat("BIC diff =", bic_diff, "- weak evidence for the more complex model\n")
} else {
  cat("BIC diff =", bic_diff, "- no evidence for the more complex model\n")
}

# Stop redirecting output to the text file
sink()







#### Visualizing Results

tab_model(mod1acc, mod2accpmL,  mod3accpmL)
tab_model(mod1acc, mod2accdmL, mod3accdmL)


## Visualizing Models

tab_model(mod1acc, mod2accdmL, mod3accdmL, show.se = TRUE)
tab_model(mod1acc,mod2accpmL,mod3accpmL, show.se = TRUE)

tab_model(mod1rt,mod2rtdmR, mod2rtdm, mod3rtdm, show.se = TRUE)
tab_model(mod1rt,mod2rtpm, mod2rtpmR,mod3rtpm, show.se = TRUE)



dat_imp$acc_og <- sin(dat_imp$acc)^2

### Interaction Plot Left PM

library(extrafont) 

windowsFonts(Times = windowsFont("Times New Roman"))   

plot_obj <- interactions::interact_plot(mod2accdmL, pred = dmL, modx = cond, interval = TRUE,
                                        int.type = c("confidence", "prediction"),
                                        int.width = 0.95, legend.main = "Condition", colors= "CUD Bright",
                                        vary.lty= TRUE, line.size = 5) 


# Modify the axis labels and legend title using ggplot2 functions
plot_obj +
  ggplot2::labs(x = "Left PM ROI Beta",
                y = "Accuracy (Proportion Correct)",
                color = "Condition") +
  ggplot2::scale_x_continuous(limits = c(-2, 2.5), breaks = seq(-2, 2.5, by = .5)) + # Adjust x-axis limits and breaks
  ggplot2::scale_y_continuous(limits = c(0.55, 1.0), breaks = seq(0.55, 1.0, by = 0.1)) + # Adjust y-axis limits and breaks
  theme(text = element_text(size = 14),
        panel.background = element_rect(fill = "whitesmoke"),  # Set background color
        panel.grid.major = element_line(color = "white"),      # Set color of major gridlines
        panel.grid.minor = element_line(color = "white"))


### Interaction Plot Left DM

windowsFonts(Times = windowsFont("Times New Roman"))   

plot_obj <- interactions::interact_plot(mod2accpmL, pred = pmL, modx = cond, interval = TRUE,
                                        int.type = c("confidence", "prediction"),
                                        int.width = 0.95, legend.main = "Condition", colors= "Rainbow",
                                        vary.lty= TRUE, line.size = 5) 


# Modify the axis labels and legend title using ggplot2 functions
plot_obj +
  ggplot2::labs(x = "Left DM ROI Beta",
                y = "Accuracy (Proportion Correct)",
                color = "Condition") +
  ggplot2::scale_x_continuous(limits = c(-2, 2.5), breaks = seq(-2, 2.5, by = .5)) + # Adjust x-axis limits and breaks
  ggplot2::scale_y_continuous(limits = c(0.55, 1), breaks = seq(0.55, 1, by = 0.1)) + # Adjust y-axis limits and breaks
  theme(text = element_text(size = 14),
        panel.background = element_rect(fill = "whitesmoke"),  # Set background color
        panel.grid.major = element_line(color = "white"),      # Set color of major gridlines
        panel.grid.minor = element_line(color = "white"))




#########################################################################
################ Checking Assumptions for best-fitting models ###########
#########################################################################
dev.off()

#### dm Assumptions ############################################

# INDEPENDENCE (1.5-2.5)
car::dwt(resid(mod2dm)) # Durbin-Watson test for first-order autocorrelation

### HOMOSCEDASTICITY
plot(mod2dm) # Level-1 Residuals
anova(
  lm(
    abs(
      resid(mod2dm))^2 ~ sub, data =  dat_imp)) # Homoscedasticity

plot(ranef(mod2dm)) # Level-Two residuals

### NORMALITY
hist(ranef(mod2dm))
# QQ-Plots for fixed & random effects #
lattice::qqmath(mod2dm)
lattice::qqmath(ranef(mod2dm))

### OUTLIERS

psych::describe(resid(mod2dm))
# Influence: Level-2 #
infl1 = hlm_influence(mod2dm, level = "sub")
dotplot_diag(infl1$cooksd, name = "cooks.distance", cutoff = "internal")
barplot(infl1$cooksd, names.arg = infl1$age11, las = 2, cex.names = 0.5)

#### Assessing effect of excluding outliers
# Dataset with 1 outlier removed
dat_o1 =  dat_imp[-14,]

mod2dm_o1 = lmer(dm ~ cond + roi+ (1  + roi|sub), REML = F, data = dat_o1)

summary(mod2dm_o1)

# Calculating Pseudo R^2 to compare model fit for models with outlier removed
tab_model(mod2dm, mod2dm_o1)

## dmR Model Assumptions##############################################

car::dwt(resid(mod2dmR)) # Durbin-Watson test for first-order autocorrelation

### HOMOSCEDASTICITY
plot(mod2dmR) # Level-1 Residuals
anova(
  lm(
    abs(
      resid(mod2dmR))^2 ~ sub, data =  dat_imp)) # Homoscedasticity

plot(ranef(mod2dmR)) # Level-Two residuals

### NORMALITY
hist(ranef(mod2dmR))
# QQ-Plots for fixed & random effects #
lattice::qqmath(mod2dmR)
lattice::qqmath(ranef(mod2dmR))

### OUTLIERS

psych::describe(resid(mod2dmR))
# Influence: Level-2 #
infl1 = hlm_influence(mod2dmR, level = "sub")
dotplot_diag(infl1$cooksd, name = "cooks.distance", cutoff = "internal")
barplot(infl1$cooksd, names.arg = infl1$age11, las = 2, cex.names = 0.5)

#### Assessing effect of excluding outliers
# Dataset with 1 outlier removed
dat_o1 =  dat_imp[-105,]

mod2dmR_o1 = lmer(dmR ~ cond + roi+ (1  + roi|sub), REML = F, data = dat_o1)

summary(mod2dmR_o1)

# Calculating Pseudo R^2 to compare model fit for models with outlier removed
tab_model(mod2dmR, mod2dmR_o1)

## pm Model Assumptions###############################################

car::dwt(resid(mod2pm)) # Durbin-Watson test for first-order autocorrelation

### HOMOSCEDASTICITY
plot(mod2pm) # Level-1 Residuals
anova(
  lm(
    abs(
      resid(mod2pm))^2 ~ sub, data =  dat_imp)) # Homoscedasticity

plot(ranef(mod2pm)) # Level-Two residuals

### NORMALITY
hist(ranef(mod2pm))
# QQ-Plots for fixed & random effects #
lattice::qqmath(mod2pm)
lattice::qqmath(ranef(mod2pm))

### OUTLIERS

psych::describe(resid(mod2pm))
# Influence: Level-2 #
infl1 = hlm_influence(mod2pm, level = "sub")
dotplot_diag(infl1$cooksd, name = "cooks.distance", cutoff = "internal")
barplot(infl1$cooksd, names.arg = infl1$age11, las = 2, cex.names = 0.5)

# None identified

## pmR Model Assumptions#############################

car::dwt(resid(mod2pmR)) # Durbin-Watson test for first-order autocorrelation

### HOMOSCEDASTICITY
plot(mod2pmR) # Level-1 Residuals
anova(
  lm(
    abs(
      resid(mod2pmR))^2 ~ sub, data =  dat_imp)) # Homoscedasticity

plot(ranef(mod2pmR)) # Level-Two residuals

### NORMALITY
hist(ranef(mod2pmR))
# QQ-Plots for fixed & random effects #
lattice::qqmath(mod2pmR)
lattice::qqmath(ranef(mod2pmR))

### OUTLIERS

psych::describe(resid(mod2pmR))
# Influence: Level-2 #
infl1 = hlm_influence(mod2pmR, level = "sub")
dotplot_diag(infl1$cooksd, name = "cooks.distance", cutoff = "internal")
barplot(infl1$cooksd, names.arg = infl1$age11, las = 2, cex.names = 0.5)

#### Assessing effect of excluding outliers
# Dataset with 1 outlier removed
dat_o1 =  dat_imp[-100,]

mod2pmR_o1 = lmer(pmR ~ cond + roi+ (1  + roi|sub), REML = F, data = dat_o1)

# Calculating Pseudo R^2 to compare model fit for models with outlier removed
tab_model(mod2pmR, mod2pmR_o1)


############################################################# RQ2 ###############################

#### Distribution of Level-1 Across conds ####
accplot<- ggplot( dat_imp, aes(cond, acc)) +
  stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.05) +
  ylab("Accuracy") + geom_smooth(method = lm, se = TRUE)

rtplot<- ggplot( dat_imp, aes(cond, rtc)) +
  stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.05) +
  ylab("Reaction Time") + geom_smooth(method = lm, se = TRUE)

grid.arrange(accplot, rtplot)

# Create a histogram comparing the distributions of 'acc' across 'cond'
acc_hist <- ggplot( dat_imp, aes(x = acc, fill = factor(cond))) +
  geom_histogram(binwidth = 0.02) +
  labs(x = "Accuracy", y = "Frequency") +
  scale_fill_discrete(name = "cond") +
  facet_wrap(~cond, ncol = 2) +  # Create separate facets for each level of 'cond'
  theme_minimal()

# Print the histograms
print(acc_hist)


#### Distribution of Level-1 Across Level-2 Units ####
subs<- sample( dat_imp$sub,12, replace=FALSE)

dat_subset <- subset( dat_imp, sub %in% c(subs))

ggplot(dat_subset, aes(x = cond, y = acc)) + 
  geom_point() +
  geom_smooth(method=glm, se = FALSE) +
  facet_wrap(~sub, ncol = 4, scales = "free_x") + 
  theme_bw() + 
  xlab("cond") + 
  ylab("Accuracy")

ggplot(dat_subset, aes(x = cond, y = rt)) + 
  geom_point() +
  geom_smooth(method=glm, se = FALSE) +
  facet_wrap(~sub, ncol = 4, scales = "free_x") + 
  theme_bw() + 
  xlab("cond") + 
  ylab("Reaction Time")

# Assessing assumptions for best-fitting models
## DM Accuracy Model Assumptions##############################################

car::dwt(resid(mod3accdm)) # Durbin-Watson test for first-order autocorrelation

### HOMOSCEDASTICITY
plot(mod3accdm) # Level-1 Residuals
anova(
  lm(
    abs(
      resid(mod3accdm))^2 ~ sub, data =  dat_imp)) # Homoscedasticity

plot(ranef(mod3accdm)) # Level-Two residuals

### NORMALITY
hist(ranef(mod3accdm))
# QQ-Plots for fixed & random effects #
lattice::qqmath(mod3accdm)
lattice::qqmath(ranef(mod3accdm))

### OUTLIERS

psych::describe(resid(mod3accdm))
# Influence: Level-2 #
infl1 = hlm_influence(mod3accdm, level = "sub")
dotplot_diag(infl1$cooksd, name = "cooks.distance", cutoff = "internal")
barplot(infl1$cooksd, names.arg = infl1$age11, las = 2, cex.names = 0.5)

#### Assessing effect of excluding outliers
# Dataset with 1 outlier removed
dat_o1 =  dat_imp[-95,]

mod3accdm_o1 = lmer(acc ~ pm + pmR + cond  + handedness + female  + (1  + pm + pmR|sub), REML = F, data =  dat_imp)

summary(mod3accdm_o1)

# Calculating Pseudo R^2 to compare model fit for models with outlier removed
tab_model(mod3accdm, mod3accdm_o1)

## PM Accuracy Model Assumptions##############################################

car::dwt(resid(mod3accpm)) # Durbin-Watson test for first-order autocorrelation

### HOMOSCEDASTICITY
plot(mod3accpm) # Level-1 Residuals
anova(
  lm(
    abs(
      resid(mod3accpm))^2 ~ sub, data =  dat_imp)) # Homoscedasticity

plot(ranef(mod3accpm)) # Level-Two residuals

### NORMALITY
hist(ranef(mod3accpm))
# QQ-Plots for fixed & random effects #
lattice::qqmath(mod3accpm)
lattice::qqmath(ranef(mod3accpm))

### OUTLIERS

psych::describe(resid(mod3accpm))
# Influence: Level-2 #
infl1 = hlm_influence(mod3accpm, level = "sub")
dotplot_diag(infl1$cooksd, name = "cooks.distance", cutoff = "internal")
barplot(infl1$cooksd, names.arg = infl1$age11, las = 2, cex.names = 0.5)

#### Assessing effect of excluding outliers
# Dataset with 1 outlier removed
dat_o1 =  dat_imp[-22,]

mod3accpm_o1 = lmer(acc ~ pm + pmR + cond  + handedness + female  + (1  + pm + pmR|sub), REML = F, data =  dat_imp)

summary(mod3accpm_o1)

# Calculating Pseudo R^2 to compare model fit for models with outlier removed
tab_model(mod3accpm, mod3accpm_o1)


###### CORRELATION MATRIC ##############
library(corrplot)

cdnam <- colnames(dat_imp)
cdnam <- intersect(cdnam, colnames(dat_imp))

corrdat<- dat_imp

# Loop through the synchronized column names
for (col_name in cdnam) {
  # Extract values from lists and convert to numeric
  corrdat[[col_name]] <- sapply(corrdat[[col_name]], function(x) as.numeric(unlist(x)))
}

corrdat<- corrdat[,c(10,14,11,15,7,8,18,4,5,2,3)]

# Calculate the correlation matrix
matrix <- cor(corrdat)
lower.tri(matrix)

source("http://www.sthda.com/upload/rquery_cormat.r")
rquery.cormat(corrdat)



# Makes a colorful correlation plot to make me happy
correlation_matrix <- corr.test(corrdat)$r

# Extract p-values from the corr.test result
p_values <- corr.test(corrdat)$p

# Create a correlation matrix plot with significance levels
corrplot(correlation_matrix, method = "color", p.mat = p_values, sig.level = 0.05)

# Create a correlation matrix plot with correlation coefficients and p-values
corrplot(
  correlation_matrix,
  method = "number",
  p.mat = p_values,
  sig.level = 0.05,
  type = "lower",  # Display only the lower triangle
  diag = FALSE,    # Do not display diagonal values
)




# Load necessary libraries
library(ggplot2)
library(dplyr)

# Summarize the data to get mean and standard error for pmL and pmR
summary_dat <- dat_imp %>%
  group_by(cond) %>%
  summarise(
    pmL_mean = mean(pmL, na.rm = TRUE),
    pmL_se = sd(pmL, na.rm = TRUE) / sqrt(n()),
    pmR_mean = mean(pmR, na.rm = TRUE),
    pmR_se = sd(pmR, na.rm = TRUE) / sqrt(n())
  )

# Reshape the data to long format for ggplot2
summary_long <- summary_dat %>%
  pivot_longer(cols = c(pmL_mean, pmR_mean), names_to = "hemisphere", values_to = "mean") %>%
  mutate(se = ifelse(hemisphere == "pmL_mean", pmL_se, pmR_se),
         hemisphere = recode(hemisphere, "pmL_mean" = "Left Hemisphere (pmL)", "pmR_mean" = "Right Hemisphere (pmR)"))

# Convert 'cond' to a factor with labels
summary_long$cond <- factor(summary_long$cond, levels = c("No Error", "Easy Error", "Hard Error"),
                            labels = c("No Error", "Easy Error", "Hard Error"))

# Plot pmL and pmR together for each condition, clustered by hemisphere
p <- ggplot(summary_long, aes(x = hemisphere, y = mean, fill = cond)) +
  geom_bar(position = position_dodge(width = 0.8), stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(title = "pmL and pmR Across Conditions",
       x = "Hemisphere", y = "Mean ± SE") +
  scale_fill_manual(values = c("#63c0fa", "#fc9431", "#44ddaa"),
                    name = "Condition",
                    labels = c("No Error", "Easy Error", "Hard Error")) +
  ylim(-0.1, 0.4) +  # Set the y-axis limits 
  theme_minimal()

# Print the plot
print(p)


#####

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Summarize the data to get mean and standard error for dmL and dmR
summary_dat <- dat_imp %>%
  group_by(cond) %>%
  summarise(
    dmL_mean = mean(dmL, na.rm = TRUE),
    dmL_se = sd(dmL, na.rm = TRUE) / sqrt(n()),
    dmR_mean = mean(dmR, na.rm = TRUE),
    dmR_se = sd(dmR, na.rm = TRUE) / sqrt(n())
  )

# Reshape the data to long format for ggplot2
summary_long <- summary_dat %>%
  pivot_longer(cols = c(dmL_mean, dmR_mean), names_to = "hemisphere", values_to = "mean") %>%
  mutate(se = ifelse(hemisphere == "dmL_mean", dmL_se, dmR_se),
         hemisphere = recode(hemisphere, "dmL_mean" = "Left DM", "dmR_mean" = "Right DM"))

# Convert 'cond' to a factor with labels
summary_long$cond <- factor(summary_long$cond, levels = c("No Error", "Easy Error", "Hard Error"),
                            labels = c("No Error", "Easy Error", "Hard Error"))

# Plot dmL and dmR together for each condition, clustered by hemisphere
p <- ggplot(summary_long, aes(x = hemisphere, y = mean, fill = cond)) +
  geom_bar(position = position_dodge(width = 0.8), stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, position = position_dodge(width = 0.8)) +
  labs(title = "dmL and dmR Across Conditions",
       x = "Hemisphere", y = "Mean ± SE") +
  scale_fill_manual(values = c("#8f4596", "#8dc07a", "#dd4241"),
                    name = "Condition",
                    labels = c("No Error", "Easy Error", "Hard Error")) +
  ylim(-0.1, 0.4) +  # Set the y-axis limits 
  theme_minimal()

# Print the plot
print(p)
