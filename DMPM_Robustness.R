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

# libraries for plotting
library(lubridate)
library(ggplot2)
library(scales)
library(gridExtra)
library(ggthemes)
library(gridExtra)
library("Hmisc")


# Importing data for RQ1 (DV= dm/pm, L1= cond, L2= sub)
dat = read.csv("C:/Users/anama/OneDrive/Documents/DMPM_Paper/WB_G.csv")

# Transforming variables to numeric/continuous
dat$handedness<- as.numeric(as.character(dat$handedness))
dat$wj_wid<- as.numeric(as.character(dat$wj_wid))

# Interpolating missing values linearly (only 2 per variable)
dat$handedness <- na_interpolation(dat,option = 'linear')$handedness
dat$wj_wid <- na_interpolation(dat,option = 'linear')$wj_wid
dat$rt <- na_interpolation(dat,option = 'linear')$rt

# Converting accuracy into a percent value
#dat$acc<- dat$acc*100

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

set.seed(103)
init = mice(dat[8:13], maxit = 0) 
meth = init$method
predM = init$predictorMatrix
imputed = mice(dat, method = meth, predictorMatrix = predM, m = 5)


# Setting variables to not be considered as a predictors during imputation
# In this case, we don't want timestamps to be factored into the calculation for imputing EEG data
predM[1:7] <- 0

# Setting method for imputing (using "norm" for normal imputation)
cols_to_impute <- names(dat)[8:13]

# Use lapply to set the imputation method for each column
meth[cols_to_impute] <- lapply(meth[cols_to_impute], function(x) "norm")

wb_dat <- complete(imputed)
wb_dat$cond<-as.factor(wb_dat$cond)
wb_dat$cond <- factor(wb_dat$cond, levels = c("0", "1", "2"), labels = c("No Error", "Easy Error", "Hard Error"))

wb_dat$acc <- asin(sqrt(wb_dat$acc))

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



########################################################################################


wbplot<- ggplot(wb_dat, aes(cond, mean)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "text", aes(label = after_stat(round(y, digits=2))), color = "blue", vjust = -0.5) +
  stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.5) +
  ylab("Whole Brain Beta") 

wbplot

# Individual growth plot
# randomly selecting 12 subject ids to plot
subs<- sample(wb_dat$sub,12, replace=FALSE)

dat_subset <- subset(wb_dat, sub %in% c(subs))

ggplot(dat_subset, aes(x = cond, y = mean)) + 
  geom_point() +
  geom_smooth(method=glm, se = FALSE) +
  facet_wrap(~sub, ncol = 4, scales = "free_x") + 
  theme_bw() + 
  xlab("cond") + 
  ylab("Whole Brain Beta")


## RQ1 ###########################################
#### Model 1 DM: Unconditional Model ####
wbmod1 = lmer(mean ~ 1 + (1|sub), REML = F, data =  wb_dat)

summary(wbmod1)

tab_model(wbmod1)

wbicc <- quantile(bootMer(wbmod1, calc.icc, nsim = 1000)$t, c(0.025, 0.975))
cat("\nModel 1 WB ICC:", wbicc, "\n")

# Perform the bootstrapping
boot_intercept <- bootMer(wbmod1, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for wbmod1 | sub':", intercept_ci, "\n")



#### Model 2 DM: Level-1 Model  ####
wbmod2mean = lmer(mean ~ cond+  (1|sub), REML = F, data =  wb_dat)

summary(wbmod2mean)

# Perform the bootstrapping
boot_intercept <- bootMer(wbmod2mean, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 dmL '1 | sub':", intercept_ci, "\n")


#### Model 3 DM: Level-1 and Level-2 Model  ####
wbmod3mean = lmer(mean ~ cond + handedness + female  +  (1 |sub), REML = F, data =  wb_dat)

summary(wbmod3mean)


# Perform the bootstrapping
boot_intercept <- bootMer(wbmod3mean, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 dmL '1 | sub':", intercept_ci, "\n")



lrtest(wbmod1,wbmod2mean)
lrtest(wbmod2mean,wbmod3mean)

tab_model(wbmod1, wbmod2mean, wbmod3mean)


### RQ2 ###########################################

#### Model 1: Unconditional Model ####
wbmod1acc = lmer(acc~ 1 + (1|sub), REML = F, data =  wb_dat)

summary(wbmod1acc)

tab_model(wbmod1acc)

wbicc <- quantile(bootMer(wbmod1acc, calc.icc, nsim = 1000)$t, c(0.025, 0.975), na.rm= TRUE)
cat("\nModel 1 WB ICC:", wbicc, "\n")

# Perform the bootstrapping
boot_intercept <- bootMer(wbmod1acc, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for wbmod1acc | sub':", intercept_ci, "\n")



#### Model 2: Level-1 Model Linear ####

wbmod2acc = lmer(acc ~ mean*cond+ (1|sub), REML = F, data =  wb_dat)

summary(wbmod2acc)

# Perform the bootstrapping
boot_intercept <- bootMer(wbmod1acc, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for wbmod1acc | sub':", intercept_ci, "\n")




wbmod2acc_s = lmer(acc ~ mean*cond+ (1 + mean|sub), REML = F, data =  wb_dat)

summary(wbmod2acc_s)

# Perform the bootstrapping
boot_intercept <- bootMer(wbmod2acc_s, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for wbmod2acc_s | sub':", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(wbmod2acc_s, extract_random_slope, nsim = 1000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope Model 2 wb slope 'mean act | sub':", slope_ci, "\n")




#### Model 3: Level-1 and Level-2 Model Random Intercept Only ####
wbmod3acc = lmer(acc ~ mean*cond + handedness + female  + (1 |sub), REML = F, data =  wb_dat)

summary(wbmod3acc)

# Perform the bootstrapping
boot_intercept <- bootMer(wbmod3acc, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for wbmod1acc | sub':", intercept_ci, "\n")



lrtest(wbmod1acc,wbmod2acc)
lrtest(wbmod2acc,wbmod2acc_s)
lrtest(wbmod2acc,wbmod3acc)
# wbmod 3 is better








####### LOADING NEW DATA FOR OTHER CHECKS #######
# Importing data for RQ1 (DV= dm/pm, L1= cond, L2= sub)
dat = read.csv("C:/Users/anama/OneDrive/Documents/DMPM_Paper/Gram_ROIsep.csv")

# Transforming variables to numeric/continuous
dat$handedness<- as.numeric(as.character(dat$handedness))
dat$kbit2<- as.numeric(as.character(dat$kbit2))

# Interpolating missing values linearly (only 2 per variable)
dat$handedness <- na_interpolation(dat,option = 'linear')$handedness
dat$kbit2 <- na_interpolation(dat,option = 'linear')$kbit2
dat$rt <- na_interpolation(dat,option = 'linear')$rt

# Converting accuracy into a percent value
#dat$acc<- dat$acc*100

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
init = mice(dat[10:13], maxit = 0) 
meth = init$method
predM = init$predictorMatrix

# Setting variables to not be considered as a predictors during imputation
# In this case, we don't want timestamps to be factored into the calculation for imputing EEG data
predM[1:9] <- 0

# Setting method for imputing (using "norm" for normal imputation)
cols_to_impute <- names(dat)[10:13]

# Use lapply to set the imputation method for each column
meth[cols_to_impute] <- lapply(meth[cols_to_impute], function(x) "norm")

set.seed(103)
imputed = mice(dat, method = meth, predictorMatrix = predM, m = 5)

dat_imp <- complete(imputed)
dat_imp$cond<-as.factor(dat_imp$cond)
dat_imp$cond <- factor(dat_imp$cond, levels = c("0", "1", "2"), labels = c("No Error", "Easy Error", "Hard Error"))

dat_imp$acc <- asin(sqrt(dat_imp$acc))




########################
##### CELF-5 Check #####
########################

## RQ1 ###########################################
#### Model 1 DM: Unconditional Model ####
mod1dmL = lmer(dmL ~ 1 + celf_sts + (1|sub), REML = F, data = dat_imp)

summary(mod1dmL)

dmLicc <- quantile(bootMer(mod1dmL, calc.icc, nsim = 1000)$t, c(0.025, 0.975))
cat("\nModel 1 DM CELF - dmL ICC:", dmLicc, "\n")

# Perform the bootstrapping
boot_intercept <- bootMer(mod1dmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for wbmod1acc | sub':", intercept_ci, "\n")




mod1dmR = lmer(dmR ~ 1 + celf_sts + (1|sub), REML = F, data = dat_imp)

summary(mod1dmR)

dmRicc <- quantile(bootMer(mod1dmR, calc.icc, nsim = 1000)$t, c(0.025, 0.975))
cat("\nModel 1 DM - dmR ICC:", dmRicc, "\n")

# Perform the bootstrapping
boot_intercept <- bootMer(mod1dmR, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for wbmod1acc | sub':", intercept_ci, "\n")

tab_model(mod1dmL, mod1dmR)


#### Model 2 DM: Level-1 Model ####
mod2dmL = lmer(dmL ~ cond  + celf_sts + (1 + cond|sub), REML = F, data = dat_imp)

summary(mod2dmL)

mod2dmR = lmer(dmR ~ cond  + celf_sts + (1|sub), REML = F, data = dat_imp)

summary(mod2dmR)

lrtest(mod1dmL, mod2dmL) # mod 2 is better for both
lrtest(mod1dmR, mod2dmR)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2dmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 dmL '1 | sub':", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(mod2dmL, extract_random_slope, nsim = 1000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope Model 2 dmL 'condition | sub':", slope_ci, "\n")




#### Model 1 PM: Unconditional Model ####
mod1pmL = lmer(pmL ~ 1 + celf_sts + (1|sub), REML = F, data = dat_imp)

summary(mod1pmL)

pmLicc <- quantile(bootMer(mod1pmL, calc.icc, nsim = 1000)$t, c(0.025, 0.975))
cat("\nModel 1 PM L ICC:", pmLicc, "\n")



mod1pmR = lmer(pmR ~ 1 + celf_sts + (1|sub), REML = F, data = dat_imp)

summary(mod1pmR)

pmRicc <- quantile(bootMer(mod1pmR, calc.icc, nsim = 1000)$t, c(0.025, 0.975))
cat("\nModel 1 PM R ICC:", pmRicc, "\n")


#### Model 2 PM: Level-1 Model Random Intercept Only ####

mod2pmL = lmer(pmL ~ celf_sts + cond  + (1|sub), REML = F, data = dat_imp)

summary(mod2pmL)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2pmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI", intercept_ci, "\n")


mod2pmR = lmer(pmR ~ celf_sts + cond +  (1|sub), REML = F, data = dat_imp)

summary(mod2pmR)

lrtest(mod1pmL, mod2pmL)
lrtest(mod1pmR, mod2pmR) 

# Perform the bootstrapping
boot_intercept <- bootMer(mod2pmR, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 pmR '1 | sub':", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(mod2pmL, extract_random_slope, nsim = 1000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope Model 2 pmL 'condition | sub':", slope_ci, "\n")








### RQ2 ###########################################

#### Model 1: Unconditional Model ####
mod1acc = lmer(acc~ 1 + celf_sts + (1|sub), REML = F, data = dat_imp)

summary(mod1acc)

mod1accicc <- quantile(bootMer(mod1acc, calc.icc, nsim = 1000)$t, c(0.025, 0.975), na.rm=TRUE)
cat("\nModel 1 PM R ICC:", mod1accicc, "\n")

tab_model(mod1acc)

#### Model 2: Level-1 Model Linear ####

mod2accdmL = lmer(acc ~ celf_sts + dmL*cond  + (1|sub), REML = F, data = dat_imp)

summary(mod2accdmL)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accdmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI", intercept_ci, "\n")


mod2accdmR = lmer(acc ~ celf_sts + dmR*cond  + (1 |sub), REML = F, data = dat_imp)

summary(mod2accdmR)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accdmR, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for DM R | sub':", intercept_ci, "\n")




mod2accpmL = lmer(acc ~ celf_sts + pmL*cond  + (1 |sub), REML = F, data = dat_imp)

summary(mod2accpmL)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accpmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI", intercept_ci, "\n")



mod2accpmR = lmer(acc ~ celf_sts + pmR*cond  +(1 |sub), REML = F, data = dat_imp)

summary(mod2accpmR)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accpmR, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI", intercept_ci, "\n")

lrtest(mod1acc,mod2accdmL)
lrtest(mod1acc,mod2accdmR)
lrtest(mod1acc,mod2accpmL)
lrtest(mod1acc,mod2accpmR)






########################
##### KBIT2 Check #####
########################

## RQ1 ###########################################
#### Model 1 DM: Unconditional Model ####
mod1dmL = lmer(dmL ~ 1 + kbit2 + (1|sub), REML = F, data = dat_imp)

summary(mod1dmL)

dmLicc <- quantile(bootMer(mod1dmL, calc.icc, nsim = 1000)$t, c(0.025, 0.975))
cat("\nModel 1 DM CELF - dmL ICC:", dmLicc, "\n")

# Perform the bootstrapping
boot_intercept <- bootMer(mod1dmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for wbmod1acc | sub':", intercept_ci, "\n")




mod1dmR = lmer(dmR ~ 1 + kbit2 + (1|sub), REML = F, data = dat_imp)

summary(mod1dmR)

dmRicc <- quantile(bootMer(mod1dmR, calc.icc, nsim = 1000)$t, c(0.025, 0.975))
cat("\nModel 1 DM - dmR ICC:", dmRicc, "\n")

# Perform the bootstrapping
boot_intercept <- bootMer(mod1dmR, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for wbmod1acc | sub':", intercept_ci, "\n")

tab_model(mod1dmL, mod1dmR)


#### Model 2 DM: Level-1 Model ####
mod2dmL = lmer(dmL ~ cond  + kbit2 + (1 + cond|sub), REML = F, data = dat_imp)

summary(mod2dmL)

mod2dmR = lmer(dmR ~ cond  + kbit2 + (1|sub), REML = F, data = dat_imp)

summary(mod2dmR)

lrtest(mod1dmL, mod2dmL) # mod 2 is better for both
lrtest(mod1dmR, mod2dmR)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2dmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 dmL '1 | sub':", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(mod2dmL, extract_random_slope, nsim = 1000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope Model 2 dmL 'condition | sub':", slope_ci, "\n")




#### Model 1 PM: Unconditional Model ####
mod1pmL = lmer(pmL ~ 1 + kbit2 + (1|sub), REML = F, data = dat_imp)

summary(mod1pmL)

pmLicc <- quantile(bootMer(mod1pmL, calc.icc, nsim = 1000)$t, c(0.025, 0.975))
cat("\nModel 1 PM L ICC:", pmLicc, "\n")



mod1pmR = lmer(pmR ~ 1 + kbit2 + (1|sub), REML = F, data = dat_imp)

summary(mod1pmR)

pmRicc <- quantile(bootMer(mod1pmR, calc.icc, nsim = 1000)$t, c(0.025, 0.975))
cat("\nModel 1 PM R ICC:", pmRicc, "\n")


#### Model 2 PM: Level-1 Model Random Intercept Only ####

mod2pmL = lmer(pmL ~ kbit2 + cond  + (1|sub), REML = F, data = dat_imp)

summary(mod2pmL)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2pmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI", intercept_ci, "\n")


mod2pmR = lmer(pmR ~ kbit2 + cond +  (1|sub), REML = F, data = dat_imp)

summary(mod2pmR)

lrtest(mod1pmL, mod2pmL)
lrtest(mod1pmR, mod2pmR) 

# Perform the bootstrapping
boot_intercept <- bootMer(mod2pmR, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for Model 2 pmR '1 | sub':", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(mod2pmL, extract_random_slope, nsim = 1000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope Model 2 pmL 'condition | sub':", slope_ci, "\n")








### RQ2 ###########################################

#### Model 1: Unconditional Model ####
mod1acc = lmer(acc~ 1 + kbit2 + (1|sub), REML = F, data = dat_imp)

summary(mod1acc)

mod1accicc <- quantile(bootMer(mod1acc, calc.icc, nsim = 1000)$t, c(0.025, 0.975), na.rm=TRUE)
cat("\nModel 1 PM R ICC:", mod1accicc, "\n")

tab_model(mod1acc)

#### Model 2: Level-1 Model Linear ####

mod2accdmL = lmer(acc ~ kbit2 + dmL*cond  + (1|sub), REML = F, data = dat_imp)

summary(mod2accdmL)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accdmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI", intercept_ci, "\n")


mod2accdmR = lmer(acc ~ kbit2 + dmR*cond  + (1 |sub), REML = F, data = dat_imp)

summary(mod2accdmR)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accdmR, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI for DM R | sub':", intercept_ci, "\n")




mod2accpmL = lmer(acc ~ kbit2 + pmL*cond  + (1 |sub), REML = F, data = dat_imp)

summary(mod2accpmL)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accpmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI", intercept_ci, "\n")



mod2accpmR = lmer(acc ~ kbit2 + pmR*cond  +(1 |sub), REML = F, data = dat_imp)

summary(mod2accpmR)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accpmR, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI", intercept_ci, "\n")

lrtest(mod1acc,mod2accdmL)
lrtest(mod1acc,mod2accdmR)
lrtest(mod1acc,mod2accpmL)
lrtest(mod1acc,mod2accpmR)








######################## MRI QC Correlations #####################################################

# Load necessary packages
library("ggpubr")
library("ggplot2")
library(sjPlot)
library(car)
library(corrplot)


# Importing data for RQ1 (DV= dm/pm, L1= cond, L2= sub)
qc_dat = read.csv("C:/Users/anama/OneDrive/Documents/DMPM_Paper/qc_Gram.csv")

# Transforming variables to numeric/continuous
qc_dat$handedness<- as.numeric(as.character(qc_dat$handedness))
qc_dat$kbit2<- as.numeric(as.character(qc_dat$kbit2))

# Interpolating missing values linearly (only 2 per variable)
qc_dat$handedness <- na_interpolation(qc_dat,option = 'linear')$handedness
qc_dat$kbit2 <- na_interpolation(qc_dat,option = 'linear')$kbit2
qc_dat$rt <- na_interpolation(qc_dat,option = 'linear')$rt


# Find average RT per condition
aggregate(qc_dat$rt, list(qc_dat$cond), FUN=mean) 

# Centering reaction time by cond
qc_dat$rtc <- ifelse(qc_dat$cond == 0, qc_dat$rt - 2.32683, 
                  ifelse(qc_dat$cond == 1, qc_dat$rt - 2.256074, 
                         ifelse(qc_dat$cond == 2, qc_dat$rt - 3.638021, qc_dat$rt)))

# Turn from character to numeric
# Loop through each column
for (col in names(qc_dat)) {
  qc_dat[[col]] <- as.numeric(qc_dat[[col]])
}


library(data.table)
library(mice)
# Here is where you would input your "save qc_data set" command of choice. I wasn't 
# sure if .txt or .csv would be moreinit = mice(qc_dat[2:nrow(qc_dat)], maxit = 0) 
init = mice(qc_dat[10:18], maxit = 0) 
meth = init$method
predM = init$predictorMatrix

# Setting variables to not be considered as a predictors during imputation
# In this case, we don't want timestamps to be factored into the calculation for imputing EEG qc_data
predM[1:9] <- 0

# Setting method for imputing (using "norm" for normal imputation)
cols_to_impute <- names(qc_dat)[10:18]

# Use lapply to set the imputation method for each column
meth[cols_to_impute] <- lapply(meth[cols_to_impute], function(x) "norm")

set.seed(103)
imputed = mice(qc_dat, method = meth, predictorMatrix = predM, m = 5)

qc_dat_imp <- complete(imputed)

qc_dat_imp$acc <- asin(sqrt(qc_dat_imp$acc))

qc<- qc_dat_imp[,c(1,6:7,9:10, 13:20)]

# List of variables to standardize
vars_to_standardize <- c( "num_repaired","efc","gcor","aor","aqi","snr","tsnr","fber")

# Loop to standardize variables
for (var in vars_to_standardize) {
  qc_dat_imp[[var]] <- scale(qc_dat_imp[[var]])
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
  as.numeric(VarCorr(model)$sub[2, 2])
}

# Define a function to extract the random intercept term variance component
extract_random_intercept <- function(model) {
  as.numeric(VarCorr(model)$sub[1, 1])
}






correlation_matrix <- corr.test(qc)$r

# Extract p-values from the corr.test result
p_values <- corr.test(qc)$p

se <- corr.test(qc)$se

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
  text.decorate = "none"  # Do not bold the text
)





# Set the output directory and file name
output_dir <- "C:/Users/anama/OneDrive/Documents/DMPM_Paper"
output_file <- file.path(output_dir, "QC_output_st.txt")

# Redirect output to a text file
sink(output_file)

# PML VS. QC

# Function to perform Pearson correlation test and format the results
perform_correlation <- function(var1, var2, data) {
  test <- cor.test(data[[var1]], data[[var2]], method = "pearson")
  result <- sprintf("#pearson (%s vs %s) -- r = %.7f, df = %d, p-value = %.6f\n", 
                    var1, var2, test$estimate, test$parameter, test$p.value)
  return(result)
}

# List of variable pairs for correlation tests
tests <- list(
  list(var1 = "pmL", var2 = "num_repaired"),
  list(var1 = "pmL", var2 = "aor"),
  list(var1 = "pmL", var2 = "snr"),
  list(var1 = "pmL", var2 = "fber"),
  list(var1 = "dmL", var2 = "num_repaired"),
  list(var1 = "dmL", var2 = "aor"),
  list(var1 = "dmL", var2 = "snr"),
  list(var1 = "dmL", var2 = "fber"),
  list(var1 = "cond", var2 = "num_repaired"),
  list(var1 = "cond", var2 = "aor"),
  list(var1 = "cond", var2 = "snr"),
  list(var1 = "cond", var2 = "fber"),
  list(var1 = "sub", var2 = "num_repaired"),
  list(var1 = "sub", var2 = "aor"),
  list(var1 = "sub", var2 = "snr"),
  list(var1 = "sub", var2 = "fber")
  
)

# Perform the tests and collect results
results <- lapply(tests, function(test) {
  perform_correlation(test$var1, test$var2, qc_dat_imp)
})

cat("\nCorrelations Between QC Indicators and DV/IV\n")

# Write the results to a text file
writeLines(unlist(results), output_file)

qc_dat_imp$cond<-as.factor(qc_dat_imp$cond)
qc_dat_imp$cond <- factor(qc_dat_imp$cond, levels = c("0", "1", "2"), labels = c("No Error", "Easy Error", "Hard Error"))

### Num_Repaired ##########################
### RQ1 ###
cat("\nNumber of Volumes Repaired: Lower is Better\n")
cat("\nModel 1 DM: Unconditional Model num_repaired\n")
mod1dmL <- lmer(dmL ~ 1 + num_repaired + (1 | sub), REML = FALSE, data = qc_dat_imp)
print(summary(mod1dmL))

# Perform the bootstrapping
boot_intercept <- bootMer(mod1dmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI num_repaired:", intercept_ci, "\n")



# Model 2 PM: Level-1 Model
cat("\nModel 2 PML: Level-1 Model num_repaired\n")
mod2pmL <- lmer(pmL ~ cond + num_repaired + (1 | sub), REML = FALSE, data = qc_dat_imp)
print(summary(mod2pmL))

# Perform the bootstrapping
boot_intercept <- bootMer(mod2pmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI num_repaired:", intercept_ci, "\n")


#### RQ2 ####

mod2accdmL <- lmer(acc ~ dmL * cond + num_repaired + (1 + dmL | sub), REML = FALSE, data = qc_dat_imp)
mod2accdmL_summary <- summary(mod2accdmL)

cat("\n Model 2: Level-1 Model \n")
print(mod2accdmL_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accdmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI num_repaired:", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(mod2accdmL, extract_random_slope, nsim = 1000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope CI num_repaired:", slope_ci, "\n")




mod2accpmL <- lmer(acc ~ pmL * cond + num_repaired + (1 | sub), REML = FALSE, data = qc_dat_imp)
mod2accpmL_summary <- summary(mod2accpmL)

cat("\n Model 2: Level-1 Model\n")
print(mod2accpmL_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accpmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI num_repaired", intercept_ci, "\n")



### Artifact to Original Ratio (AOR) ##########################
### RQ1 ###
cat("\nArtifact to Original Ratio (AOR): Lower is Better\n")
cat("\nModel 1 DM: Unconditional Model aor\n")
mod1dmL <- lmer(dmL ~ 1 + aor + (1 | sub), REML = FALSE, data = qc_dat_imp)
print(summary(mod1dmL))

# Perform the bootstrapping
boot_intercept <- bootMer(mod1dmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI aor:", intercept_ci, "\n")

mod1dmLicc <- quantile(bootMer(mod1pmL, calc.icc, nsim = 1000)$t, c(0.025, 0.975), na.rm=TRUE)
cat("\nModel ICC:", mod1dmLicc, "\n")



# Model 2 PM: Level-1 Model
cat("\nModel 2 PML: Level-1 Model aor\n")
mod2pmL <- lmer(pmL ~ cond + aor + (1 | sub), REML = FALSE, data = qc_dat_imp)
print(summary(mod2pmL))

# Perform the bootstrapping
boot_intercept <- bootMer(mod2pmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI aor:", intercept_ci, "\n")


#### RQ2 ####

mod2accdmL <- lmer(acc ~ dmL * cond + aor + (1 + dmL | sub), REML = FALSE, data = qc_dat_imp)
mod2accdmL_summary <- summary(mod2accdmL)

cat("\n Model 2: Level-1 Model \n")
print(mod2accdmL_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accdmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI aor:", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(mod2accdmL, extract_random_slope, nsim = 1000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope CI aor:", slope_ci, "\n")




mod2accpmL <- lmer(acc ~ pmL * cond + aor + (1 | sub), REML = FALSE, data = qc_dat_imp)
mod2accpmL_summary <- summary(mod2accpmL)

cat("\n Model 2: Level-1 Model\n")
print(mod2accpmL_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accpmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI aor", intercept_ci, "\n")



### Signal to Noise Ratio (SNR): Higher is Better ##########################
### RQ1 ###
cat("\nSignal to Noise Ratio (SNR): Higher is Better\n")
cat("\nModel 1 DM: Unconditional Model snr\n")
mod1dmL <- lmer(dmL ~ 1 + snr + (1 | sub), REML = FALSE, data = qc_dat_imp)
print(summary(mod1dmL))

# Perform the bootstrapping
boot_intercept <- bootMer(mod1dmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI snr:", intercept_ci, "\n")



# Model 2 PM: Level-1 Model
cat("\nModel 2 PML: Level-1 Model snr\n")
mod2pmL <- lmer(pmL ~ cond + snr + (1 | sub), REML = FALSE, data = qc_dat_imp)
print(summary(mod2pmL))

# Perform the bootstrapping
boot_intercept <- bootMer(mod2pmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI snr:", intercept_ci, "\n")


#### RQ2 ####

mod2accdmL <- lmer(acc ~ dmL * cond + snr + (1 + dmL | sub), REML = FALSE, data = qc_dat_imp)
mod2accdmL_summary <- summary(mod2accdmL)

cat("\n Model 2: Level-1 Model \n")
print(mod2accdmL_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accdmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI snr:", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(mod2accdmL, extract_random_slope, nsim = 1000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope CI snr:", slope_ci, "\n")




mod2accpmL <- lmer(acc ~ pmL * cond + snr + (1 | sub), REML = FALSE, data = qc_dat_imp)
mod2accpmL_summary <- summary(mod2accpmL)

cat("\n Model 2: Level-1 Model\n")
print(mod2accpmL_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accpmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI snr", intercept_ci, "\n")



### Foreground-to-Background Energy Ratio (FBER) ##########################
### RQ1 ###
cat("\nForeground-to-Background Energy Ratio (FBER): Higher is Better\n")
cat("\nModel 1 DM: Unconditional Model fber\n")
mod1dmL <- lmer(dmL ~ 1 + fber + (1 | sub), REML = FALSE, data = qc_dat_imp)
print(summary(mod1dmL))

# Perform the bootstrapping
boot_intercept <- bootMer(mod1dmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI fber:", intercept_ci, "\n")



# Model 2 PM: Level-1 Model
cat("\nModel 2 PML: Level-1 Model fber\n")
mod2pmL <- lmer(pmL ~ cond + fber + (1 | sub), REML = FALSE, data = qc_dat_imp)
print(summary(mod2pmL))

# Perform the bootstrapping
boot_intercept <- bootMer(mod2pmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI fber:", intercept_ci, "\n")


#### RQ2 ####

mod2accdmL <- lmer(acc ~ dmL * cond + fber + (1 + dmL | sub), REML = FALSE, data = qc_dat_imp)
mod2accdmL_summary <- summary(mod2accdmL)

cat("\n Model 2: Level-1 Model \n")
print(mod2accdmL_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accdmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI fber:", intercept_ci, "\n")

# Perform the bootstrapping
boot_slope <- bootMer(mod2accdmL, extract_random_slope, nsim = 1000)

# Calculate the quantiles for the random slope term
slope_ci <- quantile(boot_slope$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random slope term
cat("Random slope CI fber:", slope_ci, "\n")


mod2accpmL <- lmer(acc ~ pmL * cond + fber + (1 | sub), REML = FALSE, data = qc_dat_imp)
mod2accpmL_summary <- summary(mod2accpmL)

cat("\n Model 2: Level-1 Model\n")
print(mod2accpmL_summary)

# Perform the bootstrapping
boot_intercept <- bootMer(mod2accpmL, extract_random_intercept, nsim = 1000)

# Calculate the quantiles for the random intercept term
intercept_ci <- quantile(boot_intercept$t, c(0.025, 0.975), na.rm = TRUE)

# Print the confidence interval for the random intercept term
cat("Random Intercept Term CI fber", intercept_ci, "\n")



# Stop redirecting output to the text file
sink()
