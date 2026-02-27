#####################################################################
# Tchetgen Tchetgen et al. double negative control calibration method
#####################################################################

# Purpose:
# This script implements the regression-based double negative control method proposed by Tchetgen Tchetgen et al. (DOI: 10.1214/23-STS911) and applied 
# by Li et al. (DOI: 10.1097/EDE.0000000000001884).
# These methods are adapted to align with a clone-censor weight approach.
# In our application, having a wellness visit is the negative control exposure, and hip fracture is the negative control outcome. 

# Notes:
# Due to the lack of guidance on implementing double negative controls within 
# the clone-censor weight framework, results of this analysis should be interpreted cautiously. 

# Set the directory
setwd("P:/your data folder path")

## Load the required libraries
library(survival)
library(survminer)
library (tidyverse)
library(stringr)
library(ggplot2)
library(survminer)
library(ggsurvfit)
library (dplyr)
library(KMsurv)
library(tidyr)
library(furrr)
library(future)
library(here)
library(lubridate)
library(data.table)
library(haven)
library(here())

# Use the here package to create a reproducible file path.
# It constructs a path pointing to the R path:
here('src', 'setup.R')

# Prints text to console.
cat('loading canon packages....')

# Package List Definition
pkg_toload <- c('tidyverse', 
                'lubridate', 'here', 'knitr', 'quarto',
                'survival', 'future', 'progressr', 
                'ggpubr', 'survminer', 'furrr',
                'data.table', 'parglm')

hold_del <- sapply(pkg_toload, require, 
                   warn.conflicts=F, quietly=T,
                   character.only=T)
here()

cat('Done', '\n')

# ------------------------------------------------------------
# Step 1: Define negative control exposures (wellness visit) as 
# as new treatment strategies 
# ------------------------------------------------------------

# Read in panel dataset in which residents (uncloned) are followed until a wellness visit or censoring
library(haven)
expanded_wellness_0 <- read_sas("P:/your folder path/expanded_wellness_cov_tv.sas7bdat",  NULL)


# Increase global variable size for bootstrapping procedure
options(future.globals.maxSize = +Inf)

# Set analysis options for bootstrapping procedure and store in list and saved
d_output = list(runtime = Sys.time(),
                #params = l_tte_params,
                runplan = list(boots = 1,
                               workers = 1)
                )

# Calculate predicted probability of a wellness visit by interval as a function of
# baseline and time-varying predictors

d_glm_wt1_wellness = glm(outcome_wellness ~ poly(t_intrv, 2, raw=T) + as.factor(BL_SEX) + poly(BL_AGE, 2) + as.factor(BL_ADL_CAT) + BL_ANEMIA + BL_ANALGESICS
               + BL_ANTICOAGULANTS + BL_ANTICONVULSANTS + BL_ANTIDEPRESSANTS + BL_ANTIPSYCHOTICS + BL_ANTIVIRALS + BL_ARTHRITIS +
                 BL_AFIB + BL_BENZODIAZEPINES + BL_CANCER + as.factor(BL_M3_CFS) + BL_CAD + BL_COVIDVAX + BL_DVT + BL_DIABETES + 
                 as.factor(BL_ED_VISITS) + BL_FLUVAX + BL_GABAPENTANOIDS + as.factor(BL_GAGNE_CAT) + BL_HD + BL_HERPESZOSTER + 
                 BL_HCHOL + as.factor(BL_HOSP_VISITS) + BL_HTN + BL_INSULIN + BL_MENTALHEALTH + as.factor(BL_MONTH) + BL_OSTEOFRAC +
                 BL_PARKINSONS + BL_PVD + BL_PNEUMOVAX + BL_NONHIPFX + as.factor(BL_RACE) + BL_RENALF + BL_STROKE + 
                 as.factor(BL_REGION) + as.factor (BL_YEAR) + BL_ZOSTAVAX + BL_PARTBCLAIM + 
                 TV_ANTICOAGULANTS + TV_ANTIDEPRESSANTS + TV_ANTIPSYCHOTICS + TV_ANTIVIRALS + TV_CANCER + TV_COVIDVAX + TV_DVT + TV_DIABETES + 
                 TV_ED_VISITS + TV_FLUVAX + TV_HERPESZOSTER + TV_HOSP_VISITS + TV_INSULIN + TV_DISCHARGE + TV_PNEUMOVAX + TV_RENALF + TV_STROKE 
                 + TV_ZOSTAVAX + TV_SHINGRIX, 
               data=expanded_wellness_0, family=binomial())

# Store coefficients 

summary_d_glm_wt1_wellness<- summary(d_glm_wt1_wellness)
coef_df_wellness <- as.data.frame((summary_d_glm_wt1_wellness$coefficients))

write.csv(coef_df_wellness,file='wellness_vis_model.csv', row.names=TRUE)

expanded_wellness_0$pr_treat_wellness = d_glm_wt1_wellness$fitted.values  
setDT(expanded_wellness_0)
expanded_wellness_0[, cumpr_notreat_wellness := cumprod(1-pr_treat_wellness), by = .(clientid)]

# ------------------------------------------------------------
# Step 2: Define negative control outcome (hip fracture)
# ------------------------------------------------------------

# Import hip fracture outcome panel dataset containing both sets of clones

hipfrac_outcome_1 <- read_sas("P:/your folder path/expanded_hipfrac_cov_tv.sas7bdat", NULL)

# Join wellness visit probabilities to outcome data set
hipfrac_outcome1_join = left_join(hipfrac_outcome_1,
                                   select(expanded_wellness_0, clientid, t_intrv, cumpr_notreat_wellness),
                                   by=c('clientid', 't_intrv'))

# Calculate wellness visit inverse probability weights
setDT(hipfrac_outcome1_join)

hipfrac_outcome1_join[, ipw_wellness := fcase(
  Treat==0, 1 / cumpr_notreat_wellness,
  Treat==1 & t_intrv < 12, 1, 
  Treat==1 & t_intrv == 12 & t_treat < 12, 1, 
  Treat==1 & t_intrv==12 & t_treat==12, 1 / (1-cumpr_notreat_wellness), 
  Treat==1 & t_intrv==12 & t_treat>12, 0, 
  Treat==1 & t_intrv > 12, 1, 
  Treat==1 & t_intrv == 12 & is.na(t_treat), 0, 
  Treat==1 & t_intrv > 12 & is.na(t_treat), 1 
)]

# For assigned clones (Treat==1), carry the weight assigned at time 12 forward to end of follow-up

hipfrac_outcome1_join[Treat==1, ipw1_wellness := cumprod(ipw_wellness), by=list(clientid, Treat)]
hipfrac_outcome1_join[Treat==0, ipw1_wellness := ipw_wellness, by=list(clientid, Treat)]

# Check IPW summary statistics

summary(hipfrac_outcome1_join$ipw1_wellness)

hipfrac_outcome1_join %>%
  filter(ipw1_wellness != 0) %>%
  group_by(Treat) %>%
  summarise(summary_stats = list(as.data.frame(as.list(summary(ipw1_wellness))))) %>%
  unnest(summary_stats)

quantile(hipfrac_outcome1_join$ipw1_wellness[hipfrac_outcome1_join$ipw1_wellness != 0],
         probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.995, 0.997, 0.998),
         na.rm = TRUE)

hipfrac_outcome1_join %>%
  filter(ipw1_wellness != 0 & !is.na(ipw1_wellness)) %>%
  group_by(Treat) %>%
  summarise(quantiles = list(quantile(ipw1_wellness, probs = c(0.01, 0.05, 0.25, 0.5, 
                                                      0.75, 0.9, 0.95, 0.99, 
                                                      0.995, 0.997, 0.998), 
                                      na.rm = TRUE))) %>%
  unnest_wider(quantiles)

# Run hip fracture outcome model 

d_glm_pe2_hipfrac = glm(outcome_hip_frac ~ poly(t_intrv, 2, raw=T)*Treat,data=hipfrac_outcome1_join,
                family= binomial(), weights = ipw1_wellness)


# Retrieve coefficients 

summary_d_glm_pe2_hipfrac<- summary(d_glm_pe2_hipfrac)
coef_df_hipfrac <- as.data.frame((summary_d_glm_pe2_hipfrac$coefficients))

write.csv(coef_df_hipfrac,file='hipfrac_model.csv', row.names=TRUE)

# Retrieve predicted probabilities of hip fracture 
hipfrac_outcome1_join$hipfrac_hat = d_glm_pe2_hipfrac$fitted.values


# Calculate cumulative risk of hip fracture by exposure group and interval

hipfrac_outcome1_join[, pr_surv_hipfrac := cumprod(1 - hipfrac_hat), by = list(clientid, Treat)]

d_res1 = hipfrac_outcome1_join %>%
  group_by(Treat, t_intrv) %>%
  summarize(pr_ev_hipfrac = mean(1-pr_surv_hipfrac), .groups = 'drop') %>%
  ungroup %>%
  pivot_wider(., id_cols =c('t_intrv'),
              names_from = Treat,
              names_prefix = 'pr_ev_hipfrac',
              values_from = pr_ev_hipfrac
  ) %>%
  mutate(cid = pr_ev_hipfrac1 - pr_ev_hipfrac0,
         
         cir = pr_ev_hipfrac1 / pr_ev_hipfrac0)

write.csv(d_res1,file='d_res1_hipfrac.csv',
          row.names=FALSE)

# -----------------------------------------------------------------------
# Step 3: Perform IPCW analysis of vaccine exposure and incident dementia
# adjusted for hip fracture predicted probability 
# -----------------------------------------------------------------------

# Read in panel dataset in which residents (uncloned) are followed until RZV or censoring

library(haven)
expanded_vaccine_0 <- read_sas("P:/your folder path/expanded_vaccine_cov.sas7bdat",  NULL)

# Calculate predicted probabilities of vaccination

d_glm_wt1 = glm(outcome_vaccine ~ poly(t_intrv, 2, raw=T) + as.factor(BL_SEX) + poly(BL_AGE, 2) + as.factor(BL_ADL_CAT) + BL_ANEMIA + BL_ANALGESICS
                + BL_ANTICOAGULANTS + BL_ANTICONVULSANTS + BL_ANTIDEPRESSANTS + BL_ANTIPSYCHOTICS + BL_ANTIVIRALS + BL_ARTHRITIS +
                  BL_AFIB + BL_BENZODIAZEPINES + BL_CANCER + as.factor(BL_M3_CFS) + BL_CAD + BL_COVIDVAX + BL_DVT + BL_DIABETES + 
                  as.factor(BL_ED_VISITS) + BL_FLUVAX + BL_GABAPENTANOIDS + as.factor(BL_GAGNE_CAT) + BL_HD + BL_HERPESZOSTER + 
                  BL_HCHOL + as.factor(BL_HOSP_VISITS) + BL_HTN + BL_INSULIN + BL_MENTALHEALTH + as.factor(BL_MONTH) + BL_OSTEOFRAC +
                  BL_PARKINSONS + BL_PVD + BL_PNEUMOVAX + BL_NONHIPFX + as.factor(BL_RACE) + BL_RENALF + BL_STROKE + 
                  as.factor(BL_REGION) + as.factor (BL_YEAR) + BL_ZOSTAVAX + BL_PARTBCLAIM + 
                  TV_ANTICOAGULANTS + TV_ANTIDEPRESSANTS + TV_ANTIPSYCHOTICS + TV_ANTIVIRALS + TV_CANCER + TV_COVIDVAX + TV_DVT + TV_DIABETES + 
                  TV_ED_VISITS + TV_FLUVAX + TV_HERPESZOSTER + TV_HOSP_VISITS + TV_INSULIN + TV_DISCHARGE + TV_PNEUMOVAX + TV_RENALF + TV_STROKE + TV_ZOSTAVAX, 
                data=expanded_vaccine_0, family=binomial())

# Retrieve coefficients 

summary_d_glm_wt1<- summary(d_glm_wt1)
coef_df_vaccine <- as.data.frame((summary_d_glm_wt1$coefficients))

write.csv(coef_df_vaccine,file='vaccine_model.csv', row.names=TRUE)

# Get cumulative vaccine probabilities 

expanded_vaccine_0$pr_treat = d_glm_wt1$fitted.values  
setDT(expanded_vaccine_0)
expanded_vaccine_0[, cumpr_notreat := cumprod(1-pr_treat), by = .(clientid)]

# Import the dementia outcome panel dataset

Dementia_outcome_1 <- read_sas("P:/hzvnh/shared/Data/Dementia Outcome 1/dementia_outcome1_final.sas7bdat", NULL)

# Join vaccine probabilities to outcome data set 
dementia_outcome1_join = left_join(Dementia_outcome_1,
                                   select(expanded_vaccine_0, clientid, t_intrv, cumpr_notreat),
                                   by=c('clientid', 't_intrv'))

# Calculate IPWs

setDT(dementia_outcome1_join)

dementia_outcome1_join[, ipw := fcase(
  Treat==0, 1 / cumpr_notreat,
  Treat==1 & t_intrv < 12, 1, 
  Treat==1 & t_intrv == 12 & t_treat < 12, 1, 
  Treat==1 & t_intrv==12 & t_treat==12, 1 / (1-cumpr_notreat), 
  Treat==1 & t_intrv==12 & t_treat>12, 0, 
  Treat==1 & t_intrv > 12, 1, 
  Treat==1 & t_intrv == 12 & is.na(t_treat), 0, 
  Treat==1 & t_intrv > 12 & is.na(t_treat), 1 
)]

dementia_outcome1_join[Treat==1, ipw1 := cumprod(ipw), by=list(clientid, Treat)]
dementia_outcome1_join[Treat==0, ipw1 := ipw, by=list(clientid, Treat)]

# Check IPW distribution

summary(dementia_outcome1_join$ipw1)

dementia_outcome1_join %>%
  filter(ipw1 != 0) %>%
  group_by(Treat) %>%
  summarise(summary_stats = list(as.data.frame(as.list(summary(ipw1))))) %>%
  unnest(summary_stats)

quantile(dementia_outcome1_join$ipw1[dementia_outcome1_join$ipw1 != 0],
         probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.995, 0.997, 0.998),
         na.rm = TRUE)

dementia_outcome1_join %>%
  filter(ipw1 != 0 & !is.na(ipw1)) %>%
  group_by(Treat) %>%
  summarise(quantiles = list(quantile(ipw1, probs = c( 0.5, 
                                                       0.75, 0.9, 0.95, 0.99, 
                                                       0.995, 0.997, 0.998), 
                                      na.rm = TRUE))) %>%
  unnest_wider(quantiles)


# Merge hip fracture predicted probabilities to dementia outcome dataset

dementia_hipfrac_join = left_join(dementia_outcome1_join,
                                   select(hipfrac_outcome1_join, clientid, t_intrv, Treat, hipfrac_hat),
                                   by=c('clientid', 't_intrv', 'Treat'))

# Check for missing hip fracture predicted probabilities

sum(is.na(dementia_hipfrac_join$hipfrac_hat))

miss_hipfrac_hat <- dementia_hipfrac_join %>% 
  filter(is.na(hipfrac_hat))


# Sort data and carry forward last hip fracture predicted probability available

dementia_hipfrac_join_sort <- dementia_hipfrac_join  %>% 
  arrange(clientid, Treat, t_intrv)

install.packages('zoo')
library(zoo)

setDT(dementia_hipfrac_join_sort)

dementia_hipfrac_join_sort[, hipfrac_hat := na.locf(hipfrac_hat), by=c('clientid','Treat')]

sum(is.na(dementia_hipfrac_join_sort$hipfrac_hat))


# Outcome model, dementia 

d_glm_pe_dem = glm(outcome_dementia ~ poly(t_intrv, 2, raw=T)*Treat + hipfrac_hat,data=dementia_hipfrac_join_sort,
                family= binomial(), weights = ipw1)

dementia_hipfrac_join_sort$pr_ev_dem = d_glm_pe_dem$fitted.values

# Retrieve coefficients 

summary_d_glm_pe_dem <- summary(d_glm_pe_dem)
coef_df_dem <- as.data.frame((summary_d_glm_pe_dem$coefficients))

write.csv(coef_df_dem,file='dementia_model.csv', row.names=TRUE)

# Calculative cumulative risk of dementia by interval and exposure group

dementia_hipfrac_join_sort[, pr_surv_dem := cumprod(1 - pr_ev_dem), by = list(clientid, Treat)]

d_res1_dem = dementia_hipfrac_join_sort %>%
  group_by(Treat, t_intrv) %>%
  summarize(pr_ev_dem = mean(1-pr_surv_dem), .groups = 'drop') %>%
  ungroup %>%
  pivot_wider(., id_cols =c('t_intrv'),
              names_from = Treat,
              names_prefix = 'pr_ev_dem',
              values_from = pr_ev_dem
  ) %>%
  mutate(cid = pr_ev_dem1 - pr_ev_dem0,
         
         cir = pr_ev_dem1 / pr_ev_dem0)

write.csv(d_res1_dem,file='d_res1_dementia.csv',
          row.names=FALSE)

# -----------------------------------------------------------------------
# Step 4: Implement Poisson percentile-based bootstrapping procedure
# -----------------------------------------------------------------------

# Increase global variable size for bootstrapping procedure
options(future.globals.maxSize = +Inf)

# Set analysis options for bootstrapping procedure and store in list and saved
d_output = list(runtime = Sys.time(),
                #params = l_tte_params,
                runplan = list(boots = 1,
                               workers = 1)
                )

# Define bootstrapping function

d_fun_getplrwt = function(dta_outc, dta_treat, ...) {
  
  # For Poisson bootstrap - cluster by person
  # Generate a frequency weight for each person
  # Wt ~ Poisson(Lambda=1)
  d_freqwt = distinct(dta_outc, clientid) %>%
    mutate(., freqwt = rpois(n(), 1L))
  
  dta_outc = left_join(dta_outc, d_freqwt, by='clientid')
  
  dta_2 = left_join(dta_treat, d_freqwt, by='clientid')
  
  d_glm_wt = glm(outcome_vaccine ~ poly(t_intrv, 2, raw=T) + as.factor(BL_SEX) + poly(BL_AGE, 2) + as.factor(BL_ADL_CAT) + BL_ANEMIA + BL_ANALGESICS
                    + BL_ANTICOAGULANTS + BL_ANTICONVULSANTS + BL_ANTIDEPRESSANTS + BL_ANTIPSYCHOTICS + BL_ANTIVIRALS + BL_ARTHRITIS +
                      BL_AFIB + BL_BENZODIAZEPINES + BL_CANCER + as.factor(BL_M3_CFS) + BL_CAD + BL_COVIDVAX + BL_DVT + BL_DIABETES + 
                      as.factor(BL_ED_VISITS) + BL_FLUVAX + BL_GABAPENTANOIDS + as.factor(BL_GAGNE_CAT) + BL_HD + BL_HERPESZOSTER + 
                      BL_HCHOL + as.factor(BL_HOSP_VISITS) + BL_HTN + BL_INSULIN + BL_MENTALHEALTH + as.factor(BL_MONTH) + BL_OSTEOFRAC +
                      BL_PARKINSONS + BL_PVD + BL_PNEUMOVAX + BL_NONHIPFX + as.factor(BL_RACE) + BL_RENALF + BL_STROKE + 
                      as.factor(BL_REGION) + as.factor (BL_YEAR) + BL_ZOSTAVAX + BL_PARTBCLAIM + 
                      TV_ANTICOAGULANTS + TV_ANTIDEPRESSANTS + TV_ANTIPSYCHOTICS + TV_ANTIVIRALS + TV_CANCER + TV_COVIDVAX + TV_DVT + TV_DIABETES + 
                      TV_ED_VISITS + TV_FLUVAX + TV_HERPESZOSTER + TV_HOSP_VISITS + TV_INSULIN + TV_DISCHARGE + TV_PNEUMOVAX + TV_RENALF + TV_STROKE + TV_ZOSTAVAX, 
                    data=dta_2, family=binomial(), weights=freqwt)
  
  dta_2$pr_treat = d_glm_wt$fitted.values
  setDT(dta_2)
  dta_2[, cumpr_notreat := cumprod(1-pr_treat), by = .(clientid)]
  dta_3 = left_join(dta_outc, select(dta_2, clientid, t_intrv, cumpr_notreat), by=join_by(clientid, t_intrv))    
  setDT(dta_3)
  dta_3[, ipw := fcase(
    Treat==0, 1 / cumpr_notreat,
    Treat==1 & t_intrv < 12, 1,
    Treat==1 & t_intrv == 12 & t_treat < 12, 1,
    Treat==1 & t_intrv==12 & t_treat==12, 1 / (1-cumpr_notreat),
    Treat==1 & t_intrv==12 & t_treat>12, 0,
    Treat==1 & t_intrv > 12, 1,
    Treat==1 & t_intrv == 12 & is.na(t_treat), 0,
    Treat==1 & t_intrv > 12 & is.na(t_treat), 1
  )]
  
  # Create weights (ipw1)
  dta_3[Treat == 1, ipw1 := cumprod(ipw), by = .(clientid, Treat)]
  dta_3[Treat == 0, ipw1 := ipw, by = .(clientid, Treat)]
  
  d_glm_pe = glm(outcome_dementia ~ poly(t_intrv, 2, raw=T)*Treat + hipfrac_hat, data=dta_3,
                    family=binomial(), weights = ipw1*freqwt)
  dta_3$pr_ev = d_glm_pe$fitted.values
  dta_3[, pr_surv := cumprod(1 - pr_ev), by=list(clientid, Treat)]
  
  d_res = dta_3 %>%
    group_by(Treat, t_intrv) %>%
    summarize(pr_ev = weighted.mean(1-pr_surv, w = freqwt),
              .groups = 'drop') %>%
    ungroup %>%
    pivot_wider(., id_cols =c('t_intrv'),
                names_from = Treat,
                names_prefix = 'pr_ev_',
                values_from = pr_ev
    ) %>%
    mutate(cid = pr_ev_1 - pr_ev_0,
           cir = pr_ev_1 / pr_ev_0)
  
  
  return(d_res)
}

d_ids_1= distinct(dementia_outcome1_join, clientid)

# Create empty list to store results of each bootstrap replicate

results_d_bs <- list()

# Create do loop to run bootstrap using future_map command 1..x times

for (i in 1:500) {
  d_bs = future_map(.x = 1:d_output$runplan$boots,
                    .f = ~d_fun_getplrwt(
                      dementia_hipfrac_join_sort,
                      expanded_vaccine_0, d_ids_1,
                      .x),
                    .options = furrr_options(seed = T))
  
  d_bs_48 <- bind_rows(d_bs) %>% filter(t_intrv == 48) # Store estimate from Interval #48 (4 years)
  
  results_d_bs[[i]] <- d_bs_48
  
  write.csv(d_bs_48, file=paste0('d_bs_', i, '.csv'), row.names=FALSE)
}

d_surv = bind_rows(results_d_bs)

write.csv(d_surv,file='bs_reps.csv', row.names=FALSE)

# Calculate 95% CIs
d_summ_surv1 = d_surv %>%
  group_by(t_intrv) %>%
  summarize(pr_ev_0_lc = quantile(pr_ev_0, 0.025, na.rm=T),
            pr_ev_1_lc = quantile(pr_ev_1, 0.025, na.rm=T),
            cir_lc = quantile(cir, 0.025, na.rm=T),
            cid_lc = quantile(cid, 0.025, na.rm=T),
            pr_ev_0_uc = quantile(pr_ev_0, 0.975, na.rm=T),
            pr_ev_1_uc = quantile(pr_ev_1, 0.975, na.rm=T),
            cir_uc = quantile(cir, 0.975, na.rm=T),
            cid_uc = quantile(cid, 0.975, na.rm=T)) 

write.csv(d_summ_surv1,file='dementia1_bl_tv_95ci.csv', row.names=FALSE)

# Calculate 99% CIs
d_summ_surv2 = d_surv %>%
  group_by(t_intrv) %>%
  summarize(pr_ev_0_lc = quantile(pr_ev_0, 0.005, na.rm=T),
            pr_ev_1_lc = quantile(pr_ev_1, 0.005, na.rm=T),
            cir_lc = quantile(cir, 0.005, na.rm=T),
            cid_lc = quantile(cid, 0.005, na.rm=T),
            pr_ev_0_uc = quantile(pr_ev_0, 0.995, na.rm=T),
            pr_ev_1_uc = quantile(pr_ev_1, 0.995, na.rm=T),
            cir_uc = quantile(cir, 0.995, na.rm=T),
            cid_uc = quantile(cid, 0.995, na.rm=T)) 

write.csv(d_summ_surv2,file='dementia1_bl_tv_99ci.csv', row.names=FALSE)

