
#################################################################
# Project: Shingrix vaccine Exposure and Dementia Risk Study

# DESCRIPTION:
# This script estimates the association between shingrix vaccination and incident dementia using longitudinal administrative 
# claims data organized as person-time (month) panels.

# DATA SOURCES:

# DATASET 1: Treatment (exposure) panel dataset (expanded_vaccine_cov.sas7bdat)
# -----------------------------------------------------------------------------
# Structure:
#   - Person-month panel dataset capturing vaccination exposure
#   - Each row represents one time interval of follow-up (t_intrv)
#   - Includes baseline and time-varying covariates

# Follow-up ends at first occurrence of:
#   1. Shringrix Vaccination 
#   2. Death
#   3. Disenrollment
#   4. End of 4-year after admission assessment 
#   5. End of data (December, 2022)
#   6. Hospice entry

# PURPOSE:
#   - Used to model probability of vaccination over time
#   - Supports estimation of treatment assignment mechanisms
#   - Used to generate inverse probability weights

# DATASET 2: Outcome Panel Dataset (dementia_outcome1_final.sas7bdat)
# ------------------------------------------------------------
# Structure:
#   - Person month panel dataset 
#   - Each row represents one time interval of follow-up (t_intrv)
#   - Contains t_treat variable which is the monthly interval for having vaccinated
#   - Contains baseline and time-varying covariates

# Follow-up ends at first occurrence of:
#   1. Dementia diagnosis
#   2. Death
#   3. Disenrollment
#   4. End of 4-year after admission assessment 
#   5. End of data (December, 2022)
#   6. Hospice entry

#    PURPOSE:
#   - Used to estimate dementia incidence risk
#   - Allows evaluation of vaccine exposure effect on dementia

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

# -----------------------------------------------------------
# Package List Definition
# -----------------------------------------------------------
# Purpose:
#   Defines a centralized list of R packages required for this analysis. Storing package names in a vector allows
#   consistent loading, installation, and dependency management across scripts.

# Why this approach:
#   - Improves reproducibility
#   - Avoids repeated library() calls scattered across scripts
#   - Makes dependency tracking easier for collaborators
#   - Supports automated installation/loading workflows

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

# Read exposure (vaccine)treatment) SAS dataset expanded_vaccine_cov.sas7bdat into R
# NULL indicates that all variables from the SAS file should be imported.
expanded_vaccine_0 <- read_sas("P:/your folder path/expanded_vaccine_cov.sas7bdat",  NULL)

# --------------------------------------------
# Parallel / bootstrap configuration options
# --------------------------------------------
# Increase maximum allowed size of global variables used in parallel processing
# Setting to +Inf removes memory size limits, allowing large objects to be passed
# to parallel workers when using the future package.
options(future.globals.maxSize = +Inf) 

# Define custom function for printing survival plots
# grid.draw.ggsurvplot() is a wrapper around survminer’s internal print method.
# It prints ggsurvplot objects without creating a new graphics page, which is useful when embedding plots into existing layouts or reports.
grid.draw.ggsurvplot = function(x) {
  survminer:::print.ggsurvplot(x, newpage=F)
}

# Analysis runtime settings used throughout (including bootstrapping).
d_output = list(runtime = Sys.time(),
                runplan = list(boots = 500,
                               workers = 2,
                               seed = as.integer(ymd('2025-04-16'))
                ))

# Set the random number generator seed using the value stored in the run settings.
# This ensures that any stochastic procedures (e.g., bootstrapping, sampling) produce reproducible results across runs.
set.seed(d_output$runplan$seed)

# ------------------------------------------------------------
# Treatment model: pooled logistic regression for vaccination
# ------------------------------------------------------------
# Fits a pooled logistic regression model to estimate the probability of receiving vaccination ( vaccine as an outcome) 
# during each time interval.

# Model specification:
#   outcome_vaccine = binary indicator for vaccination during interval
#   Predictors include:
#       • Baseline and Time-varying covariates

# poly(t_intrv, 2, raw = TRUE):
#   Includes quadratic time terms to flexibly model changes in
#   vaccination probability across follow-up intervals.

# as.factor():
#   Converts categorical variables into factor variables so that glm() estimates category-specific effects.

# family = binomial():
#   Specifies logistic regression for binary outcome modeling.

d_glm_wt1 = glm(outcome_vaccine ~ poly(t_intrv, 2, raw=T) + as.factor(BL_SEX) + poly(BL_AGE, 2) + as.factor(BL_ADL_CAT) + BL_ANEMIA +
BL_ANALGESICS + BL_ANTICOAGULANTS + BL_ANTICONVULSANTS + BL_ANTIDEPRESSANTS + BL_ANTIPSYCHOTICS + BL_ANTIVIRALS + BL_ARTHRITIS +
BL_AFIB + BL_BENZODIAZEPINES + BL_CANCER + as.factor(BL_M3_CFS) + BL_CAD + BL_COVIDVAX + BL_DVT + BL_DIABETES + as.factor(BL_ED_VISITS) + 
BL_FLUVAX + BL_GABAPENTANOIDS + as.factor(BL_GAGNE_CAT) + BL_HD + BL_HERPESZOSTER + BL_HCHOL + as.factor(BL_HOSP_VISITS) + BL_HTN + 
BL_INSULIN + BL_MENTALHEALTH + as.factor(BL_MONTH) + BL_OSTEOFRAC + BL_PARKINSONS + BL_PVD + BL_PNEUMOVAX + BL_NONHIPFX + 
as.factor(BL_RACE) + BL_RENALF + BL_STROKE + as.factor(BL_REGION) + as.factor (BL_YEAR) + BL_ZOSTAVAX + BL_PARTBCLAIM + TV_ANTICOAGULANTS +
TV_ANTIDEPRESSANTS + TV_ANTIPSYCHOTICS + TV_ANTIVIRALS + TV_CANCER + TV_COVIDVAX + TV_DVT + TV_DIABETES + TV_ED_VISITS + TV_FLUVAX + 
TV_HERPESZOSTER + TV_HOSP_VISITS + TV_INSULIN + TV_DISCHARGE + TV_PNEUMOVAX + TV_RENALF + TV_STROKE + TV_ZOSTAVAX, 
data=expanded_vaccine_0, family=binomial())

# ------------------------------------------------------------
# Export Odds Ratio Results to CSV File
# ------------------------------------------------------------
write.csv(
  data.frame(
    variable   = names(odds_ratio),
    odds_ratio = as.numeric(odds_ratio)
  ),
  file = "P:/file path/filename.csv",
  row.names = FALSE
)

# Add predicted probabilities to exposure dataset 
# Extracts fitted (predicted) probabilities from the logistic model to the exposure dataset and stores the values as a new variable (pr_treat)
expanded_vaccine_0$pr_treat = d_glm_wt1$fitted.values  

# Convert dataset to data.table format and converts the object in-place to a data.table, 
# enabling faster data manipulation and memory-efficient operations.
setDT(expanded_vaccine_0)

# Compute cumulative probability of remaining untreated up to each time interval:
# cumpr_notreat contains the cumulative probability that an individual has remained untreated up to each time interval.
# For each clientid group, computes the cumulative product of (1 - pr_treat) and stores the result as a new column named cumpr_notreat.
expanded_vaccine_0[, cumpr_notreat := cumprod(1-pr_treat), by = .(clientid)]

# Imports dementia outcome panel dataset which contains both unassigned (not to have vaccine) treat = 0, and 
# assigned (to have vaccine) treat = 1
Dementia_outcome_1 <- read_sas("P:/your folder path/dementia_outcome1_final.sas7bdat", NULL)

# Bring the cumulative probability from the treatment (exposure) dataset to the dementia outcome dataset and 
# Include only necessary variables from exposure dataset. Match by clientid and time intervals across datasets
dementia_outcome1_join = left_join(Dementia_outcome_1,
                                   select(expanded_vaccine_0, clientid, t_intrv, cumpr_notreat),
                                   by=c('clientid', 't_intrv'))

# Converts dementia_outcome1_join into a data.table object
# Conversion happens in-place (no data copy created) which enables faster filtering, grouping, and column updates
setDT(dementia_outcome1_join)

# ---------------------------------------------------------------------
# Construct Inverse Probability Weights (IPW) with Grace-Period Rules
# ---------------------------------------------------------------------
# Below block creates an IPW variable (`ipw`) in the outcome panel using cumulative probability from exposure dataset (`cumpr_notreat`) 

# Key inputs used here:
#   - Treat == 0: untreated strategy 
#   - Treat == 1: treated strategy with a 12 month grace period
#   - t_intrv: follow-up interval index (e.g., month number)
#   - t_treat: interval when vaccination actually occurs 
#   - cumpr_notreat: cumulative probability from treatment (exposure) dataset (computed earlier as cumprod(1 - pr_treat) 

# Weights construction:

# NOTE: This is project-specific logic; changes to grace period rules or treatment

#   - For Treat==0 observations, we weight by the inverse of the probability of remaining unvaccinated up to time t (1 / cumpr_notreat).
#   - For Treat==1 observations, a grace period is applied (here, t_intrv < 12)
#   - Individuals are not censored for not yet receiving treatment during the grace window, so weights are set to 1 in that period (cant censor at grace)
#   - At the end of the grace period (t_intrv == 12): weights/censoring depend on whether vaccination occurred by the grace endpoint.
#   - After grace (t_intrv > 12), weights are generally set to 1 because the treated strategy no longer censors based on "not yet treated" (per this design)

#   - fcase() is a function from the data.table package used to create variables based on multiple conditional rules

dementia_outcome1_join[, ipw := fcase(
  Treat==0, 1 / cumpr_notreat,
  Treat==1 & t_intrv < 12, 1, # trt - cant censor prior to grace
  Treat==1 & t_intrv == 12 & t_treat < 12, 1, # treat=1 & treat time < grace end, cant censor at grace
  Treat==1 & t_intrv==12 & t_treat==12, 1 / (1-cumpr_notreat), # treated at grace, then weight: "sensitivity analysis is done using the weight of 1/pr_treat"
  Treat==1 & t_intrv==12 & t_treat>12, 0, # treat=1, and not treated before grace, censor
  Treat==1 & t_intrv > 12, 1, # treat=1, after grace period, assign all weights as 1
  Treat==1 & t_intrv == 12 & is.na(t_treat), 0, # treat=1, for infinite time, assign all weights as 0
  Treat==1 & t_intrv > 12 & is.na(t_treat), 1 # treat=1, for infinite time, assign all weights as 1
)]

# ------------------------------------------------------------
# Construct Final Weight Variable (ipw1)
# ------------------------------------------------------------
# For treated strategy observations above (Treat == 1):
# Carry forward the assigned inverse probability weight over time using cumulative product 
# This ensures that once a weight is assigned, it applies to all subsequent time intervals.

# For untreated strategy observations (Treat == 0):
# No carry-forward is required, so the final weight equals ipw.

dementia_outcome1_join[Treat==1, ipw1 := cumprod(ipw), by=list(clientid, Treat)]
dementia_outcome1_join[Treat==0, ipw1 := ipw, by=list(clientid, Treat)]

# Summarize distribution of final weights (ipw1)
summary(dementia_outcome1_join$ipw1)

# Summarize distribution excluding zero weights 
summary(dementia_outcome1_join$ipw1[dementia_outcome1_join$ipw1 != 0])

# ------------------------------------------------------------
# Generate Weight Distribution Summary by Treatment Group
# ------------------------------------------------------------
# Filters dataset, calculates descriptive statistics for final weights (ipw1) separately for each treatment group
dementia_outcome1_join %>%
  filter(ipw1 != 0) %>%
  group_by(Treat) %>%
  summarise(summary_stats = list(as.data.frame(as.list(summary(ipw1))))) %>%
  unnest(summary_stats)

# ------------------------------------------------------------
# Generate Weight Quantiles by Treatment Group
# ------------------------------------------------------------
# Calculates selected percentiles of the final weights (ipw1) separately for each treatment group s
dementia_outcome1_join %>%
  filter(ipw1 != 0 & !is.na(ipw1)) %>%
  group_by(Treat) %>%
  summarise(quantiles = list(quantile(ipw1, probs = c(0.01, 0.05, 0.25, 0.5, 
                                                      0.75, 0.9, 0.95, 0.99, 
                                                      0.995, 0.997, 0.998), 
                                      na.rm = TRUE))) %>%
unnest_wider(quantiles)

# Outcome model 
# ------------------------------------------------------------
# Outcome model: weighted pooled logistic regression
# ------------------------------------------------------------
# Fits a weighted pooled logistic regression model to estimate dementia risk over time with
  #   - quadratic time trend
  #   - Treat as an main effect
  #   - time-by-treatment interaction
# Incorporating inverse probability weights (ipw1).
# Models non-linear changes in outcome risk over time.

d_glm_pe2 = glm(outcome_dementia ~ poly(t_intrv, 2, raw=T)*Treat,data=dementia_outcome1_join,
                family= binomial(), weights = ipw1)

# Add fitted event probabilities to the outcome dataset:
# Extracts fitted (predicted) probabilities from the weighted outcome model and stores them as a new variable (pr_ev) in the dataset.
dementia_outcome1_join$pr_ev = d_glm_pe2$fitted.values

# After assigning pr_ev, calculate pr_surv by clientid and Treat

# Compute predicted probability of remaining free of incident dementia (dementia - free) across follow-up intervals
# Stores the result as pr_surv.
dementia_outcome1_join[, pr_surv := cumprod(1 - pr_ev), by = list(clientid, Treat)]

# ------------------------------------------------------------
# Summarize Predicted Cumulative Incidence by Treatment strategy over Time
# ------------------------------------------------------------
# This pipeline:
#   - Aggregates individual-level predictions by (Treat, t_intrv) level
#   - Computes the mean cumulative probability of incident dementia (i.e., 1 - pr_surv) at each interval
#   - Reshapes results to wide format to place treated and untreated estimates side-by-side
#   - Computes effect measures:
#        - cid: cumulative incidence difference (treated - untreated)
#        - cir: cumulative incidence ratio (treated / untreated)
d_res1 = dementia_outcome1_join %>%
  group_by(Treat, t_intrv) %>%
  summarize(pr_ev = mean(1-pr_surv), .groups = 'drop') %>%
  ungroup %>%
  pivot_wider(., id_cols =c('t_intrv'),
              names_from = Treat,
              names_prefix = 'pr_ev_',
              values_from = pr_ev
  ) %>%
  mutate(cid = pr_ev_1 - pr_ev_0,
         cir = pr_ev_1 / pr_ev_0)

# Saves summarized estimates of cumulative incidence of incident dementia to an external CSV file
write.csv(d_res1,file='P:your folder path/d_res1.csv', row.names=FALSE)

# ------------------------------------------------------------
# Configure Bootstrapping and Analysis Runtime Settings
# ------------------------------------------------------------
# Increase allowed size of global objects for parallel processing.
# This prevents memory errors when large datasets or model objects are passed to workers during bootstrapping.

options(future.globals.maxSize = +Inf) 

grid.draw.ggsurvplot = function(x) {
  survminer:::print.ggsurvplot(x, newpage=F)
}

d_output = list(runtime = Sys.time(),
                runplan = list(boots = 500,
                               workers = 2, 
                               seed = as.integer(ymd('2025-04-16'))
                ))

# Set up parallel processing using the future package.
# Multisession launches separate R sessions to run tasks in parallel. The number of parallel workers is taken from the analysis run plan.
plan(multisession, workers = d_output$runplan$workers)

# Set random seed for reproducibility. Ensures that stochastic procedures (e.g., bootstrapping) produce identical results when the analysis is re-run.
set.seed(d_output$runplan$seed)

# ==========================================================
# 1) Create Poisson bootstrap (clustered by person)
# ==========================================================

d_fun_getplrwt = function(dta_outc, dta_treat, ...) {
  
  # Person-month level Poisson(1) frequency weights for the bootstrap replicate
  d_freqwt = distinct(dta_outc, clientid) %>%
    mutate(., freqwt = rpois(n(), 1L))
  
  # Attach person-month level bootstrap weights to outcome dataset
  dta_outc = left_join(dta_outc, d_freqwt, by='clientid')
  
  # Attach person- month level bootstrap weights to treatment dataset
  dta_2 = left_join(dta_treat, d_freqwt, by='clientid')
  
  # ==========================================================
  # 2) Fit weighted treatment (exposure) model using glm()
  # ==========================================================
  # Fits a pooled logistic regression for outcome_vaccine using baseline and time-varying covariates. 
  
  d_glm_wt = glm(outcome_vaccine ~ poly(t_intrv, 2, raw=T) + as.factor(BL_SEX) + poly(BL_AGE, 2) + 
  as.factor(BL_ADL_CAT) + BL_ANEMIA + BL_ANALGESICS + BL_ANTICOAGULANTS + BL_ANTICONVULSANTS + BL_ANTIDEPRESSANTS + 
  BL_ANTIPSYCHOTICS + BL_ANTIVIRALS + BL_ARTHRITIS + BL_AFIB + BL_BENZODIAZEPINES + BL_CANCER + as.factor(BL_M3_CFS) + 
  BL_CAD + BL_COVIDVAX + BL_DVT + BL_DIABETES + as.factor(BL_ED_VISITS) + BL_FLUVAX + BL_GABAPENTANOIDS + as.factor(BL_GAGNE_CAT) + 
  BL_HD + BL_HERPESZOSTER + BL_HCHOL + as.factor(BL_HOSP_VISITS) + BL_HTN + BL_INSULIN + BL_MENTALHEALTH + as.factor(BL_MONTH) + 
  BL_OSTEOFRAC + BL_PARKINSONS + BL_PVD + BL_PNEUMOVAX + BL_NONHIPFX + as.factor(BL_RACE) + BL_RENALF + BL_STROKE + 
  as.factor(BL_REGION) + as.factor (BL_YEAR) + BL_ZOSTAVAX + BL_PARTBCLAIM + TV_ANTICOAGULANTS + TV_ANTIDEPRESSANTS + 
  TV_ANTIPSYCHOTICS + TV_ANTIVIRALS + TV_CANCER + TV_COVIDVAX + TV_DVT + TV_DIABETES + TV_ED_VISITS + TV_FLUVAX + 
  TV_HERPESZOSTER + TV_HOSP_VISITS + TV_INSULIN + TV_DISCHARGE + TV_PNEUMOVAX + TV_RENALF + TV_STROKE + TV_ZOSTAVAX, 
  data=dta_2, family=binomial(), weights=freqwt)  

  # ======================================================================================================
  # 3) Generate predicted treatment probabilities and compute cumulative probability of remaining untreated
  # ======================================================================================================
  # pr_treat: Model-fitted probability of vaccination in each interval.
  # cumpr_notreat: Cumulative probability of remaining untreated up to each interval, computed as cumprod(1 - pr_treat) within clientid.

  dta_2$pr_treat = d_glm_wt$fitted.values
  setDT(dta_2)
  dta_2[, cumpr_notreat := cumprod(1-pr_treat), by = .(clientid)]
  
  # ==========================================================
  # 4) Merge cumulative non-treatment probability into the outcome dataset
  # ==========================================================
  # Aligns treatment-history-derived variable (cumpr_notreat) with the outcome panel at the same (clientid, t_intrv).
  dta_3 = left_join(dta_outc, select(dta_2, clientid, t_intrv, cumpr_notreat), by=join_by(clientid, t_intrv))    
  setDT(dta_3)
  
  # ==========================================================
  # 5) Construct inverse probability weights (ipw) with project-specific rules
  # ==========================================================
  # fcase() assigns ipw based on mutually exclusive conditions.
  
  # Key logic:
  #   - Treat==0: ipw = 1 / cumpr_notreat
  #   - Treat==1: apply grace-period and censoring rules at t_intrv == 12

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
  # ==========================================================
  # 6) Create final weights (ipw1) by carrying forward treated-arm weights
  # ==========================================================
  # For Treat==1:
  #   ipw1 = cumprod(ipw) within (clientid, Treat), which carries forward a weight
  #   assigned at the grace boundary and ensures censoring (ipw==0) persists.
  
  # For Treat==0:
  #   ipw1 = ipw directly (no carry-forward rule applied here)
  
  dta_3[Treat == 1, ipw1 := cumprod(ipw), by = .(clientid, Treat)]
  dta_3[Treat == 0, ipw1 := ipw, by = .(clientid, Treat)]
  
  # ==========================================================
  # 7) Fit weighted outcome model using glm()
  # ==========================================================
  # Fits a pooled logistic model for outcome_dementia with:
  #   - quadratic time trend
  #   - Treat as an main effect
  #   - time-by-treatment interaction
  #
  # weights = ipw1 * freqwt:
  #   - ipw1 adjusts for treatment assignment/censoring per design rules
  #   - freqwt applies the Poisson bootstrap replicate weights
  
  d_glm_pe = glm(outcome_dementia ~ poly(t_intrv, 2, raw=T)*Treat, data=dta_3,
                 family=binomial(), weights = ipw1*freqwt)
  
  # ==========================================================
  # 8) Generate predicted event probabilities and predicted survival (remaining free of incident dementia)
  # ==========================================================
  # pr_ev:
  #   Fitted probability of incident dementia in each interval.
  # pr_surv:
  #   Predicted probability of remaining free of incident dementia up to each interval,
  #   computed as cumprod(1 - pr_ev) within (clientid, Treat).
                 
  dta_3$pr_ev = d_glm_pe$fitted.values
  dta_3[, pr_surv := cumprod(1 - pr_ev), by=list(clientid, Treat)]
  
  # ==================================================================================================
  # 9) Aggregate to treatment strategy level cumulative incidence and compute CID/CIR
  # ===================================================================================================
  # (1 - pr_surv) is the predicted cumulative incidence of incident dementia up to time t.
  # weighted.mean(..., w=freqwt) applies the Poisson bootstrap weights when aggregating.
  
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
  
  # ==========================================================
  # 10) Return bootstrap replicate results
  # ==========================================================
  
  return(d_res)
}

# Extract unique individual identifiers from the outcome dataset.
# This creates a one-row-per-person dataset that defines the resampling unit for the bootstrap procedure.

d_ids_1= distinct(dementia_outcome1_join, clientid)

# Run Poisson bootstrap replicates in parallel using furrr.
# For each bootstrap iteration:
#   - A Poisson(1) frequency weight is generated at the person month level inside d_fun_getplrwt()
#   - Treatment and outcome models are re-fitted
#   - Predicted cumulative incidence of incident dementia is returned
# future_map() executes iterations in parallel according to the previously defined future::plan().
d_bs = future_map(.x = 1:d_output$runplan$boots,
                  .f = ~d_fun_getplrwt(
                    select(dementia_outcome1_join, -cumpr_notreat),
                    expanded_vaccine_0, d_ids_1,
                    .x),
                  .options = furrr_options(seed = T))

# ------------------------------------------------------------
# Summarize bootstrap replicates into 95% confidence intervals
# ------------------------------------------------------------
# Below code:
#   1) Combines the list of bootstrap replicate results (d_bs) into a single long table
#   2) Computes 2.5th and 97.5th percentiles across bootstrap replicates (by time interval) to form 95% bootstrap uncertainty intervals for:
#        - pr_ev_0, pr_ev_1: cumulative incidence of incident dementia by arm
#        - cid: cumulative incidence difference
#        - cir: cumulative incidence ratio
#   3) Appends point estimates from d_res1 to the bootstrap interval summary


# ==========================================================
# 1) Stack bootstrap replicate outputs into one dataset
# ==========================================================
# d_bs is a list where each element is a data frame for one bootstrap replicate.
# bind_rows(.id='boot') stacks them into a single data frame and adds a 'boot' column indicating which replicate each row came from.
# mutate() converts boot from character to numeric for easier sorting/plotting.

d_surv = d_bs %>%
  bind_rows(.id = 'boot') %>%
  mutate(boot = as.numeric(boot))

# ==========================================================
# 2) Compute bootstrap 95% intervals by follow-up interval
# ==========================================================
# For each follow-up interval (t_intrv), compute:
#   - lower confidence bound (lc): 2.5th percentile across bootstrap replicates
#   - upper confidence bound (uc): 97.5th percentile across bootstrap replicates
#   - na.rm=TRUE ensures missing values do not break quantile calculations.
# Bind_cols() combines:
  #   - point estimates from d_res1 (pr_ev_0, pr_ev_1, cir, cid)

d_summ_surv1 = d_surv %>%
  group_by(t_intrv) %>%
  summarize(pr_ev_0_lc = quantile(pr_ev_0, 0.025, na.rm=T),
            pr_ev_1_lc = quantile(pr_ev_1, 0.025, na.rm=T),
            cir_lc = quantile(cir, 0.025, na.rm=T),
            cid_lc = quantile(cid, 0.025, na.rm=T),
            pr_ev_0_uc = quantile(pr_ev_0, 0.975, na.rm=T),
            pr_ev_1_uc = quantile(pr_ev_1, 0.975, na.rm=T),
            cir_uc = quantile(cir, 0.975, na.rm=T),
            cid_uc = quantile(cid, 0.975, na.rm=T)) %>%
  bind_cols(select(d_res1, pr_ev_0, pr_ev_1, cir, cid),
            .)

# ------------------------------------------------------------
# Export Final Bootstrap Summary Results to CSV
# ------------------------------------------------------------
write.csv(d_summ_surv1,file='P:/your folder path/dementia1_bl_tv.csv', row.names=FALSE)

# ------------------------------------------------------------
# Reload Results for Plotting
# ------------------------------------------------------------
# This step is typically used to:
#   - ensure plots are generated from saved outputs
d_summ_surv1 <- read.csv('P:/your folder path/dementia1_bl_tv.csv')

# ------------------------------------------------------------
# Plot Cumulative Incidence of Incident Dementia with 95% CIs
# ------------------------------------------------------------
# Creates a line plot of cumulative incidence of incident dementia over follow-up time for two treatment strategies, with shaded color
d_gg_ci = d_summ_surv1 %>%
  ggplot(aes(x=t_intrv)) +
  geom_line(aes(y = pr_ev_0, color='No RZV'), linewidth=1.2) +
  geom_line(aes(y = pr_ev_1, color='1+ RZV'), linewidth=1.2, linetype=2) +
  scale_color_manual(name='Treatment strategy', values=c('No RZV'='red','1+ RZV'='blue' )) +
  geom_ribbon(aes(ymin = pr_ev_0_lc, ymax = pr_ev_0_uc),
              fill='red', alpha=0.2) +
  geom_ribbon(aes(ymin = pr_ev_1_lc, ymax = pr_ev_1_uc),
              fill='blue', alpha=0.2) +
  scale_x_continuous(breaks = seq(0, 49, 6),
                     limits = c(0, 49)) +
  theme_bw() +
  labs(x = 'Follow-up', y = 'Cumulative incidence', color = "Legend")


# ------------------------------------------------------------
# Plot Cumulative Incidence Ratio (CIR) with 95% Bootstrap CIs
# ------------------------------------------------------------
# Creates a line plot of the cumulative incidence ratio (cir) over follow-up time, with a shaded color showing the 95% bootstrap percentile

d_gg_rr = d_summ_surv1 %>%
  ggplot(aes(x=t_intrv)) +
  geom_line(aes(y = cir), color='green', linewidth=1.2) +
  geom_ribbon(aes(ymin = cir_lc, ymax = cir_uc),
              fill='green', alpha=0.2) +
  scale_y_continuous(limits = c(0.5, 1.1)) +
  scale_x_continuous(breaks = seq(0, 49, 6),
                     limits = c(0, 49)) +
  theme_bw() +
  labs(x = 'Follow-up', y = 'Relative Risk', color = "Legend")

d_gg_1 = ggarrange(d_gg_ci, d_gg_rr,
                   nrow=1)

d_gg_1

# ------------------------------------------------------------
# Save Final Combined Figure to folder
# ------------------------------------------------------------
# Exports the combined figure (d_gg_1), which displays:
#   - Cumulative incidence of incident dementia by treatment strategy
#   - Cumulative incidence ratio with 95% bootstrap confidence intervals

# Figures are saved in both PNG (raster) and EPS (vector) formats or use in presentations and publication-quality documents.

ggsave("P:/your folder path/dementia1_bl_tv_v2.png",
       plot=d_gg_1, width = 12, height=8, dpi=300)

ggsave("P:/our folder path/dementia1_bl_tv.eps",
       plot=d_gg_1, width = 12, height=8, dpi=300, device = cairo_ps)


