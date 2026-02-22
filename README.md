# IN PROGRESS: Vaccination and incident dementia: Clone-censor-weights-and-bootstrapping-in-R
### R code to construct IPW, fit weighted pooled logistic models, and estimate dementia cumulative incidence with bootstrap confidence intervals.

## Overview
This repository contains code to estimate the association between vaccination and incident dementia using longitudinal administrative claims data.
The analysis uses person-time panel datasets, inverse probability weighting (IPW), pooled logistic regression, and bootstrap methods to estimate:
* Cumulative incidence of incident dementia
* Cumulative incidence difference (CID)
* Cumulative incidence ratio (CIR)
* 95% bootstrap confidence intervals

Throughout this project:
* Incident dementia refers to the diagnosis of dementia.
* Remaining free of incident dementia (dementia-free) refers to survival without dementia.
* Cumulative incidence is computed as 1 - pr\_surv
* 95% bootstrap confidence intervals are based on the 2.5th and 97.5th percentiles.

## Data Structure:
### Dataset 1: Treatment Panel (Exposure dataset)
**File:** `expanded_vaccine_cov.sas7bdat`

**Structure**
* Included only unassigned clone
* Person-month panel dataset
* Each row represents one follow-up interval (t_intrv)
* Includes baseline and time-varying covariates
  
Follow-up ends at first occurrence of
1. Shingrix vaccination
2. Death
3. Disenrollment
4. End of 4-year follow-up window 
5. End of data (December 31, 2022)
6. Hospice entry
   
Purpose
* Model vaccination probability over time
* Generate inverse probability weights

### Dataset 2: Outcome Panel
**File:** `dementia_outcome1_final.sas7bdat`

**Structure**
* Included both set of clones assigned and unassigned 
* Person-month panel dataset
* Contains treatment strategy indicator (Treat)
* Contains vaccination timing (t_treat)
* Includes baseline and time-varying covariates
  
Purpose
* Estimate dementia incidence
* Evaluate vaccine exposure effect

### Clone–Censor–Weight (CCW) design for 2 treatment startegies for 'outcome data panel' below:
### Overview of the Approach

 The CCW method consists of following steps:
1. Clone individuals into treatment strategies
2. Censor individuals when they deviate from their assigned strategy

### Step 1: Clone
Each individual is duplicated and assigned to both strategies:
* Strategy 1: No vaccination
* Strategy 2: Vaccination within the grace period

This ensures:
* Both strategies start at the same baseline
* Follow-up begins identically

This is why the outcome dataset contains:

Treat = 0, Treat = 1

### Step 2: Censor
Artificial censoring occurs when an individual deviates from the assigned strategy.

### Censoring Rules Implemented
Strategy: No Vaccination (Treat = 0)

Individuals are censored at minimum of:
* Death
* Disenrollment
* End of 4-year follow-up window 
* Hospice
* End of data (December 31, 2022)
* Dementia diagnosis
* Vaccination date 

Strategy: Vaccination (Treat = 1). A grace period is allowed (12 months)

For individuals assigned to have vaccine but they didnt have the vaccine at grace period.

Individuals are censored at minimum of:
* End of 12 month grace period
* Death
* Disenrollment
* End of 4-year follow-up window 
* End of data (December 31, 2022)
* Hospice
* Dementia diagnosis
* Vaccination date (individuals can have vaccine after 12 months)

For individuals assigned to have vaccine and have the vaccine at grace period.

Individuals are censored at minimum of:
* Death
* Disenrollment
* End of 4-year follow-up window 
* End of data (December 31, 2022)
* Hospice
* Dementia diagnosis
   
If dementia event occurs before censoring:
* Dementia event is counted

## Methods Overview:

The analytic workflow proceeds as follows:
### 1. Treatment Model
   
A pooled logistic regression estimates the probability of vaccination at each time interval using baseline and time-varying covariates.

This produces:
* Predicted treatment probabilities (pr_treat)
* Cumulative probability of remaining untreated (cumpr_notreat)
  
### 2. Inverse Probability Weights (IPW)
  
Weights are constructed to account for:
* Time-varying treatment assignment
* Grace-period design
* Artificial censoring under treatment strategies
  
Final weights:
* ipw1 (carried forward where appropriate)
  
### 3. Outcome Model
  
A weighted pooled logistic regression estimates:
P(incident dementia at time t)

Predicted values:
* pr_ev =  probability of incident dementia at each interval
* pr_surv = probability of remaining dementia-free

### 4. Cumulative Incidence Estimation 
Cumulative incidence is computed as: 1−pr_surv

Aggregated by:
* Treatment strategy
* Follow-up interval
  
Effect measures:
* CID = cumulative incidence difference
* CIR = cumulative incidence ratio

### 5. Bootstrap Inference
Uncertainty is estimated using:
* Poisson bootstrap
* Person-level resampling
* Parallel processing
  
Bootstrap results are used to compute:
* 95% confidence intervals (percentile method)

### 6. Visualization
The script produces:
1. Cumulative incidence curves
2. Cumulative incidence ratio plots
3. 95% bootstrap confidence bands
