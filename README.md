# IN PROGRESS: Clone-censor-weights-and-bootstrapping-in-R
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
* Estimate treatment assignment mechanism
* Generate inverse probability weights

### Dataset 2: Outcome Panel
**File:** `dementia_outcome1_final.sas7bdat`

**Structure**
* Person-month panel dataset
* Contains treatment strategy indicator (Treat)
* Contains vaccination timing (t_treat)
* Includes baseline and time-varying covariates
  
Follow-up ends at first occurrence of
1. Dementia diagnosis
2. Death
3. Disenrollment
4. End of 4-year follow-up window
5. End of data (December 31, 2022)
6. Hospice entry
   
Purpose
* Estimate dementia incidence
* Evaluate vaccine exposure effect

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
