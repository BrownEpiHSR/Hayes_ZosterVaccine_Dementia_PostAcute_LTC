## McGowan Negative Control Outcome–Based On Quantitative Bias Analysis (QBA)

## Overview
This repository includes a simple quantitative bias analysis (QBA) that uses negative control outcomes to assess and correct for residual confounding in the primary 
recombinant herpes zoster vaccine (RZV) –dementia analysis.

## Purpose
#### The purpose of this analysis is to evaluate:
1.	To quantify the extent of residual confounding for the vaccine-dementia associations of interest.
2.	To use the magnitudes of association between the RZV treatment strategy and negative control outcomes to calibrate the main effects.
      In other words, to correct the main estimates according to the extent of residual confounding.

Under the assumptions of negative controls, this approach provides a structured sensitivity analysis to quantify how robust the primary findings are to unmeasured or residual confounding. 

## Conceptual Framework
  * RR_observed = observed risk ratio from the primary dementia analysis
  * RR_hip_fracture = observed risk ratio for the association between RZV and hip fracture (negative control outcome 1)
  * RR_wellness_visit = observed risk ratio for the association between RZV and wellness visit (negative control outcome 2)
  * RR_calibrated_1 = original risk ratio from the primary analysis calibrated according to the extent of residual confound found in negative control 1
  * RR_calibrated_2 = original risk ratio from the primary analysis calibrated according to the extent of residual confound found in negative control 1

We estimated the strength of an unmeasured and normally distributed continuous confounding that, if controlled for, would move the negative control RRs to the null (RR ≈ 1). We then apply those inferred confounding parameters to the primary dementia RR to obtain a bias-calibrated estimate.

## Software

The QBA is implemented using the R package:

tipr: The package provides functions to adjust observed risk ratios under assumed confounding structures.

Primary Dementia Effect Estimates Used for Calibration

From the main clone–censor–weight (CCW) + IPW analysis:

* RR_observed = 0.76
* 95% confidence interval of RR_observed = 0.69 – 0.84
  
These values are treated as the “observed” effects to be calibrated.

# Negative Control Outcomes

Two negative control outcomes are evaluated:

 #### 1. Hip Fracture
 * RR_hip_fracture ) = 0.86

 * Interpretation: Individuals who received a recombinant herpes zoster vaccine within 12 months of admission to a skilled nursing
   facility had a 14% lower risk of hip fracture compared to individuals who never received the RZV. 

 #### 2. Wellness Visits
 * Point Estimate (RR) = 1.08 (95% CI 1.00–1.15)
 
 * Interpretation: Individuals who received a recombinant herpes zoster vaccine within 12 months of admission to a skilled nursing facility
     had a 8% greater risk of getting a wellness visit compared to individuals who never received the RZV. 

## Method 1: Continuous Confounder Model

We assume an unmeasured continuous and normally distributed confounder characterized by:
* Unmeasured confounder-Exposure strength of association (x)
  
  * Test X values ranged between 0.50 to 4.0
  * X represents the scaled mean difference of the unmeasured confounder between the exposed and unexposed populations.

* Unmeasured confounder –Outcome strength of association (y)

    * Test Y values ranged between 0.50 and 4.0
    * Y represents the estimated relationship between the confounder and binary outcome on the relative scale
        (i.e., a value of 1.5 would represent a 50% greater risk in the outcome per every one-unit increase in the continuous unmeasured confounder). 

1. Generate 10,000 random combinations of (x, y).
2. Pass each x-y pair of associations through the tipr package to identify which pair attenuated the RR_hip_fracture and RR_wellness_visit to 1.0. 
3. Apply those x-y bias parameter values to adjust the primary dementia RR. Under all relevant assumptions, this yields a calibrated RR that reflects
   the degree of confounding needed to fully explain the negative control association.

## Assumptions 
   1. We assume that the relationship between RZV and each negative control outcome captures all residual confounding present in the primary RZV-dementia association.
   2. We assume that the bias parameters used to nullify the negative control outcome analyses correspond with the residual confounding structure and extent in
         the primary RZV-dementia association. 

## Method 2: Binary Confounder Model

As an alternative specification, we assume a binary unmeasured confounder defined by:
* a = prevalence among vaccinated
* b = prevalence among unvaccinated
* c = confounder–outcome effect (RR scale)
  
1. Randomly sample 10,000 (a, b, c) combinations.
2. Identify parameter sets that move the negative control RR to ~1.
3. Apply those parameters to the primary dementia RR.
This provides an alternative calibration under a discrete confounding structure.

### Assumptions -- If we want to include any ???
### Limitations -- If we want to include any ???


