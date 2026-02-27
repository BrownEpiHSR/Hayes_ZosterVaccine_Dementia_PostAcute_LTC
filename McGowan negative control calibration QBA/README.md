# McGowan Negative Control Outcome–Based On Quantitative Bias Analysis (QBA)

## Overview
This repository includes a quantitative bias analysis (QBA) that uses negative control outcomes to assess the potential impact of residual confounding in the primary vaccine–dementia analysis.

## Purpose
#### The purpose of this analysis is to evaluate:
If residual confounding is strong enough to produce the observed association between vaccination and a negative control outcome, how much could that same confounding shift the primary dementia risk estimate?

This approach does not prove causality. Instead, it provides a structured sensitivity analysis to quantify how robust the primary findings are to unmeasured or residual confounding.

## Conceptual Framework
  * RRobs = observed risk ratio from the primary dementia analysis
  * RRnc = observed risk ratio for a negative control outcome  

We estimate the strength of confounding required to move the negative control RR to the null (RR ≈ 1). We then apply those inferred confounding parameters to the primary dementia RR to obtain a bias-calibrated estimate.

## Software

The QBA is implemented using the R package:

tipr: The package provides functions to adjust observed risk ratios under assumed confounding structures.

Primary Dementia Effect Estimates Used for Calibration

From the main clone–censor–weight (CCW) + IPW analysis:

* Point estimate: RR = 0.76
* 95% CI: 0.69 – 0.84
  
These values are treated as the “observed” effects to be calibrated.

# Negative Control Outcomes

Two negative control outcomes are evaluated:

 #### 1. Hip Fracture
 * Point Estimate (RR) = 0.86

 #### 2. Wellness Visits
 * Point Estimate (RR) = 1.08 (95% CI 1.00–1.15)

## Method 1: Continuous Confounder Model

We assume an unmeasured continuous confounder characterized by:
* Exposure–confounder association (x)
* Confounder–outcome association (y)

1. Generate 10,000 random combinations of (x, y).
2. For each pair, adjust the negative control RR.
3. Identify parameter combinations that move the negative control RR to approximately 1.
4. Apply those parameter values to adjust the primary dementia RR.
This yields a calibrated RR that reflects the degree of confounding needed to fully explain the negative control association.

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


