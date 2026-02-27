### Repository Overview

This repository contains three complementary analytical components examining the association between recombinant herpes zoster vaccination (RHZ) and incident dementia using a target trial emulation framework with clone–censor–weight methods and negative control calibration strategies.

This repository is organized into three main analytical folders. Each folder contains:

The folders are organized as follows:

1. Primary Analysis

What it contains:
* Code
     * R scripts used for clone–censor–weight implementation, model estimation, and bootstrapping
* Data Documentation
     * Excel files documenting the project details, data sources, covariates definitions used in the analysis

Purpose: 

Implements the main target trial emulation evaluating the association between RHZ vaccination and incident dementia.

2. McGowan Negative Control Calibration QBA

What it contains:
* Code
    * R scripts implementing simulation-based bias parameter exploration
    * Calibration of observed risk ratios under specified confounding scenarios
    * Graphical and numerical sensitivity analyses

* README.md
    * Description of the QBA framework
    * Explanation of bias parameters
 
Purpose: 

Implements the McGowan quantitative bias analysis (QBA) negative control calibration approach.

3. Tchetgen Tchetgen et al. Double Negative Control Calibration

What it contains:
* Code
    * R scripts for negative control exposure modeling (wellness visit)
    * Hip fracture outcome modeling
    * Construction of predicted hip fracture probabilities (hipfrac_hat)
    * Inclusion of predicted negative control outcome in the dementia model
    * Cluster-level Poisson bootstrap
    * 95% and 99% confidence interval estimation

* README.md
    * Methodological explanation of the double negative control approach
    * Description of the clone–censor adaptation

Purpose:
 Implements the regression-based double negative control method proposed by Tchetgen Tchetgen et al., adapted to the clone–censor–weight framework.



