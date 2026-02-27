# Double Negative Control Calibration (Tchetgen Tchetgen et al.)

This repository implements a regression-based double negative control calibration analysis to evaluate and adjust for residual confounding in a target trial emulation of recombinant herpes zoster vaccination (rHZ) and incident dementia among skilled nursing facility residents.

To probe and calibrate potential unmeasured confounding, we implemented the double negative control method proposed by Tchetgen Tchetgen et
al. (DOI: 10.1214/23-STS911), and applied by Li et al. (https://projecteuclid.org/journals/statistical-science/volume-39/issue-3/An-Introduction-to-Proximal-Causal-Inference/10.1214/23-STS911.short), adapting it to the clone–censor weighted framework.

We acknowledge that there is currently limited methodological guidance for implementing double negative control methods in clone–censor structured longitudinal data with treatment strategies. This repository documents our adaptation for transparency.

### Double Negative Control Framework
In addition to the McGowan calibration approach, we implemented the regression-based double negative control method described by Tchetgen Tchetgen et al. and operationalized by Li et al.
This framework uses:

*	A negative control exposure
*	A negative control outcome
*	Regression calibration using predicted probabilities for the negative control outcome

### Negative Control Exposure

* Treatment Strategy

Receipt of a wellness visit within 12 months of admission to the skilled nursing facility vs. never receiving a wellness visit over follow-up.

* Clone–Censor Implementation

Each individual was cloned and assigned to:

1.	Receive a wellness visit within 12 months of admission to the skilled nursing facility
2.	Never receive a wellness visit over the course of follow-up

Clones were censored if they deviated from their assigned strategy:
* Assigned to receive within 12 months but did not → censored immediately after the 12 month grace window
* Assigned to never receive but received one → censored at visit date

Inverse probability of censoring weights were estimated using:
* Baseline covariates
* Time-varying covariates

### Negative Control Outcome
* Hip Fracture

Comparison:
   * Wellness visit within 12 months days
      vs.
   * No wellness visit over follow-up

### Calibration Mechanism
From the negative control model:

1.	A weighted pooled logistic regression was fit for hip fracture.
2.	The time-varying predicted probability of hip fracture (hipfrac_hat) was generated.
3.	These predicted probabilities were retained for each person-time interval to be included as a covariate in the main analysis outcome model. 

#### Correcting the original RZV-dementia risk ratio: 

glm(outcome_dementia ~ poly(t_intrv,2)*Treat + hipfrac_hat,
    family = binomial(),
    weights = ipw1)

This follows the double negative control logic:

  * The negative control outcome captures residual confounding structure.
  * Including predicted hip fracture probability adjusts for residual open backdoor pathways related to healthy vaccinee biases.

#### Bootstrap Procedure
We implemented a cluster-level Poisson bootstrap:
*	Bootstrap at individual level
*	Re-estimate exposure model
*	Recalculate weights
* Refit dementia model including hipfrac_hat
*	Extract month 48 cumulative incidence
*	Repeat for 100 bootstrap replicates

Confidence intervals computed using percentile method:
  * 95% CI → 2.5th / 97.5th percentiles
  *	99% CI → 0.5th / 99.5th percentiles

#### Conceptual Interpretation
The double negative control approach:
*	Tests for residual confounding
* Quantifies potential bias
*	Calibrates the primary estimate
*	Strengthens causal interpretation under clone–censor framework

#### Methodological Assumptions and Considerations
*	We assume that the negative control exposure and outcome model sufficiently captures the residual confounding of interest, 
  which cannot be empirically confirmed. We note that compared to the negative control analysis that used McGowan, the double 
  negative control analysis resulted in a smaller attenuation of the main effects, suggesting the existence of residual bias. 
*	There is limited published guidance exists for double negative control analyses implementation in clone–censor weighted. 
  It is possible that different definitions of our negative control exposure would lend to different bias corrected results. 


