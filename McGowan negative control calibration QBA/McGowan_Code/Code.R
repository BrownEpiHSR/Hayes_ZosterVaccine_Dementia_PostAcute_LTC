#########################################
# Simple Quantitative Bias Analysis (QBA) 
# Using Negative Control Outcomes
#########################################

# Purpose:
# - This script performs a calibration-style QBA using negative control outcomes.
# - The idea is to use an observed association between the exposure (RZV vaccination)
# - and an outcome that should not plausibly be causally affected by RZV (negative control),
# - to infer the strength of residual confounding (or bias) that could remain in the primary study

# Approach (high-level):
#   1) Take the calculated Risk Ratio (RR) for a negative control outcome (hip fracture andwellness visit)
#   2) Search for bias parameter values that would "explain away" that RR (i.e., adjust it to ~1)
#   3) Apply those inferred bias parameters to the primary dementia RR to obtain a calibrated RR

# Notes:
#   - It is a sensitivity analysis that asks:
#       "If residual confounding is strong enough to create the observed association 
#        with a negative control outcome, how much could it shift our primary effect?"
#   - The tipr package provides helper functions to adjust observed RRs under assumed confounding structures.

install.packages("tipr")
library(tipr)

# Print options: reduce scientific notation for readability
options(scipen = 999)

# -----------------------------------------------------------
# Primary study effect estimates to be calibrated
# -----------------------------------------------------------
# These are the dementia primary results from the main CCW/IPW analysis.
rzv_rr <-0.76
rzv_lb <-0.69
rzv_ub <-0.84
  
##############################################
# Negative Control Outcome 1: Hip Fracture
# Point estimate = 0.86
##############################################

set.seed(1234) 

# Number of random parameter combinations to search over
n_points <- 10000

# Observed RR between RZV and hip fracture (negative control)
hip_fract_rr <-0.86

# -----------------------------------------------------------
# Continuous-confounding QBA: parameter search
# -----------------------------------------------------------
# We search over two continuous bias parameters:
#   x = "exposure–confounder effect" (how strongly the confounder is associated with vaccination)
#   y = "confounder–outcome effect"  (how strongly the confounder is associated with the outcome)

# Goal:
#   Find (x, y) such that the adjusted RR for hip fracture is ~1.0.

# Generate 10000 points with x and y between 0.01 and 5.0
random_grid_hip <- data.frame(
  x = runif(n_points, min = 0.50, max = 4),
  y = runif(n_points, min = 0.50, max = 4)
)

# For each (x, y) combination, compute the adjusted RR for the hip fracture association.
# adjust_rr_with_continuous() returns a list; and retain the adjusted RR

random_grid_hip$rr_adjusted <- mapply(
  function(x, y) {
    adjust_rr_with_continuous(
      effect_observed = hip_fract_rr,
      exposure_confounder_effect = x,
      confounder_outcome_effect = y
    )$rr_adjusted
  },
  random_grid_hip$x,
  random_grid_hip$y
)

# Identify the values of X and Y where adjusted RR is close to 1.0 and retain the confounder strengths that 
# render the NCO point estimate 1.0

null_param_hip <- subset(random_grid_hip, rr_adjusted > 0.999 & rr_adjusted < 1.001, select = c(x, y,rr_adjusted)) 
null_param_hip$diff <-abs(1- null_param_hip$rr_adjusted) # Choose the best pair among near-null solutions (closest to 1).
null_param_hip <- subset(null_param_hip, diff == min(diff),select = c(x, y, rr_adjusted))

#Pass the bias parameters through the tipr package on the primary analysis's RR, LB and UB of the CIs
hip_frac_pe <- adjust_rr_with_continuous(effect_observed = rzv_rr, exposure_confounder_effect = null_param_hip$x, confounder_outcome_effect = null_param_hip$y)
hip_frac_lb <- adjust_rr_with_continuous(effect_observed = rzv_lb, exposure_confounder_effect = null_param_hip$x, confounder_outcome_effect = null_param_hip$y)
hip_frac_up <- adjust_rr_with_continuous(effect_observed = rzv_ub, exposure_confounder_effect = null_param_hip$x, confounder_outcome_effect = null_param_hip$y)

# Final calibrated estimates and bounds
hip_frac_pe$rr_adjusted
hip_frac_lb$rr_adjusted
hip_frac_up$rr_adjusted

##############################################
# Negative Control Outcome 2: Wellness Visit
# # Association = 1.08 (95% CI 1.00–1.15)
##############################################

set.seed(1234) 

wellness_rr <-1.08

# Generate 10000 points with x and y between 0.01 and 5.0
random_grid_wellness <- data.frame(
  x = runif(n_points, min = 0.50, max = 4),
  y = runif(n_points, min = 0.50, max = 4)
)

# For each X and Y combination, retain the adjusted RR
random_grid_wellness$rr_adjusted <- mapply(
  function(x, y) {
    adjust_rr_with_continuous(
      effect_observed = wellness_rr,
      exposure_confounder_effect = x,
      confounder_outcome_effect = y
    )$rr_adjusted
  },
  random_grid_wellness$x,
  random_grid_wellness$y
)

# Identify the values of X and Y where adjusted RR is close to 1.0. 
null_param_wellness <- subset(random_grid_wellness, rr_adjusted > 0.999 & rr_adjusted < 1.001, select = c(x, y,rr_adjusted))
null_param_wellness$diff <- abs(1- null_param_wellness$rr_adjusted)
null_param_wellness <- subset(null_param_wellness, diff == min(diff),select = c(x, y, rr_adjusted))
null_param_wellness$y <- 1/null_param_wellness$y # Take the inverse because the main study RR is protective

# IMPORTANT:
#   The primary analysis RR is protective (<1), while the wellness RR is >1.
#   To map the bias direction consistently (protective vs harmful),
#   this code inverts y (confounder–outcome effect) so the calibration aligns
#   with a protective direction in the primary analysis.

null_param_wellness$y <- 1 / null_param_wellness$y

null_param_wellness

wellness_pe <- adjust_rr_with_continuous(effect_observed = rzv_rr, exposure_confounder_effect = null_param_wellness$x, confounder_outcome_effect = null_param_wellness$y) 
wellness_lb <- adjust_rr_with_continuous(effect_observed = rzv_lb, exposure_confounder_effect = null_param_wellness$x, confounder_outcome_effect = null_param_wellness$y) 
wellness_up <- adjust_rr_with_continuous(effect_observed = rzv_ub, exposure_confounder_effect = null_param_wellness$x, confounder_outcome_effect = null_param_wellness$y) 

wellness_pe$rr_adjusted
wellness_lb$rr_adjusted
wellness_up$rr_adjusted


##############################
# Final corrected results
##############################

hip_frac_pe$rr_adjusted
hip_frac_lb$rr_adjusted
hip_frac_up$rr_adjusted

print(c(hip_frac_pe$rr_adjusted,hip_frac_lb$rr_adjusted,hip_frac_up$rr_adjusted))

wellness_pe$rr_adjusted
wellness_lb$rr_adjusted
wellness_up$rr_adjusted

print(c(wellness_pe$rr_adjusted,wellness_lb$rr_adjusted,wellness_up$rr_adjusted))

####################################################################################################################################################################################

# Purpose:
#   Repeat the calibration using tipr's binary-confounder formulation.
#   Here, the confounder is assumed to be binary, characterized by:
#     a = prevalence of confounder among exposed (vaccinated)
#     b = prevalence of confounder among unexposed (unvaccinated)
#     c = confounder–outcome effect (RR scale)

# Goal:
#   Find (a, b, c) such that the adjusted hip fracture RR is ~1,
#   then apply those parameters to the primary analysis RR.

set.seed(1234) 

n_points <- 10000

hip_fract_rr <-0.86
rzv_rr <-0.76
rzv_lb <-0.69
rzv_ub <-0.84

# Generate 10000 points with x and y between 0.01 and 5.0
random_grid3 <- data.frame(
  a = runif(n_points, min = 0.001, max = .999),
  b = runif(n_points, min = 0.001, max = .999),
  c = runif(n_points, min = 0.10, max = 4)
)

# For each X and Y combination, retain the adjusted RR
random_grid3$rr_adjusted <- mapply(
  function(a, b, c) {
    adjust_rr_with_binary(
      effect_observed = hip_fract_rr,
      exposed_confounder_prev = a,
      unexposed_confounder_prev = b,
      confounder_outcome_effect = c
    )$rr_adjusted
  },
  random_grid3$a,
  random_grid3$b,
  random_grid3$c
)

# Identify the values of X and Y where adjusted RR is close to 1.0
null_parameters <- subset(random_grid3, rr_adjusted > 0.999 & rr_adjusted < 1.001, select = c(a, b,c, rr_adjusted))
null_parameters$diff <- abs(1-null_parameters$rr_adjusted)
null_parameters <- subset(null_parameters, diff == min(diff),select = c(a, b, c, rr_adjusted))

hip_frac_pe <- adjust_rr_with_binary(effect_observed = rzv_rr, exposed_confounder_prev = null_parameters$a, unexposed_confounder_prev = null_parameters$b, confounder_outcome_effect= null_parameters$c)
hip_frac_lb <- adjust_rr_with_binary(effect_observed = rzv_lb, exposed_confounder_prev = null_parameters$a, unexposed_confounder_prev = null_parameters$b, confounder_outcome_effect= null_parameters$c)
hip_frac_up <- adjust_rr_with_binary(effect_observed = rzv_ub, exposed_confounder_prev = null_parameters$a, unexposed_confounder_prev = null_parameters$b, confounder_outcome_effect= null_parameters$c)

hip_frac_pe$rr_adjusted
hip_frac_lb$rr_adjusted
hip_frac_up$rr_adjusted





