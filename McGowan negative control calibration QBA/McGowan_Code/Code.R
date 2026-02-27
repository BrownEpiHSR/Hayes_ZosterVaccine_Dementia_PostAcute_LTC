#########################################
# Simple Quantitative Bias Analysis (QBA) 
# Using Negative Control Outcomes
#########################################

# Purpose:

#   This script performs a calibration-style QBA using negative control outcomes.
#   The idea is to use an observed association between the exposure (RZV vaccination)
#   and an outcome that should not plausibly be causally affected by RZV (negative control),
#   to infer the strength of residual confounding (or bias) that could remain in the primary study.

install.packages("tipr")
library(tipr)

options(scipen = 999)

#######################################
#Simple QBA for a continuous confounder
#######################################

#Primary study's risk ratios and lower/upper confidence bounds to calibrate
rzv_rr <-0.76
rzv_lb <-0.69
rzv_ub <-0.84


##############################################
#Hip fracture outcome
#Point estimate = 0.86
##############################################

set.seed(1234) 

n_points <- 10000

hip_fract_rr <-0.86

# Generate 10000 points with x and y between 0.01 and 5.0
random_grid_hip <- data.frame(
  x = runif(n_points, min = 0.50, max = 4),
  y = runif(n_points, min = 0.50, max = 4)
)

# For each X and Y combination, retain the adjusted RR
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

#Identify the values of X and Y where adjusted RR is close to 1.0 and retain the confounder strengths that render the NCO point estimate 1.0
null_param_hip <- subset(random_grid_hip, rr_adjusted > 0.999 & rr_adjusted < 1.001, select = c(x, y,rr_adjusted))
null_param_hip$diff <-abs(1- null_param_hip$rr_adjusted)
null_param_hip <- subset(null_param_hip, diff == min(diff),select = c(x, y, rr_adjusted))

#Pass the bias parameters through the tipr package on the main study's RR, LB and UB of the CIs
hip_frac_pe <- adjust_rr_with_continuous(effect_observed = rzv_rr, exposure_confounder_effect = null_param_hip$x, confounder_outcome_effect = null_param_hip$y)
hip_frac_lb <- adjust_rr_with_continuous(effect_observed = rzv_lb, exposure_confounder_effect = null_param_hip$x, confounder_outcome_effect = null_param_hip$y)
hip_frac_up <- adjust_rr_with_continuous(effect_observed = rzv_ub, exposure_confounder_effect = null_param_hip$x, confounder_outcome_effect = null_param_hip$y)

#Final calibrated estimates
hip_frac_pe$rr_adjusted
hip_frac_lb$rr_adjusted
hip_frac_up$rr_adjusted


##############################################
#Wellness visit
#Association = 1.08 (95%CI=1.00-1.15)
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

#Identify the values of X and Y where adjusted RR is close to 1.0. 
null_param_wellness <- subset(random_grid_wellness, rr_adjusted > 0.999 & rr_adjusted < 1.001, select = c(x, y,rr_adjusted))
null_param_wellness$diff <- abs(1- null_param_wellness$rr_adjusted)
null_param_wellness <- subset(null_param_wellness, diff == min(diff),select = c(x, y, rr_adjusted))
null_param_wellness$y <- 1/null_param_wellness$y #Take the inverse because the main study RR is protective

null_param_wellness

wellness_pe <- adjust_rr_with_continuous(effect_observed = rzv_rr, exposure_confounder_effect = null_param_wellness$x, confounder_outcome_effect = null_param_wellness$y) 
wellness_lb <- adjust_rr_with_continuous(effect_observed = rzv_lb, exposure_confounder_effect = null_param_wellness$x, confounder_outcome_effect = null_param_wellness$y) 
wellness_up <- adjust_rr_with_continuous(effect_observed = rzv_ub, exposure_confounder_effect = null_param_wellness$x, confounder_outcome_effect = null_param_wellness$y) 

wellness_pe$rr_adjusted
wellness_lb$rr_adjusted
wellness_up$rr_adjusted


##############################
#Final corrected results
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


################################
#Simple QBA for a binary confounder
################################

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

#Identify the values of X and Y where adjusted RR is close to 1.0
null_parameters <- subset(random_grid3, rr_adjusted > 0.999 & rr_adjusted < 1.001, select = c(a, b,c, rr_adjusted))
null_parameters$diff <- abs(1-null_parameters$rr_adjusted)
null_parameters <- subset(null_parameters, diff == min(diff),select = c(a, b, c, rr_adjusted))

hip_frac_pe <- adjust_rr_with_binary(effect_observed = rzv_rr, exposed_confounder_prev = null_parameters$a, unexposed_confounder_prev = null_parameters$b, confounder_outcome_effect= null_parameters$c)
hip_frac_lb <- adjust_rr_with_binary(effect_observed = rzv_lb, exposed_confounder_prev = null_parameters$a, unexposed_confounder_prev = null_parameters$b, confounder_outcome_effect= null_parameters$c)
hip_frac_up <- adjust_rr_with_binary(effect_observed = rzv_ub, exposed_confounder_prev = null_parameters$a, unexposed_confounder_prev = null_parameters$b, confounder_outcome_effect= null_parameters$c)

hip_frac_pe$rr_adjusted
hip_frac_lb$rr_adjusted
hip_frac_up$rr_adjusted





