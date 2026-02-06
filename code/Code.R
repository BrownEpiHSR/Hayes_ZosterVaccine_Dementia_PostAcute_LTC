###############################################################
# Project: Vaccine Exposure and Dementia Risk Study

# DESCRIPTION:
# This script analyzes the association between vaccination and incident dementia using longitudinal administrative claims data.

# This analysis uses two longitudinal panel datasets organized at the person-time level with baseline and time-varying covariates.

# DATA SOURCES:

# DATASET 1: Treatment (exposure) panel dataset (expanded_vaccine_cov.sas7bdat)
# ------------------------------------------------------------

# Structure:
#   - Person-time monthly panel dataset capturing vaccination exposure
#   - Each row represents one time interval of follow-up (t_intrv)
#   - Includes baseline and time-varying covariates

# Follow-up ends at first occurrence of:
#   1. Shringrix Vaccination 
#   2. Death
#   3. Disenrollment
#   4. End of 4-year after admission assessment 
#   5. End of data (December, 2022)
#   6. Hospice entry
