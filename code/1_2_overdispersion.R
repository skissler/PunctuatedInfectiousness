library(tidyverse) 
source("code/0_utils.R")

# ==============================================================================
# Overdispersion from punctuated infectiousness x periodic contacts 
# 
# Corresponds to writeup sections: 
#   - "The impact of the individual infectiousness profile on uncontrolled 
#      epidemic dynamics 
# 
# - Computes offspring distribution for a bunch of one-step transmissions 
# - Plots a heatmap of k as a function of punctuation (kappa) and contact
#   amplitude
# - Plots a line plot of k as a function of punctuation (kappa) for a few 
#   contact amplitudes 
# - Plots histograms of offspring distributions across punctuation (kappa)
#   values with a Poisson reference 
# - Plots epidemic trajectories and computes extinction probability, final 
#   size, and time to 100 infected/time to peak 
# ==============================================================================

