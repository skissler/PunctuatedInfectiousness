library(tidyverse)
source("code/utils.R")

# ==============================================================================
# Identifiability of the punctuation parameter kappa
#
# Corresponds to writeup section:
#   "Identifiability of punctuated infectiousness"
#
# Simulates household / contact-tracing study designs and asks how well
# kappa can be estimated from realistic data. Key questions:
#   - Given a household study with N households, can we distinguish
#     kappa = 1 from kappa = 9?
#   - What data signatures (e.g., clustering of secondary infection times)
#     are most informative about kappa?
#   - How does sample size affect identifiability?
# ==============================================================================
