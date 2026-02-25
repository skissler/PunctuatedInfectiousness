library(tidyverse)
source("code/utils.R")

# ==============================================================================
# Impact of punctuated infectiousness on test-based screening and isolation
#
# Corresponds to writeup section:
#   "The impact of punctuated infectiousness on test-based screening and
#    isolation countermeasures"
#
# TODO: Implement analysis of how test sensitivity interacts with punctuated
# individual infectiousness profiles. Key questions:
#   - How does the timing of testing relative to infection affect detection?
#   - Does punctuated infectiousness make periodic screening more or less
#     effective at reducing transmission?
#   - What is the optimal screening interval as a function of kappa?
# ==============================================================================
