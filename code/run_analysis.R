library(tidyverse) 
library(odin) 

# Setup
source('code/utils.R') 
source('code/parameters.R')

# Uncontrolled epidemics 
source('code/episims_seir.R') 
source('code/episims_gamma.R')
source('code/survival_theory.R')
source('code/variable_contacts.R')

# Controlled epidemics 
source('code/isolation.R')
source('code/gathering_sizes.R')

# Inference
source('code/identifiability.R')

