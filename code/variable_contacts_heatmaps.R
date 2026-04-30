# ==============================================================================
# Import
# ==============================================================================

library(tidyverse)
source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

pars <- parslist[[2]]
pathogen <- pars$pathogen
T        <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

c_per             <- 1
c_amp_vals_sim    <- seq(0, 1, length.out=25)
psi_vals_sim      <- seq(0, 0.5, length.out=25)
c_amp_vals_theory <- seq(0, 1, length.out=100)
psi_vals_theory   <- seq(0, 0.5, length.out=100)

n_index <- 10000
k_cap   <- 50

# ==============================================================================
# Periodic contacts: heatmaps
# ==============================================================================

# Define contact modulator: c_i(t) = 1 - c_amp*cos(2*pi*t/c_per)
# Returns R0*c_i(t), for use with gen_inf_attempts_gamma_contacts
make_contact_fn_periodic <- function(R0, c_amp, c_per) {
	stopifnot(c_amp >= 0, c_amp <= 1)
	function(t) R0*(1 - c_amp*cos(2*pi*t/c_per))
}

# Initialize a data frame for simulated k-values 
sim_grid_periodic <- expand_grid(psi=psi_vals_sim, c_amp=c_amp_vals_sim) %>%
	mutate(k_sim = NA_real_)

# Generate the simulated k-values
for(idx in 1:nrow(sim_grid_periodic)){
	psi_val   <- sim_grid_periodic$psi[idx]
	c_amp_val <- sim_grid_periodic$c_amp[idx]
	z     <- make_contact_fn_periodic(R0, c_amp_val, c_per)
	z_max <- R0*(1 + c_amp_val)
	gfun  <- gen_inf_attempts_gamma_contacts(T=T, z=z, z_max=z_max,
	                                         alpha=alpha, psi=psi_val)
	tinfs      <- c_per * runif(n_index)
	noffspring <- sapply(lapply(tinfs, gfun), length)
	m_s <- mean(noffspring)
	v_s <- var(noffspring)
	sim_grid_periodic$k_sim[idx] <- if(v_s > m_s){ m_s^2 / (v_s - m_s) } else { Inf }
}

# Generate the theoretical k-values 
theory_grid_periodic <- expand_grid(
	    psi = psi_vals_theory, 
	    c_amp = c_amp_vals_theory) %>%
	mutate(k_theory = (2/c_amp^2)*(1 + (2*pi/c_per)^2 / beta^2)^(psi*alpha))

# Add capped k column
sim_grid_periodic <- sim_grid_periodic %>%
	mutate(k_capped = pmin(k_sim, k_cap))

theory_grid_periodic <- theory_grid_periodic %>% 
	mutate(k_capped = pmin(k_theory, k_cap))

# Generate the heat map
fig_heatmap_periodic <- ggplot() +
	geom_tile(data = sim_grid_periodic, aes(x=psi, y=c_amp, fill=k_capped)) +
	geom_contour(data = theory_grid_periodic,
		aes(x=psi, y=c_amp, z=k_theory),
		col = "black", linewidth = 1.0, breaks = c(1, 2, 5, 10, 20, 50)) +
	geom_contour(data = theory_grid_periodic,
		aes(x=psi, y=c_amp, z=k_theory),
		col = "white", linewidth = 0.4, breaks = c(1, 2, 5, 10, 20, 50)) +
	scale_fill_viridis_c(option = "inferno", name = sprintf("k\n(capped\nat %d)", k_cap), limits = c(0, k_cap)) +
	theme_classic() +
	labs(x = expression(psi),
	     y = expression("Contact amplitude " * italic(A)),
	     title = sprintf("%s (R0 = %g)", pathogen, R0))


# ==============================================================================
# Periodic contacts: epidemic simulations
# ==============================================================================



# -----

# PERIODIC CONTACTS 
# - Define the contact function (e.g. `make_contact_fn_periodic()`)
#    * Returns R_0 * (1 - c_amp * cos(2*pi*t/c_per)), where R_0, c_amp, and 
#      c_per are input parameters 
# - Create a grid of psi/c_amp pairs. Start with psi in [0,1] and A in [0,1]. 
#    * Start with a 20x20 grid; we can refine later. 
# - For each psi/c_amp pair, calclulate k using the theoretical result
#    * Specifically, I think this should be: 
#      k = (2 / c_amp^2) * (1 + omega^2 / beta^2)^(psi * alpha)
#      where omega = 2*pi/c_per; where alpha and beta are the parameters of the
#      generation interval distribution (Gamma(alpha, beta)); and psi is the 
#      punctuation parameter. 
# - For each psi/c_amp pair, simulate secondary infections for 10000 index
#   cases 
#    * I think this can be done by simply running 
#      gen_inf_attempts_gamma_contacts() with the appropriate inputs. 
# - Using the simulated secondary infections, calculate k empirically for each 
#   psi/c_amp pair 
#    * I think this calculation is k = m_s^2 / (v_s - m_s), where 
#      m_s is the mean number of offspring and v_s is the variance in the number
#      of offspring. This equation only holds when v_s > m_s; otherwise, it 
#      should return Inf. 
# - Generate a heat map with the empirical k values for each psi/c_amp pair as 
#   the patches. Overlay contours reflecting the theoretical k values across 
#   psi/c_amp. 
# - Simulate 5,000 epidemics in 10,000 people for psi=0, psi=0.5, and psi=1, 
#   holding contact amplitude fixed at 0.5. Calculate the number of epidemics 
#   that don't establish, i.e., that fail to reach 500 infections. Report these 
#   failure proportions in a table. 
#    * I think this can be done using sim_stochastic_fast() with
#      gen_inf_attempts_gamma_contacts() and the appropriate inputs. 

# GAMMA CONTACTS 
# - Define the contact function (e.g. `make_contact_fn_gamma()`)
#    * This is genuinely new. See lines 87-94 of writeup/model.tex for a 
#      description of this contact model. Basically, for each person, contact 
#      levels change at arrivals of a Poisson process with rate lambda, so that
#      the amount of time spent at each level is exponentially distributed with
#      mean 1/lambda. The contact levels are drawn independently from a 
#      Gamma(k_c, k_c) distribution. As before, if a level is 
#      l ~ Gamma(k_c, k_c), then the function should output R_0*l for all the 
#      levels l. 
# - Create a grid of psi/c_amp pairs. Start with psi in [0,1] and k_c in 
#   [0.1,1000], using a log scale for k_c. 
#    * Start with a 20x20 grid; we can refine later. 
# - For each psi/c_amp pair, simulate secondary infections for 10000 index 
#   cases 
# - Using the simulated secondary infections, calculate k empirically for each 
#   psi/c_amp pair 
# - Generate a heat map with the empirical k values for each psi/c_amp pair as 
#   the patches. 
# - Simulate 5,000 epidemics in 10,000 people for psi=0, psi=0.5, and psi=1, 
#   holding contact amplitude fixed at 0.5. Calculate the number of epidemics 
#   that don't establish, i.e., that fail to reach 500 infections. Report these 
#   failure proportions in a table. 


