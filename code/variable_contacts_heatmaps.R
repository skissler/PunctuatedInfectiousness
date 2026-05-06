# ==============================================================================
# Import
# ==============================================================================

library(tidyverse)
source("code/utils.R")
source("code/global_parameters.R")
source("code/parameters.R")

# c_per             <- 1
c_per             <- 7
c_amp_vals_sim    <- seq(0, 1, length.out=25)
# psi_vals_sim      <- seq(0, 0.5, length.out=25)
psi_vals_sim      <- seq(0, 1, length.out=25)
c_amp_vals_theory <- seq(0, 1, length.out=100)
# psi_vals_theory   <- seq(0, 0.5, length.out=100)
psi_vals_theory   <- seq(0, 1, length.out=100)

# lambda_gp        <- 2     # Poisson switching rate (switches/day)
lambda_gp        <- 1     # Poisson switching rate (switches/day)
k_c_vals_sim     <- exp(seq(log(0.1), log(1000), length.out=25))
k_c_vals_theory  <- exp(seq(log(0.1), log(1000), length.out=100))

n_index <- 10000
k_cap   <- 50

# Define contact modulator: c_i(t) = 1 - c_amp*cos(2*pi*t/c_per)
# Returns R0*c_i(t), for use with gen_inf_attempts_gamma_contacts
make_contact_fn_periodic <- function(R0, c_amp, c_per) {
	stopifnot(c_amp >= 0, c_amp <= 1)
	function(t) R0*(1 - c_amp*cos(2*pi*t/c_per))
}

# ==============================================================================
# Start pathogen loop
# ==============================================================================

for (pars in parslist) {

pathogen <- pars$pathogen
Tgen     <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

cat(sprintf("\n===== %s: T=%.2f, alpha=%.2f, R0=%.1f =====\n",
            pathogen, Tgen, alpha, R0))

# ==============================================================================
# Periodic contacts: heatmaps
# ==============================================================================

# Initialize a data frame for simulated k-values
sim_grid_periodic <- expand_grid(psi=psi_vals_sim, c_amp=c_amp_vals_sim) %>%
	mutate(k_sim = NA_real_)

# Generate the simulated k-values
for(idx in 1:nrow(sim_grid_periodic)){
	psi_val   <- sim_grid_periodic$psi[idx]
	c_amp_val <- sim_grid_periodic$c_amp[idx]
	z     <- make_contact_fn_periodic(R0, c_amp_val, c_per)
	z_max <- R0*(1 + c_amp_val)
	gfun  <- gen_inf_attempts_gamma_contacts(T=Tgen, z=z, z_max=z_max,
	                                         alpha=alpha, psi=psi_val)
	tinfs      <- c_per * runif(n_index)
	noffspring <- sapply(lapply(tinfs, gfun), length)
	m_s <- mean(noffspring)
	v_s <- var(noffspring)
	sim_grid_periodic$k_sim[idx] <- if(v_s > m_s){ m_s^2 / (v_s - m_s) } else { Inf }
	if(idx %% 50 == 0) cat(sprintf("  [%s] heatmap: %d / %d\n", pathogen, idx, nrow(sim_grid_periodic)))
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

save_fig(fig_heatmap_periodic, sprintf("fig_heatmap_periodic_%s", pathogen), width=5, height=5)
cat(sprintf("  Saved fig_heatmap_periodic_%s\n", pathogen))

# ==============================================================================
# Periodic contacts: epidemic simulations
# ==============================================================================

c_amp_epi <- 0.5
psi_epi   <- c(0, 0.5, 1)

epi_results <- vector("list", nsim_small * length(psi_epi))
epi_idx     <- 0L

for (psi_val in psi_epi) {
	z     <- make_contact_fn_periodic(R0, c_amp_epi, c_per)
	z_max <- R0 * (1 + c_amp_epi)
	gfun  <- gen_inf_attempts_gamma_contacts(T=Tgen, z=z, z_max=z_max,
	                                         alpha=alpha, psi=psi_val)
	for (sim in 1:nsim_small) {
		tinf <- sim_stochastic_fast(n=popsize, gen_inf_attempts=gfun,
		                            maxinf=establishment_threshold)
		n_infected <- sum(is.finite(tinf))

		epi_idx <- epi_idx + 1L
		epi_results[[epi_idx]] <- tibble(
			sim = sim, psi = psi_val, n_infected = n_infected,
			established = as.integer(n_infected >= establishment_threshold))

		if (sim %% 500 == 0) cat(sprintf("  [%s] psi=%.1f: sim %d/%d\n",
		                                  pathogen, psi_val, sim, nsim_small))
	}
}

epi_df <- bind_rows(epi_results)

pest_table <- epi_df %>%
	group_by(psi) %>%
	summarise(p_extinct = 1 - mean(established), .groups = "drop")

cat(sprintf("  %s: P(extinction) by psi (c_amp = %g):\n", pathogen, c_amp_epi))
print(pest_table)


# ==============================================================================
# Gamma/Poisson contacts: heatmaps
# ==============================================================================

# Initialize a data frame for simulated k-values
sim_grid_gp <- expand_grid(psi=psi_vals_sim, k_c=k_c_vals_sim) %>%
	mutate(k_sim = NA_real_)

# Generate the simulated k-values
for(idx in 1:nrow(sim_grid_gp)){
	psi_val <- sim_grid_gp$psi[idx]
	k_c_val <- sim_grid_gp$k_c[idx]
	gfun <- gen_inf_attempts_gammapoisson_contacts(Tgen, R0, alpha, psi_val, k_c_val, lambda_gp)
	noffspring <- sapply(rep(0, n_index), function(t) length(gfun(t)))
	m_s <- mean(noffspring)
	v_s <- var(noffspring)
	sim_grid_gp$k_sim[idx] <- if(v_s > m_s){ m_s^2 / (v_s - m_s) } else { Inf }
	if(idx %% 50 == 0) cat(sprintf("  [%s] GP heatmap: %d / %d\n", pathogen, idx, nrow(sim_grid_gp)))
}

# Add capped k column
sim_grid_gp <- sim_grid_gp %>%
	mutate(k_capped = pmin(k_sim, k_cap))

# Generate theoretical k-values for contour overlay
theory_grid_gp <- expand_grid(psi=psi_vals_theory, k_c=k_c_vals_theory) %>%
	mutate(k_theory = mapply(k_theory_gammapoisson, k_c, psi, alpha, lambda_gp, beta)) %>%
	mutate(k_capped = pmin(k_theory, k_cap))

# Generate the heat map with theoretical contours
fig_heatmap_gp <- ggplot() +
	geom_tile(data = sim_grid_gp, aes(x=psi, y=k_c, fill=k_capped)) +
	geom_contour(data = theory_grid_gp,
		aes(x=psi, y=k_c, z=k_theory),
		col = "black", linewidth = 1.0, breaks = c(1, 2, 5, 10, 20, 50)) +
	geom_contour(data = theory_grid_gp,
		aes(x=psi, y=k_c, z=k_theory),
		col = "white", linewidth = 0.4, breaks = c(1, 2, 5, 10, 20, 50)) +
	scale_y_log10() +
	scale_fill_viridis_c(option = "inferno", name = sprintf("k\n(capped\nat %d)", k_cap), limits = c(0, k_cap)) +
	theme_classic() +
	labs(x = expression(psi),
	     y = expression("Contact shape " * italic(k)[c]),
	     title = sprintf("%s: Gamma/Poisson contacts (R0 = %g, λ = %g)", pathogen, R0, lambda_gp))

save_fig(fig_heatmap_gp, sprintf("fig_heatmap_gammapoisson_%s", pathogen), width=5, height=5)
cat(sprintf("  Saved fig_heatmap_gammapoisson_%s\n", pathogen))

# ==============================================================================
# Gamma/Poisson contacts: epidemic simulations
# ==============================================================================

k_c_epi <- 1     # moderate heterogeneity in contact levels
psi_epi_gp <- c(0, 0.5, 1)

epi_results_gp <- vector("list", nsim_small * length(psi_epi_gp))
epi_idx_gp     <- 0L

for (psi_val in psi_epi_gp) {
	gfun <- gen_inf_attempts_gammapoisson_contacts(Tgen, R0, alpha, psi_val, k_c_epi, lambda_gp)
	for (sim in 1:nsim_small) {
		tinf <- sim_stochastic_fast(n=popsize, gen_inf_attempts=gfun,
		                            maxinf=establishment_threshold)
		n_infected <- sum(is.finite(tinf))

		epi_idx_gp <- epi_idx_gp + 1L
		epi_results_gp[[epi_idx_gp]] <- tibble(
			sim = sim, psi = psi_val, n_infected = n_infected,
			established = as.integer(n_infected >= establishment_threshold))

		if (sim %% 500 == 0) cat(sprintf("  [%s] GP psi=%.1f: sim %d/%d\n",
		                                  pathogen, psi_val, sim, nsim_small))
	}
}

epi_df_gp <- bind_rows(epi_results_gp)

pest_table_gp <- epi_df_gp %>%
	group_by(psi) %>%
	summarise(p_extinct = 1 - mean(established), .groups = "drop")

cat(sprintf("  %s: P(extinction) by psi (Gamma/Poisson, k_c = %g, lambda = %g):\n",
            pathogen, k_c_epi, lambda_gp))
print(pest_table_gp)



} # end pathogen loop



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


# THEORY
# I wasn't able to pin down exact page numbers from web searches — the book's content isn't indexed
#   at that granularity online. But based on the 4th edition's chapter structure (852 pages total),
#   here's what I can tell you with confidence:

#   Papoulis & Pillai, 4th ed. (2002) — the relevant chapters are:

#   - Chapter 9: General Concepts (~p. 363) — introduces stochastic processes, stationarity,
#   autocorrelation functions
#   - Chapter 11: Spectral Representation (~p. 465) — defines the power spectral density, proves the
#   Wiener–Khinchin theorem (PSD ↔ autocovariance as Fourier pair), and then derives the filtering
#   result: if WSS input X(t) passes through a linear system with transfer function H(ω), the output
#   PSD is S_YY(ω) = |H(ω)|² S_XX(ω). The output variance formula Var[Y] = ∫ |H(ω)|² S_XX(ω) dω is an
#   immediate corollary (set lag = 0 in the output autocovariance).

#   The key equation you want is typically labeled something like Eq. (11-49) or nearby (different
#   printings vary), in the section titled "System Response — Power Spectral Density" within Chapter
#   11. The result is stated as:

#   $$S_{YY}(\omega) = |H(\omega)|^2 S_{XX}(\omega)$$

#   and the variance follows by integration: Var[Y] = R_YY(0) = ∫ S_YY(ω) dω = ∫ |H|² S_XX dω.

#   I'd recommend looking at the index entry for "Linear systems, response to random input" or "Power
#   spectrum, output" — those should point you directly to the right pages. The treatment is
#   self-contained and quite readable; Papoulis builds it up from the convolution integral with clear
#   physical interpretation.

#   For a second, more compact treatment, the
#   https://ocw.mit.edu/courses/6-011-introduction-to-communication-control-and-signal-processing-sprin
#   g-2010/8075041184d566103ce7c3f69afc5e75_MIT6_011S10_chap10.pdf cover the same material in about 20
#   pages and are freely available.

#   Sources:
#   -
#   https://books.google.com/books/about/Probability_Random_Variables_and_Stochas.html?id=YYwQAQAAIAAJ
#   - https://searchworks.stanford.edu/view/4708590
#   - https://www.dsprelated.com/books/1633.php
#   - https://ocw.mit.edu/courses/6-011-introduction-to-communication-control-and-signal-processing-spr
#   ing-2010/8075041184d566103ce7c3f69afc5e75_MIT6_011S10_chap10.pdf