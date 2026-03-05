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

# ==============================================================================
# 1. Parameters
# ==============================================================================

T <- 5              # Mean generation interval 
popshape <- 10      # Shape parameter for pop. generation interval distribution
r <- popshape/T     # Rate for generation interval distribution 
z_mean <- 2         # Contact function mean value (~R0)
z_per <- 7          # Contact function period 

# Parameter grid (kappa values and contact function amplitudes)
# kappa_vals <- seq(from=0.1, to=0.9, by=0.1)   # Punctuation
# z_amp_vals <- seq(from=0, to=1.2, by=0.1)       # Contact function amplitude

make_contact_fn <- function(z_mean, z_amp, z_per){
	stopifnot(z_amp <= z_mean)
	function(t) z_mean - z_amp*cos(2*pi*t/z_per)
}

# ==============================================================================
# 2. Overdispersion (k) vs kappa for three values of z_mean
# ==============================================================================

omega <- 2*pi/z_per # pre-compute period factor 
n_index <- 1000     # number of index cases per simulation replicate
n_reps  <- 10       # number of replicates per (kappa, z_mean) pair

z_mean_vals <- c(2, 5, 12)
kappa_vals <- seq(from=0.1, to=0.9, by=0.1)
kappa_vals_fine <- seq(0.01, 0.99, length.out = 200)

theory_df <- expand_grid(kappa = kappa_vals_fine, z_mean = z_mean_vals) %>%
	mutate(
		z_amp = z_mean - 1,
		eps = z_amp / z_mean,
		rho = (r / sqrt(r^2 + omega^2))^(kappa * popshape),
		var_nu = z_mean^2 * eps^2 * rho^2 / 2,
		k_theory = z_mean^2 / var_nu
	)

# --- Simulation points ---
sim_k_df <- expand_grid(
	kappa = kappa_vals,
	z_mean_sim = z_mean_vals,
	rep = 1:n_reps
) %>% mutate(k_sim = NA_real_)

for(idx in seq_len(nrow(sim_k_df))){
	kappa   <- sim_k_df$kappa[idx]
	z_mean    <- sim_k_df$z_mean_sim[idx]
	z_amp    <- z_mean - 1
	z  <- make_contact_fn(z_mean, z_amp, z_per)
	gfun  <- gen_inf_attempts_gamma_contacts(T=T, z=z, z_max=z_mean + z_amp,
	                                         popshape=popshape, kappa=kappa)
	tinfs <- z_per * runif(n_index)
	noffspring  <- sapply(lapply(tinfs, gfun), length)
	m_sim <- mean(noffspring)
	v_sim <- var(noffspring)
	sim_k_df$k_sim[idx] <- if(v_sim>m_sim){m_sim^2/(v_sim-m_sim)} else {Inf}
	if(idx %% 20 == 0) cat(sprintf("  %d / %d\n", idx, nrow(sim_k_df)))
}

# --- Plot ---
k_cap <- 50

fig_k_vs_kappa <- ggplot() +
	geom_line(data = theory_df,
		aes(x = kappa, y = k_theory, col = factor(z_mean)),
		linewidth = 1, alpha = 0.7) +
	geom_point(data = filter(sim_k_df, k_sim <= k_cap),
		aes(x = kappa, y = k_sim, col = factor(z_mean_sim)),
		alpha = 0.4, size = 0.8) +
	geom_point(data = filter(sim_k_df, k_sim > k_cap) %>%
			group_by(kappa, z_mean_sim) %>%
			mutate(idx = seq_along(k_sim),
			       row = (idx - 1) %/% 5,
			       col = (idx - 1) %% 5) %>%
			group_by(kappa, z_mean_sim, row) %>%
			mutate(n_in_row = n(),
			       kappa_spread = kappa + (col - (n_in_row - 1)/2) * 0.01,
			       k_plot = k_cap - row * 1.5) %>%
			ungroup(),
		aes(x = kappa_spread, y = k_plot, col = factor(z_mean_sim)),
		shape = 4, alpha = 0.5, size = 1.5) +
	scale_color_manual(
		values = c("2" = "black", "5" = "blue", "12" = "red"),
		name = expression(bar(z) ~ "(~ R"[0] * ")")) +
	scale_y_continuous(limits = c(0, k_cap), expand = expansion(mult = c(0.02, 0))) +
	coord_cartesian(clip = "off") +
	theme_classic() +
	labs(x = expression(kappa ~ "(0 = punctuated, 1 = smooth)"),
	     y = "Dispersion parameter k")

# ==============================================================================
# 3. Heatmaps of k over (kappa, z_amp) for z_mean = 2, 5, 12
# ==============================================================================

n_index_heat <- 10000
kappa_vals_heat   <- seq(0.05, 0.95, by = 0.05)

heatmap_list <- list()
sim_grid_all <- list()

for(z_mean in z_mean_vals){

	z_amp_max <- z_mean - 1
	z_amp_vals <- seq(0, z_amp_max, length.out = 20)

	# --- Simulation grid ---
	sim_grid <- expand_grid(kappa = kappa_vals_heat, z_amp = z_amp_vals) %>%
		mutate(k_sim = NA_real_)

	for(idx in seq_len(nrow(sim_grid))){
		kappa  <- sim_grid$kappa[idx]
		z_amp   <- sim_grid$z_amp[idx]
		z <- make_contact_fn(z_mean, z_amp, z_per)
		gfun <- gen_inf_attempts_gamma_contacts(T=T, z=z, z_max=z_mean + z_amp,
		                                        popshape=popshape, kappa=kappa)
		tinfs <- z_per * runif(n_index_heat)
		noffspring  <- sapply(lapply(tinfs, gfun), length)
		m_s <- mean(noffspring)
		v_s <- var(noffspring)
		sim_grid$k_sim[idx] <- if(v_s>m_s){m_s^2/(v_s-m_s)} else {Inf}
	}
	cat(sprintf("  Done simulations for z_mean = %g\n", z_mean))
	sim_grid_all[[as.character(z_mean)]] <- sim_grid

	sim_means <- sim_grid %>%
		mutate(z_mean = z_mean)

	# --- Theoretical contours ---
	theory_heat <- expand_grid(kappa = seq(0.01, 0.99, length.out = 100),
	                           z_amp = seq(0.01, z_amp_max, length.out = 100)) %>%
		mutate(
			eps = z_amp / z_mean,
			rho = (r / sqrt(r^2 + omega^2))^(kappa * popshape),
			var_nu = z_mean^2 * eps^2 * rho^2 / 2,
			k_theory = z_mean^2 / var_nu,
			z_mean = z_mean
		)

	# --- Plot ---
	heatmap_list[[as.character(z_mean)]] <- ggplot() +
		geom_tile(data = sim_means, aes(x = kappa, y = z_amp, fill = k_sim)) +
		geom_contour(data = theory_heat,
			aes(x = kappa, y = z_amp, z = k_theory),
			col = "white", alpha = 0.6, breaks = c(1, 2, 5, 10, 20, 50, 100)) +
		scale_fill_viridis_c(option = "inferno", name = "k",
		                     trans = "log10", limits = c(1, NA)) +
		theme_classic() +
		labs(x = expression(kappa),
		     y = "Contact amplitude",
		     title = bquote(bar(z) == .(z_mean)))
}

fig_heatmaps <- patchwork::wrap_plots(heatmap_list, nrow = 1) +
	patchwork::plot_annotation(
		title = "Offspring overdispersion (NB k) across punctuation and contact amplitude"
	)

# Same heatmaps with k capped at 50, linear scale
heatmap_capped_list <- list()
for(z_mean in z_mean_vals){
	z_amp_max <- z_mean - 1
	z_amp_vals <- seq(0, z_amp_max, length.out = 20)

	sim_means <- sim_grid_all[[as.character(z_mean)]] %>%
		mutate(k_capped = pmin(k_sim, 50), z_mean = z_mean)

	theory_heat <- expand_grid(kappa = seq(0.01, 0.99, length.out = 500),
	                           z_amp = seq(0.01, z_amp_max, length.out = 500)) %>%
		mutate(
			eps = z_amp / z_mean,
			rho = (r / sqrt(r^2 + omega^2))^(kappa * popshape),
			var_nu = z_mean^2 * eps^2 * rho^2 / 2,
			k_theory = pmin(z_mean^2 / var_nu, 50)
		)

	heatmap_capped_list[[as.character(z_mean)]] <- ggplot() +
		geom_tile(data = sim_means, aes(x = kappa, y = z_amp, fill = k_capped)) +
		geom_contour(data = theory_heat,
			aes(x = kappa, y = z_amp, z = k_theory),
			col = "black", linewidth = 1.0, breaks = c(1, 2, 5, 10, 20, 50)) +
		geom_contour(data = theory_heat,
			aes(x = kappa, y = z_amp, z = k_theory),
			col = "white", linewidth = 0.4, breaks = c(1, 2, 5, 10, 20, 50)) +
		scale_fill_viridis_c(option = "inferno", name = "k\n(capped\nat 50)",
		                     limits = c(0, 50)) +
		theme_classic() +
		labs(x = expression(kappa),
		     y = "Contact amplitude",
		     title = bquote(bar(z) == .(z_mean)))
}

fig_heatmaps_capped <- patchwork::wrap_plots(heatmap_capped_list, nrow = 1) +
	patchwork::plot_annotation(
		title = "Offspring overdispersion (NB k, capped at 50)"
	)

