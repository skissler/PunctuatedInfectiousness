library(tidyverse)
library(patchwork)
source("code/utils.R")
source("code/parameters.R")

# ==============================================================================
# Overdispersion from punctuated infectiousness x periodic contacts
#
# Corresponds to writeup sections:
#   - "The impact of the individual infectiousness profile on uncontrolled
#      epidemic dynamics
#
# - Computes offspring distribution for a bunch of one-step transmissions
# - Plots a heatmap of k as a function of punctuation (psi) and contact
#   amplitude
# - Plots a line plot of k as a function of punctuation (psi) for a few
#   contact amplitudes
# - Plots histograms of offspring distributions across punctuation (psi)
#   values with a Poisson reference
# - Plots epidemic trajectories and computes extinction probability, final
#   size, and time to 100 infected/time to peak
# ==============================================================================

# ==============================================================================
# 1. Global parameters (shared across pathogens)
# ==============================================================================

z_per <- 1          # Contact function period

make_contact_fn <- function(z_mean, z_amp, z_per){
	stopifnot(z_amp <= z_mean)
	function(t) z_mean - z_amp*cos(2*pi*t/z_per)
}

# Simulation parameters
n_index <- 1000     # number of index cases per simulation replicate
n_reps  <- 10       # number of replicates per (psi, z_mean) pair
n_index_heat <- 10000
psi_vals <- seq(from=0.01, to=0.2, by=0.01)
psi_vals_fine <- seq(0.01, 0.2, length.out = 200)
psi_vals_heat <- seq(0.01, 0.2, by = 0.01)
z_per_vals <- seq(0.5, 20, length.out = 20)
k_cap <- 50

# Epidemic simulation parameters
popsize <- 1000
nsim <- 200

# Storage for per-pathogen plots (assembled after loop)
line_plots <- list()
heatmap_list <- list()
heatmap_capped_list <- list()
heatmap_period_list <- list()
heatmap_period_capped_list <- list()

# ==============================================================================
# 2-4. Loop over pathogens
# ==============================================================================

for (pars in parslist) {

pathogen <- pars$pathogen
T        <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

z_mean   <- R0
z_amp    <- R0 - 1
omega    <- 2*pi/z_per

cat(sprintf("\n===== %s: T=%.2f, alpha=%.2f, beta=%.3f, R0=%.1f =====\n",
            pathogen, T, alpha, beta, R0))

# ==============================================================================
# 2. Overdispersion (k) vs psi
# ==============================================================================

theory_df <- tibble(psi = psi_vals_fine) %>%
	mutate(
		eps = z_amp / z_mean,
		rho = (beta / sqrt(beta^2 + omega^2))^(psi * alpha),
		var_nu = z_mean^2 * eps^2 * rho^2 / 2,
		k_theory = z_mean^2 / var_nu
	)

# --- Simulation points ---
sim_k_df <- expand_grid(
	psi = psi_vals,
	rep = 1:n_reps
) %>% mutate(k_sim = NA_real_)

for(idx in seq_len(nrow(sim_k_df))){

	psi_val <- sim_k_df$psi[idx]
	z       <- make_contact_fn(z_mean, z_amp, z_per)
	gfun    <- gen_inf_attempts_gamma_contacts(T=T, z=z, z_max=z_mean + z_amp,
	                                           alpha=alpha, psi=psi_val)

	tinfs       <- z_per * runif(n_index)
	noffspring  <- sapply(lapply(tinfs, gfun), length)

	m_sim <- mean(noffspring)
	v_sim <- var(noffspring)

	sim_k_df$k_sim[idx] <- if(v_sim>m_sim){m_sim^2/(v_sim-m_sim)} else {Inf}
	if(idx %% 20 == 0) cat(sprintf("  [%s] k vs psi: %d / %d\n", pathogen, idx, nrow(sim_k_df)))
}

# --- Plot ---
fig_k_vs_psi <- ggplot() +
	geom_line(data = theory_df,
		aes(x = psi, y = k_theory),
		linewidth = 1, alpha = 0.7) +
	geom_point(data = filter(sim_k_df, k_sim <= k_cap),
		aes(x = psi, y = k_sim),
		alpha = 0.4, size = 0.8) +
	geom_point(data = filter(sim_k_df, k_sim > k_cap) %>%
			group_by(psi) %>%
			mutate(idx = seq_along(k_sim),
			       row = (idx - 1) %/% 5,
			       col = (idx - 1) %% 5) %>%
			group_by(psi, row) %>%
			mutate(n_in_row = n(),
			       psi_spread = psi + (col - (n_in_row - 1)/2) * 0.01,
			       k_plot = k_cap - row * 1.5) %>%
			ungroup(),
		aes(x = psi_spread, y = k_plot),
		shape = 4, alpha = 0.5, size = 1.5) +
	scale_y_continuous(limits = c(0, k_cap), expand = expansion(mult = c(0.02, 0))) +
	coord_cartesian(clip = "off") +
	theme_classic() +
	labs(x = expression(psi ~ "(0 = punctuated, 1 = smooth)"),
	     y = "Dispersion parameter k",
	     title = sprintf("%s (R0 = %g)", pathogen, R0))

line_plots[[pathogen]] <- fig_k_vs_psi
save_fig(fig_k_vs_psi, sprintf("fig_overdispersion_lines_%s", pathogen))
cat(sprintf("Saved fig_overdispersion_lines_%s\n", pathogen))

# ==============================================================================
# 3. Heatmaps of k over (psi, z_amp)
# ==============================================================================

z_amp_max <- z_mean - 1
z_amp_vals <- seq(0, z_amp_max, length.out = 20)

# --- Simulation grid ---
sim_grid <- expand_grid(psi = psi_vals_heat, z_amp = z_amp_vals) %>%
	mutate(k_sim = NA_real_)

for(idx in seq_len(nrow(sim_grid))){
	psi_val <- sim_grid$psi[idx]
	za      <- sim_grid$z_amp[idx]
	z <- make_contact_fn(z_mean, za, z_per)
	gfun <- gen_inf_attempts_gamma_contacts(T=T, z=z, z_max=z_mean + za,
	                                        alpha=alpha, psi=psi_val)
	tinfs <- z_per * runif(n_index_heat)
	noffspring  <- sapply(lapply(tinfs, gfun), length)
	m_s <- mean(noffspring)
	v_s <- var(noffspring)
	sim_grid$k_sim[idx] <- if(v_s>m_s){m_s^2/(v_s-m_s)} else {Inf}
}
cat(sprintf("  [%s] Done heatmap simulations\n", pathogen))

# --- Theoretical contours ---
theory_heat <- expand_grid(psi = seq(0.01, 0.2, length.out = 100),
                           z_amp = seq(0.01, z_amp_max, length.out = 100)) %>%
	mutate(
		eps = z_amp / z_mean,
		rho = (beta / sqrt(beta^2 + omega^2))^(psi * alpha),
		var_nu = z_mean^2 * eps^2 * rho^2 / 2,
		k_theory = z_mean^2 / var_nu
	)

# --- Log-scale heatmap ---
heatmap_list[[pathogen]] <- ggplot() +
	geom_tile(data = sim_grid, aes(x = psi, y = z_amp, fill = k_sim)) +
	geom_contour(data = theory_heat,
		aes(x = psi, y = z_amp, z = k_theory),
		col = "white", alpha = 0.6, breaks = c(1, 2, 5, 10, 20, 50, 100)) +
	scale_fill_viridis_c(option = "inferno", name = "k",
	                     trans = "log10", limits = c(1, NA)) +
	theme_classic() +
	labs(x = expression(psi),
	     y = "Contact amplitude",
	     title = sprintf("%s (R0 = %g)", pathogen, R0))

save_fig(heatmap_list[[pathogen]], sprintf("fig_overdispersion_heatmap_%s", pathogen))
cat(sprintf("Saved fig_overdispersion_heatmap_%s\n", pathogen))

# --- Capped linear heatmap ---
sim_means_capped <- sim_grid %>%
	mutate(k_capped = pmin(k_sim, 50))

theory_heat_capped <- expand_grid(psi = seq(0.01, 0.2, length.out = 500),
                                  z_amp = seq(0.01, z_amp_max, length.out = 500)) %>%
	mutate(
		eps = z_amp / z_mean,
		rho = (beta / sqrt(beta^2 + omega^2))^(psi * alpha),
		var_nu = z_mean^2 * eps^2 * rho^2 / 2,
		k_theory = pmin(z_mean^2 / var_nu, 50)
	)

heatmap_capped_list[[pathogen]] <- ggplot() +
	geom_tile(data = sim_means_capped, aes(x = psi, y = z_amp, fill = k_capped)) +
	geom_contour(data = theory_heat_capped,
		aes(x = psi, y = z_amp, z = k_theory),
		col = "black", linewidth = 1.0, breaks = c(1, 2, 5, 10, 20, 50)) +
	geom_contour(data = theory_heat_capped,
		aes(x = psi, y = z_amp, z = k_theory),
		col = "white", linewidth = 0.4, breaks = c(1, 2, 5, 10, 20, 50)) +
	scale_fill_viridis_c(option = "inferno", name = "k\n(capped\nat 50)",
	                     limits = c(0, 50)) +
	theme_classic() +
	labs(x = expression(psi),
	     y = "Contact amplitude",
	     title = sprintf("%s (R0 = %g)", pathogen, R0))

save_fig(heatmap_capped_list[[pathogen]], sprintf("fig_overdispersion_heatmap_capped_%s", pathogen))
cat(sprintf("Saved fig_overdispersion_heatmap_capped_%s\n", pathogen))

# ==============================================================================
# 3b. Heatmaps of k over (psi, z_per)
# ==============================================================================

eps <- z_amp / z_mean

# --- Simulation grid ---
sim_grid_per <- expand_grid(psi = psi_vals_heat, z_per = z_per_vals) %>%
	mutate(k_sim = NA_real_)

for(idx in seq_len(nrow(sim_grid_per))){
	psi_val <- sim_grid_per$psi[idx]
	zp     <- sim_grid_per$z_per[idx]
	z      <- make_contact_fn(z_mean, z_amp, zp)
	gfun   <- gen_inf_attempts_gamma_contacts(T=T, z=z, z_max=z_mean + z_amp,
	                                          alpha=alpha, psi=psi_val)
	tinfs      <- zp * runif(n_index_heat)
	noffspring <- sapply(lapply(tinfs, gfun), length)
	m_s <- mean(noffspring)
	v_s <- var(noffspring)
	sim_grid_per$k_sim[idx] <- if(v_s > m_s){ m_s^2 / (v_s - m_s) } else { Inf }
}
cat(sprintf("  [%s] Done period heatmap simulations\n", pathogen))

# --- Theoretical contours ---
theory_heat_per <- expand_grid(psi = seq(0.01, 0.2, length.out = 100),
                               z_per = seq(0.5, 20, length.out = 100)) %>%
	mutate(
		omega_p  = 2 * pi / z_per,
		rho      = (beta / sqrt(beta^2 + omega_p^2))^(psi * alpha),
		var_nu   = z_mean^2 * eps^2 * rho^2 / 2,
		k_theory = z_mean^2 / var_nu
	)

# --- Log-scale heatmap ---
heatmap_period_list[[pathogen]] <- ggplot() +
	geom_tile(data = sim_grid_per, aes(x = psi, y = z_per, fill = k_sim)) +
	geom_contour(data = theory_heat_per,
		aes(x = psi, y = z_per, z = k_theory),
		col = "white", alpha = 0.6, breaks = c(1, 2, 5, 10, 20, 50, 100)) +
	scale_fill_viridis_c(option = "viridis", name = "k",
	                     trans = "log10", limits = c(1, NA)) +
	theme_classic() +
	labs(x = expression(psi),
	     y = "Contact period (days)",
	     title = sprintf("%s (R0 = %g)", pathogen, R0))

save_fig(heatmap_period_list[[pathogen]], sprintf("fig_overdispersion_heatmap_period_%s", pathogen))
cat(sprintf("Saved fig_overdispersion_heatmap_period_%s\n", pathogen))

# --- Capped linear heatmap ---
theory_heat_per_capped <- expand_grid(psi = seq(0.01, 0.2, length.out = 500),
                                      z_per = seq(0.5, 20, length.out = 500)) %>%
	mutate(
		omega_p  = 2 * pi / z_per,
		rho      = (beta / sqrt(beta^2 + omega_p^2))^(psi * alpha),
		var_nu   = z_mean^2 * eps^2 * rho^2 / 2,
		k_theory = pmin(z_mean^2 / var_nu, 50)
	)

sim_per_capped <- sim_grid_per %>%
	mutate(k_capped = pmin(k_sim, 50))

heatmap_period_capped_list[[pathogen]] <- ggplot() +
	geom_tile(data = sim_per_capped, aes(x = psi, y = z_per, fill = k_capped)) +
	geom_contour(data = theory_heat_per_capped,
		aes(x = psi, y = z_per, z = k_theory),
		col = "black", linewidth = 1.0, breaks = c(1, 2, 5, 10, 20, 50)) +
	geom_contour(data = theory_heat_per_capped,
		aes(x = psi, y = z_per, z = k_theory),
		col = "white", linewidth = 0.4, breaks = c(1, 2, 5, 10, 20, 50)) +
	scale_fill_viridis_c(option = "inferno", name = "k\n(capped\nat 50)",
	                     limits = c(0, 50)) +
	theme_classic() +
	labs(x = expression(psi),
	     y = "Contact period (days)",
	     title = sprintf("%s (R0 = %g)", pathogen, R0))

save_fig(heatmap_period_capped_list[[pathogen]], sprintf("fig_overdispersion_heatmap_period_capped_%s", pathogen))
cat(sprintf("Saved fig_overdispersion_heatmap_period_capped_%s\n", pathogen))

# ==============================================================================
# 4. Full epidemic simulations
# ==============================================================================

scenarios <- list(
	list(psi=0.25, z_amp=0,        label="Smooth + constant"),
	list(psi=0.05, z_amp=0,        label="Punctuated + constant"),
	list(psi=0.25, z_amp=R0 - 1,   label="Smooth + periodic"),
	list(psi=0.05, z_amp=R0 - 1,   label="Punctuated + periodic")
)

epi_df <- tibble()

for(sc in scenarios){

	z <- make_contact_fn(z_mean=z_mean, z_amp=sc$z_amp, z_per=z_per)
	z_max <- z_mean + sc$z_amp

	gfun <- gen_inf_attempts_gamma_contacts(T=T, z=z, z_max=z_max,
	                                        alpha=alpha, psi=sc$psi)

	for(sim in seq_len(nsim)){
		tinf <- sim_stochastic_fast(n=popsize, gen_inf_attempts=gfun)
		n_infected <- sum(tinf<Inf)
		epi_df <- bind_rows(epi_df, tibble(
			sim=sim,
			psi=sc$psi,
			z_mean=z_mean,
			z_amp=sc$z_amp,
			label=sc$label,
			final_size=n_infected,
			tinf_sorted=list(sort(tinf[tinf<Inf]))
		))
		if(sim %% 20 == 0) cat(sprintf("  [%s] %s: %d / %d\n", pathogen, sc$label, sim, nsim))
	}

}

epi_df <- epi_df %>%
	mutate(established=as.integer(final_size>=0.1*popsize))

fig_final_size <- epi_df %>%
	filter(established==1) %>%
	ggplot(aes(x = final_size)) +
	geom_histogram(binwidth = 1, fill = "steelblue", col = "white", alpha = 0.7) +
	facet_wrap(~label, ncol = 2) +
	theme_classic() +
	labs(x = "Final epidemic size",
	     y = "Count",
	     title = sprintf("Final size distributions (%s, R0 = %g)", pathogen, R0),
	     subtitle = sprintf("N = %d, %d simulations per scenario", popsize, nsim))

save_fig(fig_final_size, sprintf("fig_epidemic_final_size_%s", pathogen))
cat(sprintf("Saved fig_epidemic_final_size_%s\n", pathogen))

fs_table <- epi_df %>%
	group_by(label) %>%
	summarise(
		mean_fs = mean(final_size),
		sd_fs   = sd(final_size),
		median_fs = median(final_size),
		p_established = mean(established),
		.groups = "drop"
	)

fs_table_established <- epi_df %>%
	filter(established==1) %>%
	group_by(label) %>%
	summarise(
		n = n(),
		mean_fs = mean(final_size),
		sd_fs   = sd(final_size),
		.groups = "drop"
	) %>%
	print()


curve_df <- tibble()
for (idx in seq_len(nrow(epi_df))) {
	row <- epi_df[idx, ]
	if (!row$established) next
	ts <- row$tinf_sorted[[1]]
	curve_df <- bind_rows(curve_df, tibble(
		sim = row$sim,
		label = row$label,
		tinf = ts,
		cuminf = seq_along(ts)
	))
}

# Keep all established epidemics for plotting
if (nrow(curve_df) > 0) {
	plot_sims <- curve_df %>%
		group_by(label) %>%
		distinct(sim)

	fig_curves <- curve_df %>%
		semi_join(plot_sims, by = c("label", "sim")) %>%
		ggplot(aes(x = tinf, y = cuminf, group = interaction(label, sim))) +
		geom_line(alpha = 0.2, linewidth=0.2) +
		facet_wrap(~label, ncol = 2) +
		theme_classic() +
		labs(x = "Time (days)", y = "Cumulative infections",
		     title = sprintf("Epidemic trajectories (%s, R0 = %g)", pathogen, R0))

	save_fig(fig_curves, sprintf("fig_epidemic_curves_%s", pathogen))
	cat(sprintf("Saved fig_epidemic_curves_%s\n", pathogen))
}

} # end pathogen loop

# ==============================================================================
# 5. Composite figures across pathogens (patchwork)
# ==============================================================================

fig_lines_composite <- patchwork::wrap_plots(line_plots, nrow = 1) +
	patchwork::plot_annotation(
		title = "Offspring overdispersion (NB k) vs punctuation"
	)
save_fig(fig_lines_composite, "fig_overdispersion_lines", width = 14, height = 5)
cat("Saved fig_overdispersion_lines\n")

fig_heatmaps <- patchwork::wrap_plots(heatmap_list, nrow = 1) +
	patchwork::plot_annotation(
		title = "Offspring overdispersion (NB k) across punctuation and contact amplitude"
	)
save_fig(fig_heatmaps, "fig_overdispersion_heatmap", width = 14, height = 5)
cat("Saved fig_overdispersion_heatmap\n")

fig_heatmaps_capped <- patchwork::wrap_plots(heatmap_capped_list, nrow = 1) +
	patchwork::plot_annotation(
		title = "Offspring overdispersion (NB k, capped at 50)"
	)
save_fig(fig_heatmaps_capped, "fig_overdispersion_heatmap_capped", width = 14, height = 5)
cat("Saved fig_overdispersion_heatmap_capped\n")

fig_heatmaps_period <- patchwork::wrap_plots(heatmap_period_list, nrow = 1) +
	patchwork::plot_annotation(
		title = "Offspring overdispersion (NB k) across punctuation and contact period"
	)
save_fig(fig_heatmaps_period, "fig_overdispersion_heatmap_period", width = 14, height = 5)
cat("Saved fig_overdispersion_heatmap_period\n")

fig_heatmaps_period_capped <- patchwork::wrap_plots(heatmap_period_capped_list, nrow = 1) +
	patchwork::plot_annotation(
		title = "Offspring overdispersion (NB k, capped at 50) across punctuation and contact period"
	)
save_fig(fig_heatmaps_period_capped, "fig_overdispersion_heatmap_period_capped", width = 14, height = 5)
cat("Saved fig_overdispersion_heatmap_period_capped\n")
