library(tidyverse)
library(patchwork)
source("code/utils.R")
source("code/parameters.R")

# ==============================================================================
# Overdispersion from punctuated infectiousness x periodic contacts
#
# Uses the contact modulator parameterization from the paper:
#   c(t) = 1 - A*cos(omega*t)
# where A in [0,1] is the amplitude and c(t) has mean 1. The actual contact
# rate is R0*c(t), so contacts range from R0*(1-A) to R0*(1+A).
#
# Closed-form dispersion parameter (Eq. 6 in the paper):
#   k = (2/A^2) * (1 + omega^2/beta^2)^(psi*alpha)
#
# - Computes offspring distribution for one-step transmissions
# - Plots k vs psi for a fixed amplitude (A=0.5) across all pathogens
# - Plots per-pathogen k vs psi with theory and simulation
# - Plots heatmaps of k over (psi, A) and (psi, period)
# - Plots epidemic trajectories under periodic contacts
# ==============================================================================

# ==============================================================================
# 1. Global parameters
# ==============================================================================

z_per   <- 1          # Default contact period (1 day = diurnal)
A_default <- 0.5      # Default amplitude for line plots

# Contact modulator: c(t) = 1 - A*cos(2*pi*t/P), mean = 1
# Returns z(t) = R0 * c(t), the actual contact rate for use with
# gen_inf_attempts_gamma_contacts
make_contact_fn <- function(R0, A, z_per) {
	stopifnot(A >= 0, A <= 1)
	function(t) R0 * (1 - A * cos(2 * pi * t / z_per))
}

# Simulation parameters
n_index      <- 1000   # index cases per replicate for line plots
n_reps       <- 10     # replicates per (psi, A) pair
n_index_heat <- 10000  # index cases for heatmap cells
psi_vals      <- seq(0.01, 0.20, by = 0.01)        # simulation grid (narrow — k explodes fast)
psi_vals_fine <- seq(0, 1, length.out = 500)        # theory curves (fine, for log-scale plots)
psi_vals_heat      <- seq(0.01, 0.20, by = 0.01)   # heatmap grid (narrow range)
psi_vals_heat_wide <- seq(0.01, 0.50, by = 0.025)  # heatmap grid (wide range)
z_per_vals    <- seq(0.5, 20, length.out = 20)      # period grid for heatmaps
k_cap <- 50

# Epidemic simulation parameters
popsize <- 1000
nsim    <- 200

# Storage for per-pathogen plots (assembled after loop)
line_plots                      <- list()
heatmap_list                    <- list()
heatmap_capped_list             <- list()
heatmap_list_wide               <- list()
heatmap_capped_list_wide        <- list()
heatmap_period_list             <- list()
heatmap_period_capped_list      <- list()

# Storage for cross-pathogen k vs psi data
cross_pathogen_theory <- list()

# ==============================================================================
# 2-4. Loop over pathogens
# ==============================================================================

for (pars in parslist) {

pathogen <- pars$pathogen
T        <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

omega    <- 2 * pi / z_per

cat(sprintf("\n===== %s: T=%.2f, alpha=%.2f, beta=%.3f, R0=%.1f =====\n",
            pathogen, T, alpha, beta, R0))

# ==============================================================================
# 2. Overdispersion (k) vs psi — single amplitude A_default
# ==============================================================================

# --- Theory: k = (2/A^2) * (1 + omega^2/beta^2)^(psi*alpha) ---
theory_df <- tibble(psi = psi_vals_fine) %>%
	mutate(
		k_theory = (2 / A_default^2) * (1 + omega^2 / beta^2)^(psi * alpha)
	)

# Store for cross-pathogen figure
cross_pathogen_theory[[pathogen]] <- theory_df %>%
	mutate(pathogen = pathogen)

# --- Simulation points ---
sim_k_df <- expand_grid(
	psi = psi_vals,
	rep = 1:n_reps
) %>% mutate(k_sim = NA_real_)

for(idx in seq_len(nrow(sim_k_df))){

	psi_val <- sim_k_df$psi[idx]
	z       <- make_contact_fn(R0, A_default, z_per)
	z_max   <- R0 * (1 + A_default)
	gfun    <- gen_inf_attempts_gamma_contacts(T=T, z=z, z_max=z_max,
	                                           alpha=alpha, psi=psi_val)

	tinfs       <- z_per * runif(n_index)
	noffspring  <- sapply(lapply(tinfs, gfun), length)

	m_sim <- mean(noffspring)
	v_sim <- var(noffspring)

	sim_k_df$k_sim[idx] <- if(v_sim > m_sim){ m_sim^2 / (v_sim - m_sim) } else { Inf }
	if(idx %% 50 == 0) cat(sprintf("  [%s] k vs psi: %d / %d\n", pathogen, idx, nrow(sim_k_df)))
}

# --- Per-pathogen line plot ---
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
			       psi_spread = psi + (col - (n_in_row - 1)/2) * 0.02,
			       k_plot = k_cap - row * 1.5) %>%
			ungroup(),
		aes(x = psi_spread, y = k_plot),
		shape = 4, alpha = 0.5, size = 1.5) +
	scale_y_continuous(limits = c(0, k_cap), expand = expansion(mult = c(0.02, 0))) +
	coord_cartesian(clip = "off") +
	theme_classic() +
	labs(x = expression(psi),
	     y = "Dispersion parameter k",
	     title = sprintf("%s (R0 = %g, A = %g)", pathogen, R0, A_default))

line_plots[[pathogen]] <- fig_k_vs_psi
save_fig(fig_k_vs_psi, sprintf("fig_overdispersion_lines_%s", pathogen))
cat(sprintf("Saved fig_overdispersion_lines_%s\n", pathogen))

# ==============================================================================
# 3. Heatmaps of k over (psi, A)
# ==============================================================================

A_vals <- seq(0.02, 0.98, length.out = 20)

# --- Simulation grid ---
sim_grid <- expand_grid(psi = psi_vals_heat, A = A_vals) %>%
	mutate(k_sim = NA_real_)

for(idx in seq_len(nrow(sim_grid))){
	psi_val <- sim_grid$psi[idx]
	A_val   <- sim_grid$A[idx]
	z     <- make_contact_fn(R0, A_val, z_per)
	z_max <- R0 * (1 + A_val)
	gfun  <- gen_inf_attempts_gamma_contacts(T=T, z=z, z_max=z_max,
	                                         alpha=alpha, psi=psi_val)
	tinfs      <- z_per * runif(n_index_heat)
	noffspring <- sapply(lapply(tinfs, gfun), length)
	m_s <- mean(noffspring)
	v_s <- var(noffspring)
	sim_grid$k_sim[idx] <- if(v_s > m_s){ m_s^2 / (v_s - m_s) } else { Inf }
}
cat(sprintf("  [%s] Done heatmap simulations\n", pathogen))

# --- Theoretical contours ---
theory_heat <- expand_grid(psi = seq(0.01, 0.20, length.out = 100),
                           A = seq(0.01, 0.99, length.out = 100)) %>%
	mutate(
		k_theory = (2 / A^2) * (1 + omega^2 / beta^2)^(psi * alpha)
	)

# --- Log-scale heatmap ---
heatmap_list[[pathogen]] <- ggplot() +
	geom_tile(data = sim_grid, aes(x = psi, y = A, fill = k_sim)) +
	geom_contour(data = theory_heat,
		aes(x = psi, y = A, z = k_theory),
		col = "white", alpha = 0.6, breaks = c(1, 2, 5, 10, 20, 50, 100)) +
	scale_fill_viridis_c(option = "inferno", name = "k",
	                     trans = "log10", limits = c(1, NA)) +
	theme_classic() +
	labs(x = expression(psi),
	     y = expression("Contact amplitude " * italic(A)),
	     title = sprintf("%s (R0 = %g)", pathogen, R0))

save_fig(heatmap_list[[pathogen]], sprintf("fig_overdispersion_heatmap_%s", pathogen))
cat(sprintf("Saved fig_overdispersion_heatmap_%s\n", pathogen))

# --- Capped linear heatmap ---
sim_means_capped <- sim_grid %>%
	mutate(k_capped = pmin(k_sim, k_cap))

theory_heat_capped <- expand_grid(psi = seq(0.01, 0.20, length.out = 200),
                                  A = seq(0.01, 0.99, length.out = 200)) %>%
	mutate(
		k_theory = pmin((2 / A^2) * (1 + omega^2 / beta^2)^(psi * alpha), k_cap)
	)

heatmap_capped_list[[pathogen]] <- ggplot() +
	geom_tile(data = sim_means_capped, aes(x = psi, y = A, fill = k_capped)) +
	geom_contour(data = theory_heat_capped,
		aes(x = psi, y = A, z = k_theory),
		col = "black", linewidth = 1.0, breaks = c(1, 2, 5, 10, 20, 50)) +
	geom_contour(data = theory_heat_capped,
		aes(x = psi, y = A, z = k_theory),
		col = "white", linewidth = 0.4, breaks = c(1, 2, 5, 10, 20, 50)) +
	scale_fill_viridis_c(option = "inferno", name = sprintf("k\n(capped\nat %d)", k_cap),
	                     limits = c(0, k_cap)) +
	theme_classic() +
	labs(x = expression(psi),
	     y = expression("Contact amplitude " * italic(A)),
	     title = sprintf("%s (R0 = %g)", pathogen, R0))

save_fig(heatmap_capped_list[[pathogen]], sprintf("fig_overdispersion_heatmap_capped_%s", pathogen))
cat(sprintf("Saved fig_overdispersion_heatmap_capped_%s\n", pathogen))

# ==============================================================================
# 3a. Heatmaps of k over (psi, A) — WIDE psi range [0, 0.5]
# ==============================================================================

# --- Simulation grid (wider psi range) ---
sim_grid_wide <- expand_grid(psi = psi_vals_heat_wide, A = A_vals) %>%
	mutate(k_sim = NA_real_)

for(idx in seq_len(nrow(sim_grid_wide))){
	psi_val <- sim_grid_wide$psi[idx]
	A_val   <- sim_grid_wide$A[idx]
	z     <- make_contact_fn(R0, A_val, z_per)
	z_max <- R0 * (1 + A_val)
	gfun  <- gen_inf_attempts_gamma_contacts(T=T, z=z, z_max=z_max,
	                                         alpha=alpha, psi=psi_val)
	tinfs      <- z_per * runif(n_index_heat)
	noffspring <- sapply(lapply(tinfs, gfun), length)
	m_s <- mean(noffspring)
	v_s <- var(noffspring)
	sim_grid_wide$k_sim[idx] <- if(v_s > m_s){ m_s^2 / (v_s - m_s) } else { Inf }
}
cat(sprintf("  [%s] Done wide-range heatmap simulations\n", pathogen))

# --- Theoretical contours (wider psi range) ---
theory_heat_wide <- expand_grid(psi = seq(0.01, 0.50, length.out = 100),
                                A = seq(0.01, 0.99, length.out = 100)) %>%
	mutate(
		k_theory = (2 / A^2) * (1 + omega^2 / beta^2)^(psi * alpha)
	)

# --- Log-scale heatmap (wide) ---
heatmap_list_wide[[pathogen]] <- ggplot() +
	geom_tile(data = sim_grid_wide, aes(x = psi, y = A, fill = k_sim)) +
	geom_contour(data = theory_heat_wide,
		aes(x = psi, y = A, z = k_theory),
		col = "white", alpha = 0.6, breaks = c(1, 2, 5, 10, 20, 50, 100, 1000, 1e4, 1e6)) +
	scale_fill_viridis_c(option = "inferno", name = "k",
	                     trans = "log10", limits = c(1, NA)) +
	theme_classic() +
	labs(x = expression(psi),
	     y = expression("Contact amplitude " * italic(A)),
	     title = sprintf("%s (R0 = %g)", pathogen, R0))

save_fig(heatmap_list_wide[[pathogen]], sprintf("fig_overdispersion_heatmap_wide_%s", pathogen))
cat(sprintf("Saved fig_overdispersion_heatmap_wide_%s\n", pathogen))

# --- Capped linear heatmap (wide) ---
sim_means_capped_wide <- sim_grid_wide %>%
	mutate(k_capped = pmin(k_sim, k_cap))

theory_heat_capped_wide <- expand_grid(psi = seq(0.01, 0.50, length.out = 200),
                                       A = seq(0.01, 0.99, length.out = 200)) %>%
	mutate(
		k_theory = pmin((2 / A^2) * (1 + omega^2 / beta^2)^(psi * alpha), k_cap)
	)

heatmap_capped_list_wide[[pathogen]] <- ggplot() +
	geom_tile(data = sim_means_capped_wide, aes(x = psi, y = A, fill = k_capped)) +
	geom_contour(data = theory_heat_capped_wide,
		aes(x = psi, y = A, z = k_theory),
		col = "black", linewidth = 1.0, breaks = c(1, 2, 5, 10, 20, 50)) +
	geom_contour(data = theory_heat_capped_wide,
		aes(x = psi, y = A, z = k_theory),
		col = "white", linewidth = 0.4, breaks = c(1, 2, 5, 10, 20, 50)) +
	scale_fill_viridis_c(option = "inferno", name = sprintf("k\n(capped\nat %d)", k_cap),
	                     limits = c(0, k_cap)) +
	theme_classic() +
	labs(x = expression(psi),
	     y = expression("Contact amplitude " * italic(A)),
	     title = sprintf("%s (R0 = %g)", pathogen, R0))

save_fig(heatmap_capped_list_wide[[pathogen]], sprintf("fig_overdispersion_heatmap_capped_wide_%s", pathogen))
cat(sprintf("Saved fig_overdispersion_heatmap_capped_wide_%s\n", pathogen))

# ==============================================================================
# 3b. Heatmaps of k over (psi, period)
# ==============================================================================

# --- Simulation grid ---
sim_grid_per <- expand_grid(psi = psi_vals_heat, z_per = z_per_vals) %>%
	mutate(k_sim = NA_real_)

for(idx in seq_len(nrow(sim_grid_per))){
	psi_val <- sim_grid_per$psi[idx]
	zp      <- sim_grid_per$z_per[idx]
	z       <- make_contact_fn(R0, A_default, zp)
	z_max   <- R0 * (1 + A_default)
	gfun    <- gen_inf_attempts_gamma_contacts(T=T, z=z, z_max=z_max,
	                                           alpha=alpha, psi=psi_val)
	tinfs      <- zp * runif(n_index_heat)
	noffspring <- sapply(lapply(tinfs, gfun), length)
	m_s <- mean(noffspring)
	v_s <- var(noffspring)
	sim_grid_per$k_sim[idx] <- if(v_s > m_s){ m_s^2 / (v_s - m_s) } else { Inf }
}
cat(sprintf("  [%s] Done period heatmap simulations\n", pathogen))

# --- Theoretical contours ---
theory_heat_per <- expand_grid(psi = seq(0.01, 0.20, length.out = 100),
                               z_per = seq(0.5, 20, length.out = 100)) %>%
	mutate(
		omega_p  = 2 * pi / z_per,
		k_theory = (2 / A_default^2) * (1 + omega_p^2 / beta^2)^(psi * alpha)
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
	     title = sprintf("%s (R0 = %g, A = %g)", pathogen, R0, A_default))

save_fig(heatmap_period_list[[pathogen]], sprintf("fig_overdispersion_heatmap_period_%s", pathogen))
cat(sprintf("Saved fig_overdispersion_heatmap_period_%s\n", pathogen))

# --- Capped linear heatmap ---
theory_heat_per_capped <- expand_grid(psi = seq(0.01, 0.20, length.out = 200),
                                      z_per = seq(0.5, 20, length.out = 200)) %>%
	mutate(
		omega_p  = 2 * pi / z_per,
		k_theory = pmin((2 / A_default^2) * (1 + omega_p^2 / beta^2)^(psi * alpha), k_cap)
	)

sim_per_capped <- sim_grid_per %>%
	mutate(k_capped = pmin(k_sim, k_cap))

heatmap_period_capped_list[[pathogen]] <- ggplot() +
	geom_tile(data = sim_per_capped, aes(x = psi, y = z_per, fill = k_capped)) +
	geom_contour(data = theory_heat_per_capped,
		aes(x = psi, y = z_per, z = k_theory),
		col = "black", linewidth = 1.0, breaks = c(1, 2, 5, 10, 20, 50)) +
	geom_contour(data = theory_heat_per_capped,
		aes(x = psi, y = z_per, z = k_theory),
		col = "white", linewidth = 0.4, breaks = c(1, 2, 5, 10, 20, 50)) +
	scale_fill_viridis_c(option = "inferno", name = sprintf("k\n(capped\nat %d)", k_cap),
	                     limits = c(0, k_cap)) +
	theme_classic() +
	labs(x = expression(psi),
	     y = "Contact period (days)",
	     title = sprintf("%s (R0 = %g, A = %g)", pathogen, R0, A_default))

save_fig(heatmap_period_capped_list[[pathogen]], sprintf("fig_overdispersion_heatmap_period_capped_%s", pathogen))
cat(sprintf("Saved fig_overdispersion_heatmap_period_capped_%s\n", pathogen))

# ==============================================================================
# 4. Full epidemic simulations
# ==============================================================================

scenarios <- list(
	list(psi=0.25, A=0,          label="Smooth + constant"),
	list(psi=0.05, A=0,          label="Punctuated + constant"),
	list(psi=0.25, A=A_default,  label="Smooth + periodic"),
	list(psi=0.05, A=A_default,  label="Punctuated + periodic")
)

epi_df <- tibble()

for(sc in scenarios){

	if(sc$A > 0){
		z     <- make_contact_fn(R0, sc$A, z_per)
		z_max <- R0 * (1 + sc$A)
		gfun  <- gen_inf_attempts_gamma_contacts(T=T, z=z, z_max=z_max,
		                                         alpha=alpha, psi=sc$psi)
	} else {
		gfun <- gen_inf_attempts_gamma(T, R0, alpha, sc$psi)
	}

	for(sim in seq_len(nsim)){
		tinf <- sim_stochastic_fast(n=popsize, gen_inf_attempts=gfun)
		n_infected <- sum(is.finite(tinf))
		epi_df <- bind_rows(epi_df, tibble(
			sim=sim,
			psi=sc$psi,
			A=sc$A,
			label=sc$label,
			final_size=n_infected,
			tinf_sorted=list(sort(tinf[is.finite(tinf)]))
		))
		if(sim %% 20 == 0) cat(sprintf("  [%s] %s: %d / %d\n", pathogen, sc$label, sim, nsim))
	}

}

epi_df <- epi_df %>%
	mutate(established=as.integer(final_size >= 0.1 * popsize))

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
		title = sprintf("Offspring overdispersion (NB k) vs punctuation (A = %g, period = %g day)", A_default, z_per)
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
		title = sprintf("Offspring overdispersion (NB k, capped at %d)", k_cap)
	)
save_fig(fig_heatmaps_capped, "fig_overdispersion_heatmap_capped", width = 14, height = 5)
cat("Saved fig_overdispersion_heatmap_capped\n")

fig_heatmaps_wide <- patchwork::wrap_plots(heatmap_list_wide, nrow = 1) +
	patchwork::plot_annotation(
		title = expression("Offspring overdispersion (NB k) across punctuation and contact amplitude (" * psi * " in [0, 0.5])")
	)
save_fig(fig_heatmaps_wide, "fig_overdispersion_heatmap_wide", width = 14, height = 5)
cat("Saved fig_overdispersion_heatmap_wide\n")

fig_heatmaps_capped_wide <- patchwork::wrap_plots(heatmap_capped_list_wide, nrow = 1) +
	patchwork::plot_annotation(
		title = sprintf("Offspring overdispersion (NB k, capped at %d) — psi in [0, 0.5]", k_cap)
	)
save_fig(fig_heatmaps_capped_wide, "fig_overdispersion_heatmap_capped_wide", width = 14, height = 5)
cat("Saved fig_overdispersion_heatmap_capped_wide\n")

fig_heatmaps_period <- patchwork::wrap_plots(heatmap_period_list, nrow = 1) +
	patchwork::plot_annotation(
		title = sprintf("Offspring overdispersion (NB k) across punctuation and contact period (A = %g)", A_default)
	)
save_fig(fig_heatmaps_period, "fig_overdispersion_heatmap_period", width = 14, height = 5)
cat("Saved fig_overdispersion_heatmap_period\n")

fig_heatmaps_period_capped <- patchwork::wrap_plots(heatmap_period_capped_list, nrow = 1) +
	patchwork::plot_annotation(
		title = sprintf("Offspring overdispersion (NB k, capped at %d) across punctuation and contact period (A = %g)", k_cap, A_default)
	)
save_fig(fig_heatmaps_period_capped, "fig_overdispersion_heatmap_period_capped", width = 14, height = 5)
cat("Saved fig_overdispersion_heatmap_period_capped\n")

# ==============================================================================
# 6. Cross-pathogen k vs psi figure (A = 0.5, all pathogens on one plot)
# ==============================================================================

cross_df <- bind_rows(cross_pathogen_theory) %>%
	mutate(pathogen = factor(pathogen, levels = c("influenza", "omicron", "measles"),
	                         labels = c("Influenza", "Omicron", "Measles")))

# Log-scale version
fig_k_cross_log <- ggplot(cross_df, aes(x = psi, y = k_theory, col = pathogen)) +
	geom_line(linewidth = 1) +
	geom_hline(yintercept = 2 / A_default^2, lty = "dashed", col = "grey50", linewidth = 0.4) +
	annotate("text", x = 0.97, y = 0.7 * 2 / A_default^2,
	         label = sprintf("k = 2/A\u00b2 = %g", 2 / A_default^2),
	         hjust = 1, size = 3, col = "grey40") +
	scale_y_log10(breaks = 10^(0:6), labels = scales::label_comma()) +
	coord_cartesian(ylim = c(5, 1e6)) +
	scale_color_manual(values = c("Influenza" = "#E41A1C", "Omicron" = "#377EB8", "Measles" = "#4DAF4A")) +
	labs(x = expression(psi),
	     y = expression("Dispersion parameter " * italic(k) * " (log scale)"),
	     col = NULL,
	     title = expression("Overdispersion from periodic contacts (" * italic(A) * " = 0.5, period = 1 day)")) +
	theme_classic(base_size = 13) +
	theme(legend.position.inside = c(0.15, 0.85),
	      legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

save_fig(fig_k_cross_log, "fig_k_vs_psi_periodic", width = 7, height = 5)
cat("Saved fig_k_vs_psi_periodic\n")

# Linear-scale version (capped)
fig_k_cross_linear <- cross_df %>%
	mutate(k_capped = pmin(k_theory, k_cap)) %>%
	ggplot(aes(x = psi, y = k_capped, col = pathogen)) +
	geom_line(linewidth = 1) +
	geom_hline(yintercept = 2 / A_default^2, lty = "dashed", col = "grey50", linewidth = 0.4) +
	annotate("text", x = 0.97, y = 0.85 * 2 / A_default^2,
	         label = sprintf("k = 2/A\u00b2 = %g", 2 / A_default^2),
	         hjust = 1, size = 3, col = "grey40") +
	scale_y_continuous(breaks = seq(0, k_cap, by = 10), limits = c(0, k_cap)) +
	scale_color_manual(values = c("Influenza" = "#E41A1C", "Omicron" = "#377EB8", "Measles" = "#4DAF4A")) +
	labs(x = expression(psi),
	     y = expression("Dispersion parameter " * italic(k)),
	     col = NULL,
	     title = expression("Overdispersion from periodic contacts (" * italic(A) * " = 0.5, period = 1 day)")) +
	theme_classic(base_size = 13) +
	theme(legend.position.inside = c(0.85, 0.35),
	      legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

save_fig(fig_k_cross_linear, "fig_k_vs_psi_periodic_linear", width = 7, height = 5)
cat("Saved fig_k_vs_psi_periodic_linear\n")
