library(tidyverse)
source("code/utils.R")

# ==============================================================================
# Overdispersion from punctuated infectiousness x periodic contacts
#
# Key idea: when g(t) = 1 + epsilon*sin(2*pi*t/T), the realized R0 for
# individual i is nu_i = R0 * integral[b_i(tau) * g(t_i + tau) dtau].
# With punctuated profiles (small kappa), b_i is a narrow spike, so nu_i
# depends sensitively on whether the spike lands on a high- or low-contact
# phase. With smooth profiles (large kappa), b_i averages over the cycle
# and nu_i ≈ R0 for everyone.
# ==============================================================================

set.seed(42)

# ==============================================================================
# Parameters
# ==============================================================================

mu          <- 5       # mean generation time
R0          <- 2       # basic reproduction number
alpha_total <- 10      # population kernel shape
r           <- alpha_total / mu
T_period    <- 7       # contact cycle period (weekly)

# Parameter grid
kappas  <- c(0.5, 1, 2, 4, 6, 8, 9.5)
epsilons <- seq(0, 0.9, by = 0.1)

# Contact function: g(t) = 1 + epsilon * sin(2*pi*t / T)
make_contact_fn <- function(epsilon, T_period) {
	function(t) 1 + epsilon * sin(2 * pi * t / T_period)
}

# ==============================================================================
# Part A: Branching process offspring distribution
# ==============================================================================

cat("=== Part A: Branching process offspring distribution ===\n\n")

N_rep <- 10000

# Storage for results
results <- expand_grid(kappa = kappas, epsilon = epsilons) %>%
	mutate(mean_offspring = NA_real_,
	       var_offspring  = NA_real_,
	       var_mean_ratio = NA_real_,
	       nb_k           = NA_real_)

for (idx in seq_len(nrow(results))) {
	kappa   <- results$kappa[idx]
	epsilon <- results$epsilon[idx]

	alpha_shift <- alpha_total - kappa

	# Contact function for this epsilon
	g <- make_contact_fn(epsilon, T_period)

	# Precompute quadrature grid (midpoint rule avoids dgamma(0) = Inf for kappa < 1)
	u_max <- qgamma(0.9999, shape = kappa, rate = r)
	n_quad <- 500
	du <- u_max / n_quad
	u_grid <- seq(du / 2, u_max - du / 2, length.out = n_quad)
	f_kappa_vals <- dgamma(u_grid, shape = kappa, rate = r)

	offspring <- numeric(N_rep)
	for (i in seq_len(N_rep)) {
		# Random infection time (random phase)
		t_i <- runif(1, 0, T_period)
		# Onset shift for this individual
		s_i <- rgamma(1, shape = alpha_shift, rate = r)
		# Compute realized R0 via numerical integration
		# nu_i = R0 * integral[f_kappa(u) * g(t_i + s_i + u) du]
		g_vals <- g(t_i + s_i + u_grid)
		nu_i <- R0 * sum(f_kappa_vals * g_vals) * du

		# Draw offspring count
		offspring[i] <- rpois(1, nu_i)
	}

	results$mean_offspring[idx] <- mean(offspring)
	results$var_offspring[idx]  <- var(offspring)
	results$var_mean_ratio[idx] <- var(offspring) / mean(offspring)

	# Fit negative binomial dispersion parameter k
	# Method of moments: Var = mu + mu^2/k => k = mu^2 / (Var - mu)
	m <- mean(offspring)
	v <- var(offspring)
	if (v > m) {
		results$nb_k[idx] <- m^2 / (v - m)
	} else {
		results$nb_k[idx] <- Inf  # Poisson or underdispersed
	}

	if (idx %% 10 == 0) cat(sprintf("  %d / %d\n", idx, nrow(results)))
}

cat("\nDone computing offspring distributions.\n\n")

# Print summary table
cat("=== Summary: Variance/Mean ratio ===\n\n")
results %>%
	select(kappa, epsilon, mean_offspring, var_mean_ratio, nb_k) %>%
	filter(epsilon %in% c(0, 0.3, 0.6, 0.9)) %>%
	arrange(epsilon, kappa) %>%
	print(n = Inf)

# Verification: mean should be ≈ R0 everywhere
cat(sprintf("\n=== Mean offspring (should be ≈ %.1f everywhere) ===\n", R0))
cat(sprintf("  Range: [%.3f, %.3f]\n",
            min(results$mean_offspring), max(results$mean_offspring)))

# Verification: epsilon=0 should give var/mean ≈ 1 (Poisson)
cat("\n=== Var/Mean ratio at epsilon=0 (should be ≈ 1) ===\n")
results %>% filter(epsilon == 0) %>%
	select(kappa, var_mean_ratio) %>%
	print(n = Inf)

# ==============================================================================
# Part B: Visualizations
# ==============================================================================

cat("\n=== Part B: Generating figures ===\n\n")

# --- B1: Heatmap of variance/mean ratio ---

fig_heatmap <- results %>%
	ggplot(aes(x = kappa, y = epsilon, fill = var_mean_ratio)) +
	geom_tile() +
	scale_fill_viridis_c(option = "inferno", name = "Var / Mean") +
	geom_contour(aes(x = kappa, y = epsilon, z = var_mean_ratio),
	             col = "white", alpha = 0.5, inherit.aes = FALSE) +
	theme_classic() +
	labs(x = expression(kappa ~ "(punctuation: small = spiky, large = smooth)"),
	     y = expression(epsilon ~ "(contact amplitude)"),
	     title = "Offspring overdispersion: Var/Mean ratio",
	     subtitle = sprintf("R0 = %g, alpha_total = %g, weekly contact cycle (T = %g days)",
	                        R0, alpha_total, T_period))

ggsave("figures/fig_overdispersion_heatmap.pdf", fig_heatmap, width = 8, height = 6)
cat("Saved figures/fig_overdispersion_heatmap.pdf\n")

# --- B2: Line plot of variance vs kappa for several epsilon ---

fig_lines <- results %>%
	filter(epsilon %in% c(0, 0.3, 0.5, 0.7, 0.9)) %>%
	ggplot(aes(x = kappa, y = var_mean_ratio, col = factor(epsilon),
	           group = epsilon)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	geom_hline(yintercept = 1, linetype = "dashed", col = "grey50") +
	theme_classic() +
	scale_color_brewer(palette = "YlOrRd", name = expression(epsilon),
	                   direction = -1) +
	labs(x = expression(kappa),
	     y = "Var / Mean (offspring count)",
	     title = "Punctuation amplifies contact-driven overdispersion",
	     subtitle = "Dashed line = Poisson (Var/Mean = 1)")

ggsave("figures/fig_overdispersion_lines.pdf", fig_lines, width = 8, height = 5)
cat("Saved figures/fig_overdispersion_lines.pdf\n")

# --- B3: Example offspring histograms ---

example_pairs <- list(
	list(kappa = 9.5, epsilon = 0,   label = "Smooth, constant contacts"),
	list(kappa = 9.5, epsilon = 0.9, label = "Smooth, periodic contacts"),
	list(kappa = 0.5, epsilon = 0,   label = "Punctuated, constant contacts"),
	list(kappa = 0.5, epsilon = 0.9, label = "Punctuated, periodic contacts")
)

# Regenerate offspring distributions for example pairs
hist_df <- tibble()
for (ex in example_pairs) {
	kappa   <- ex$kappa
	epsilon <- ex$epsilon
	alpha_shift <- alpha_total - kappa
	g <- make_contact_fn(epsilon, T_period)

	u_max <- qgamma(0.9999, shape = kappa, rate = r)
	n_quad <- 500
	du <- u_max / n_quad
	u_grid <- seq(du / 2, u_max - du / 2, length.out = n_quad)
	f_kappa_vals <- dgamma(u_grid, shape = kappa, rate = r)

	offspring <- numeric(N_rep)
	for (i in seq_len(N_rep)) {
		t_i <- runif(1, 0, T_period)
		s_i <- rgamma(1, shape = alpha_shift, rate = r)
		g_vals <- g(t_i + s_i + u_grid)
		nu_i <- R0 * sum(f_kappa_vals * g_vals) * du
		offspring[i] <- rpois(1, nu_i)
	}

	hist_df <- bind_rows(hist_df, tibble(
		offspring = offspring,
		label = ex$label,
		kappa = kappa,
		epsilon = epsilon
	))
}

# Add Poisson reference
hist_df <- hist_df %>%
	mutate(label = factor(label, levels = sapply(example_pairs, `[[`, "label")))

# Compute Poisson reference for comparison
poisson_ref <- tibble(
	offspring = 0:max(hist_df$offspring),
	density = dpois(0:max(hist_df$offspring), lambda = R0)
)

fig_histograms <- hist_df %>%
	ggplot(aes(x = offspring)) +
	geom_histogram(aes(y = after_stat(density)), binwidth = 1,
	               fill = "steelblue", col = "white", alpha = 0.7) +
	geom_line(data = poisson_ref, aes(x = offspring, y = density),
	          col = "red", linewidth = 0.8, linetype = "dashed") +
	facet_wrap(~label, ncol = 2) +
	theme_classic() +
	labs(x = "Number of offspring",
	     y = "Density",
	     title = "Offspring distributions vs. Poisson(R0) reference",
	     subtitle = "Red dashed = Poisson(2); blue = simulated offspring distribution")

ggsave("figures/fig_overdispersion_histograms.pdf", fig_histograms, width = 10, height = 7)
cat("Saved figures/fig_overdispersion_histograms.pdf\n")

# --- B4: Heatmap of NB dispersion parameter k ---

fig_heatmap_k <- results %>%
	mutate(nb_k_capped = pmin(nb_k, 50)) %>%
	ggplot(aes(x = kappa, y = epsilon, fill = nb_k_capped)) +
	geom_tile() +
	scale_fill_viridis_c(option = "viridis", name = "NB k\n(capped at 50)",
	                     direction = -1) +
	theme_classic() +
	labs(x = expression(kappa ~ "(punctuation: small = spiky, large = smooth)"),
	     y = expression(epsilon ~ "(contact amplitude)"),
	     title = "Negative binomial dispersion parameter k",
	     subtitle = "Lower k = more overdispersed; k -> Inf = Poisson")

ggsave("figures/fig_overdispersion_nbk.pdf", fig_heatmap_k, width = 8, height = 6)
cat("Saved figures/fig_overdispersion_nbk.pdf\n")

# ==============================================================================
# Part C: Full epidemic simulations
# ==============================================================================

cat("\n=== Part C: Full epidemic simulations ===\n\n")

n_pop <- 1000
nsim  <- 200

# Selected scenarios
scenarios <- list(
	list(kappa = 9.5, epsilon = 0,   label = "Smooth + constant"),
	list(kappa = 0.5, epsilon = 0,   label = "Punctuated + constant"),
	list(kappa = 9.5, epsilon = 0.9, label = "Smooth + periodic"),
	list(kappa = 0.5, epsilon = 0.9, label = "Punctuated + periodic")
)

epi_df <- tibble()

for (sc in scenarios) {
	cat(sprintf("  Running: %s (kappa=%.1f, epsilon=%.1f)\n",
	            sc$label, sc$kappa, sc$epsilon))

	g <- make_contact_fn(sc$epsilon, T_period)
	g_max <- 1 + sc$epsilon  # supremum of 1 + epsilon*sin(...)

	gen <- make_profile_gamma_contacts(
		mu = mu, R0 = R0, alpha_total = alpha_total,
		kappa = sc$kappa, contact_fn = g, g_max = g_max
	)

	for (sim in seq_len(nsim)) {
		tinf <- sim_stochastic_fast(n = n_pop, gen_contacts = gen)
		n_infected <- sum(tinf < Inf)
		epi_df <- bind_rows(epi_df, tibble(
			sim = sim,
			kappa = sc$kappa,
			epsilon = sc$epsilon,
			label = sc$label,
			final_size = n_infected,
			tinf_sorted = list(sort(tinf[tinf < Inf]))
		))
	}
}

epi_df <- epi_df %>%
	mutate(label = factor(label, levels = sapply(scenarios, `[[`, "label")),
	       established = final_size >= 0.1 * n_pop)

# --- C1: Final size distributions ---

fig_final_size <- epi_df %>%
	ggplot(aes(x = final_size)) +
	geom_histogram(binwidth = 20, fill = "steelblue", col = "white", alpha = 0.7) +
	facet_wrap(~label, ncol = 2) +
	theme_classic() +
	labs(x = "Final epidemic size",
	     y = "Count",
	     title = "Final size distributions across scenarios",
	     subtitle = sprintf("N = %d, R0 = %g, %d simulations per scenario",
	                        n_pop, R0, nsim))

ggsave("figures/fig_epidemic_final_size.pdf", fig_final_size, width = 10, height = 7)
cat("Saved figures/fig_epidemic_final_size.pdf\n")

# --- C2: Summary statistics ---

cat("\n=== Final size summary (all epidemics) ===\n\n")
epi_df %>%
	group_by(label) %>%
	summarise(
		mean_fs = mean(final_size),
		sd_fs   = sd(final_size),
		median_fs = median(final_size),
		p_established = mean(established),
		.groups = "drop"
	) %>%
	print()

cat("\n=== Final size summary (established epidemics only) ===\n\n")
epi_df %>%
	filter(established) %>%
	group_by(label) %>%
	summarise(
		n = n(),
		mean_fs = mean(final_size),
		sd_fs   = sd(final_size),
		.groups = "drop"
	) %>%
	print()

# --- C3: Epidemic curves (established epidemics, first 10 per scenario) ---

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

# Keep first 10 established epidemics per scenario for plotting
if (nrow(curve_df) > 0) {
	plot_sims <- curve_df %>%
		group_by(label) %>%
		distinct(sim) %>%
		slice_head(n = 10)

	fig_curves <- curve_df %>%
		semi_join(plot_sims, by = c("label", "sim")) %>%
		ggplot(aes(x = tinf, y = cuminf, group = interaction(label, sim))) +
		geom_line(alpha = 0.5) +
		facet_wrap(~label, ncol = 2) +
		theme_classic() +
		labs(x = "Time (days)", y = "Cumulative infections",
		     title = "Example epidemic trajectories (established epidemics)")

	ggsave("figures/fig_epidemic_curves.pdf", fig_curves, width = 10, height = 7)
	cat("Saved figures/fig_epidemic_curves.pdf\n")
}

cat("\nDone.\n")
