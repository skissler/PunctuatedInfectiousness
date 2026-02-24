library(tidyverse)

# ==============================================================================
# Visualize the interpolated individual infectiousness profile a_i(tau)
# for different values of the punctuation parameter kappa.
#
# Construction (from notes/overview.md):
#   Draw tau_i* ~ Abar(tau)  (the generation interval distribution)
#   a_i(tau) = R0 * Abar(tau) * h(tau; tau_i*, kappa) / E[h(tau; tau_i*, kappa)]
#
# where h is a Gaussian kernel: h(tau; tau_i*, kappa) = exp(-kappa * (tau - tau_i*)^2)
# and E[h] = integral of Abar(s) * h(tau; s, kappa) ds  (a Gaussian smoothing of Abar).
#
# Limits:
#   kappa -> 0:   h ~ const, so a_i(tau) -> A(tau) for all i  (smooth/all-equal)
#   kappa -> Inf: h -> delta,  so a_i(tau) -> R0 * delta(tau - tau_i*)  (spike)
# ==============================================================================

# Parameters
e_dur <- 2
i_dur <- 3
R0 <- 2
nu <- 1 / e_dur    # latency rate
gam <- 1 / i_dur   # recovery rate

# Normalized generation interval density: Abar(tau) = A(tau) / R0
Abar <- function(tau) {
	nu * gam / (nu - gam) * (exp(-gam * tau) - exp(-nu * tau))
}

# Numerical grid
tau_max <- 25
n_grid <- 2000
tau_grid <- seq(0, tau_max, length.out = n_grid)
dtau <- tau_grid[2] - tau_grid[1]
Abar_grid <- Abar(tau_grid)

# Compute E[h(tau; tau*, kappa)] on the grid for a given kappa.
# This is the Gaussian convolution of Abar evaluated at each tau.
compute_Eh <- function(kappa) {
	# E[h(tau_j; tau*, kappa)] = sum_k Abar(tau_k) * exp(-kappa*(tau_j - tau_k)^2) * dtau
	# Vectorized via outer product (n_grid x n_grid matrix; fine for n_grid ~ 2000)
	diffs <- outer(tau_grid, tau_grid, "-")       # tau_j - tau_k
	H <- exp(-kappa * diffs^2)                     # h(tau_j; tau_k, kappa)
	as.numeric(H %*% (Abar_grid * dtau))           # integrate over tau_k
}

# Compute a_i(tau) for a given tau_i* and precomputed Eh
compute_ai <- function(tau_star, kappa, Eh) {
	h_vals <- exp(-kappa * (tau_grid - tau_star)^2)
	R0 * Abar_grid * h_vals / Eh
}

# ==============================================================================
# Generate profiles for a range of kappa values
# ==============================================================================

set.seed(42)
kappas <- c(0.01, 0.1, 0.5, 1, 2, 5, 20)
n_individuals <- 8

# Draw tau_i* values (same individuals across all kappa for comparison)
tau_stars <- rexp(n_individuals, 1 / e_dur) + rexp(n_individuals, 1 / i_dur)

# Build a data frame of all profiles
profile_df <- tibble()
integral_df <- tibble()

for(kappa in kappas) {
	Eh <- compute_Eh(kappa)

	for(j in seq_along(tau_stars)) {
		ai <- compute_ai(tau_stars[j], kappa, Eh)
		ai_integral <- sum(ai) * dtau

		profile_df <- bind_rows(profile_df, tibble(
			tau = tau_grid,
			ai = ai,
			kappa = kappa,
			individual = j,
			tau_star = tau_stars[j]
		))
		integral_df <- bind_rows(integral_df, tibble(
			kappa = kappa,
			individual = j,
			tau_star = tau_stars[j],
			integral = ai_integral
		))
	}
}

profile_df <- profile_df %>%
	mutate(kappa_label = factor(
		sprintf("kappa = %s", kappa),
		levels = sprintf("kappa = %s", sort(unique(kappa)))
	))

# ==============================================================================
# Plot 1: Individual profiles across kappa values
# ==============================================================================

Atau_ref <- expand_grid(
	tau = tau_grid,
	kappa_label = levels(profile_df$kappa_label)
) %>% mutate(Atau = R0 * Abar(tau))

fig_profiles <- profile_df %>%
	ggplot(aes(x = tau, y = ai, group = individual, col = factor(individual))) +
		geom_line(alpha = 0.7) +
		geom_line(aes(x = tau, y = Atau), col = "black", linewidth = 0.8,
		          linetype = "dashed", inherit.aes = FALSE, data = Atau_ref) +
		facet_wrap(~kappa_label, scales = "free_y") +
		theme_classic() +
		labs(x = expression(tau ~ "(days since infection)"),
		     y = expression(a[i](tau)),
		     title = "Individual infectiousness profiles across kappa",
		     subtitle = "Dashed black = population-level A(tau)") +
		theme(legend.position = "none")

ggsave("figures/fig_profiles.pdf", fig_profiles, width = 14, height = 8)
cat("Saved figures/fig_profiles.pdf\n")

# ==============================================================================
# Plot 2: Population average vs. A(tau)
# ==============================================================================

pop_avg_df <- profile_df %>%
	group_by(kappa, kappa_label, tau) %>%
	summarise(mean_ai = mean(ai), .groups = "drop")

fig_pop_avg <- pop_avg_df %>%
	ggplot(aes(x = tau)) +
		geom_line(aes(y = mean_ai), col = "red", linewidth = 0.8) +
		geom_line(aes(y = R0 * Abar(tau)), col = "black",
		          linewidth = 0.8, linetype = "dashed") +
		facet_wrap(~kappa_label) +
		theme_classic() +
		labs(x = expression(tau ~ "(days since infection)"),
		     y = expression(a[i](tau)),
		     title = "Population average (red) vs. A(tau) (dashed black)",
		     subtitle = sprintf("Average over %d individuals", n_individuals))

ggsave("figures/fig_pop_avg.pdf", fig_pop_avg, width = 14, height = 8)
cat("Saved figures/fig_pop_avg.pdf\n")

# ==============================================================================
# Plot 3: Individual total infectiousness (integral of a_i)
# ==============================================================================

fig_integrals <- integral_df %>%
	ggplot(aes(x = factor(kappa), y = integral)) +
		geom_hline(yintercept = R0, linetype = "dashed", col = "blue") +
		geom_jitter(width = 0.15, alpha = 0.6) +
		theme_classic() +
		labs(x = expression(kappa),
		     y = expression(integral(a[i](tau) * dtau)),
		     title = "Total individual infectiousness across kappa",
		     subtitle = sprintf("Dashed line = R0 = %g", R0))

ggsave("figures/fig_integrals.pdf", fig_integrals, width = 8, height = 5)
cat("Saved figures/fig_integrals.pdf\n")

# Print summary
cat("\n=== Total infectiousness per individual (should be R0 =", R0, ") ===\n\n")
integral_df %>%
	group_by(kappa) %>%
	summarise(mean = mean(integral), sd = sd(integral),
	          min = min(integral), max = max(integral), .groups = "drop") %>%
	print()

# NOTE: This construction guarantees E[a_i(tau)] = A(tau) (population average)
# but does NOT guarantee integral(a_i) = R0 for each individual. The individual
# integrals are close to R0 (within ~5%) and converge as kappa -> 0 or kappa -> Inf.
# If exact conservation is needed, each a_i can be renormalized to integrate to R0.
