library(tidyverse)

# ==============================================================================
# Visualize the shifted Gamma individual infectiousness profile
#
# Construction (from make_profile_gamma in utils.R):
#   infection_attempt_time = s_i + epsilon_j
#   s_i     ~ Gamma(alpha_total - kappa, r)   individual onset shift
#   eps_j   ~ Gamma(kappa, r)                  jitter (same for all individuals)
#   r       = alpha_total / mu
#
# Since Gamma + Gamma = Gamma (same rate), the marginal distribution of
# infection attempt times is Gamma(alpha_total, r), i.e. the population kernel is
# A(tau) = R0 * dgamma(tau; alpha_total, r), INVARIANT across kappa.
#
# Each individual profile is:
#   a_i(tau) = R0 * dgamma(tau - s_i; kappa, r)    for tau > s_i
#
# All profiles have the EXACT same shape, height, and width — just shifted.
#
# Limits:
#   kappa -> 0:            narrow spikes at random s_i (punctuated)
#   kappa -> alpha_total:  broad bumps ≈ A(tau)        (smooth)
# ==============================================================================

# Parameters
mu <- 5
R0 <- 2
alpha_total <- 10
r <- alpha_total / mu

cat(sprintf("Population kernel: Gamma(%g, %g), mean = %.1f, sd = %.2f\n",
            alpha_total, r, mu, mu / sqrt(alpha_total)))

# Numerical grid
tau_max <- 20
n_grid <- 2000
tau_grid <- seq(0, tau_max, length.out = n_grid)

# Population kernel (same for ALL kappa values)
Atau <- R0 * dgamma(tau_grid, shape = alpha_total, rate = r)

# ==============================================================================
# Generate profiles for a range of kappa values
# ==============================================================================

set.seed(42)
kappas <- c(0.5, 2, 4, 6, 8, 9.5)
n_individuals <- 8

# Build data frames
profile_df <- tibble()

for (kappa in kappas) {

	alpha_shift <- alpha_total - kappa

	# Draw onset shifts for this kappa
	s_stars <- rgamma(n_individuals, shape = alpha_shift, rate = r)

	for (j in seq_along(s_stars)) {
		# Individual profile: Gamma(kappa, r) density shifted by s_i
		ai <- R0 * dgamma(tau_grid - s_stars[j], shape = kappa, rate = r)
		ai[tau_grid < s_stars[j]] <- 0  # zero before onset

		profile_df <- bind_rows(profile_df, tibble(
			tau = tau_grid,
			ai = ai,
			kappa = kappa,
			individual = j,
			s_star = s_stars[j]
		))
	}
}

profile_df <- profile_df %>%
	mutate(kappa_label = factor(
		sprintf("kappa == %s", kappa),
		levels = sprintf("kappa == %s", sort(unique(kappa)))
	))

# ==============================================================================
# Plot 1: Individual profiles across kappa values
# ==============================================================================

Atau_ref <- expand_grid(
	tau = tau_grid,
	kappa_label = levels(profile_df$kappa_label)
) %>% mutate(Atau = R0 * dgamma(tau, shape = alpha_total, rate = r))

fig_profiles <- profile_df %>%
	ggplot(aes(x = tau)) +
		geom_line(aes(y = ai, group = individual, col = factor(individual)),
		          alpha = 0.7) +
		geom_line(aes(x = tau, y = Atau), col = "black", linewidth = 0.8,
		          linetype = "dashed", inherit.aes = FALSE, data = Atau_ref) +
		facet_wrap(~kappa_label, scales = "free_y",
		           labeller = label_parsed) +
		theme_classic() +
		labs(x = expression(tau ~ "(days since infection)"),
		     y = expression(a[i](tau)),
		     title = sprintf("Shifted Gamma profiles (alpha_total = %g, mu = %g)", alpha_total, mu),
		     subtitle = "Dashed black = A(tau), same in every panel; colored = individual profiles (identical shape, shifted)") +
		theme(legend.position = "none")

ggsave("figures/fig_profiles_gamma.pdf", fig_profiles, width = 14, height = 8)
cat("Saved figures/fig_profiles_gamma.pdf\n")

# ==============================================================================
# Plot 2: Population average vs. A(tau)
# ==============================================================================

pop_avg_df <- profile_df %>%
	group_by(kappa, kappa_label, tau) %>%
	summarise(mean_ai = mean(ai), .groups = "drop")

fig_pop_avg <- pop_avg_df %>%
	ggplot(aes(x = tau)) +
		geom_line(aes(y = mean_ai), col = "red", linewidth = 0.8) +
		geom_line(aes(y = R0 * dgamma(tau, shape = alpha_total, rate = r)),
		          col = "black", linewidth = 0.8, linetype = "dashed") +
		facet_wrap(~kappa_label, labeller = label_parsed) +
		theme_classic() +
		labs(x = expression(tau ~ "(days since infection)"),
		     y = expression(a[i](tau)),
		     title = "Sample mean of individual profiles (red) vs. A(tau) (dashed black)",
		     subtitle = sprintf("Average over %d individuals -- A(tau) is the SAME target in every panel", n_individuals))

ggsave("figures/fig_pop_avg_gamma.pdf", fig_pop_avg, width = 14, height = 8)
cat("Saved figures/fig_pop_avg_gamma.pdf\n")

# ==============================================================================
# Plot 3: Individual total infectiousness (integral of a_i)
# ==============================================================================

dtau <- tau_grid[2] - tau_grid[1]
integral_df <- profile_df %>%
	group_by(kappa, individual, s_star) %>%
	summarise(integral = sum(ai) * dtau, .groups = "drop")

fig_integrals <- integral_df %>%
	ggplot(aes(x = factor(kappa), y = integral)) +
		geom_hline(yintercept = R0, linetype = "dashed", col = "blue") +
		geom_jitter(width = 0.15, alpha = 0.6) +
		theme_classic() +
		labs(x = expression(kappa),
		     y = expression(integral(a[i](tau) * dtau)),
		     title = "Total individual infectiousness across kappa",
		     subtitle = sprintf("Dashed line = R0 = %g", R0))

ggsave("figures/fig_integrals_gamma.pdf", fig_integrals, width = 8, height = 5)
cat("Saved figures/fig_integrals_gamma.pdf\n")

# Print summary
cat(sprintf("\n=== Total infectiousness per individual (should be R0 = %g) ===\n\n", R0))
integral_df %>%
	group_by(kappa) %>%
	summarise(mean = mean(integral), sd = sd(integral),
	          min = min(integral), max = max(integral), .groups = "drop") %>%
	print()
