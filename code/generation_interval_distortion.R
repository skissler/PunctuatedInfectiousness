library(tidyverse)
source("code/utils.R")

# ==============================================================================
# Generation interval distortion under test-and-isolate
#
# Test-based screening with isolation truncates individual infectiousness
# profiles at the detection time. This reshapes the generation interval
# distribution of successful transmissions differently depending on kappa:
#
#   - Spike (small kappa): all-or-nothing removal; generation interval shape
#     among surviving transmissions is approximately unchanged.
#   - Smooth (large kappa): right-tail truncation; generation interval shifts
#     shorter, partially offsetting the R_eff reduction via Euler-Lotka.
#
# See notes/findings.md, Section 14.
# ==============================================================================

set.seed(2024)

# ==============================================================================
# Section 1: Parameters
# ==============================================================================

mu <- 5
R0 <- 2
alpha_total <- 10
r <- alpha_total / mu

kappas <- c(0.5, 1, 2, 4, 6, 8, 9.5)

# Test and screening parameters
sensitivity <- 0.85
fp_rate <- 0.02
window_before <- 3
window_after <- 3    # symmetric window for this analysis
screen_interval <- 3 # test every 3 days

n_ind <- 100000  # individuals per kappa value

cat("=== Generation interval distortion under test-and-isolate ===\n")
cat(sprintf("  mu=%g, R0=%g, alpha_total=%g, r=%g\n", mu, R0, alpha_total, r))
cat(sprintf("  screen_interval=%d, sensitivity=%.2f, window=[-%g, +%g]\n",
            screen_interval, sensitivity, window_before, window_after))
cat(sprintf("  n_ind=%d\n\n", n_ind))

# ==============================================================================
# Section 2: Simulate test-and-isolate, collect generation intervals
# ==============================================================================

get_peak_offset <- function(kappa, r) max(0, (kappa - 1) / r)

# For each kappa, simulate n_ind individuals:
#   - Draw individual profile (s_i, peak, attempt times)
#   - Apply screening with random phase
#   - Record which attempts survive (occur before detection)
#   - Collect generation intervals of surviving attempts

all_gi_natural <- list()   # generation intervals without intervention
all_gi_effective <- list() # generation intervals under test-and-isolate
summary_results <- list()

for (ki in seq_along(kappas)) {
	kappa <- kappas[ki]
	alpha_shift <- alpha_total - kappa
	peak_offset <- get_peak_offset(kappa, r)

	gi_natural <- numeric(0)
	gi_effective <- numeric(0)
	n_attempts_total <- 0
	n_attempts_surviving <- 0

	for (i in 1:n_ind) {
		# Draw individual
		s_i <- rgamma(1, shape = alpha_shift, rate = r)
		nattempts <- rpois(1, R0)
		if (nattempts == 0L) next

		# Attempt times (relative to infection at t=0)
		attempt_times <- s_i + sort(rgamma(nattempts, shape = kappa, rate = r))
		t_peak <- s_i + peak_offset

		# Record natural generation intervals
		gi_natural <- c(gi_natural, attempt_times)
		n_attempts_total <- n_attempts_total + nattempts

		# Simulate screening
		t_max_ind <- max(t_peak + window_after + screen_interval,
		                 max(attempt_times))
		phase <- runif(1, 0, screen_interval)
		test_times <- seq(phase, t_max_ind, by = screen_interval)

		# Find detection time
		t_detection <- Inf
		for (tt in test_times) {
			in_win <- (tt >= t_peak - window_before) &
				(tt <= t_peak + window_after)
			if (in_win) {
				if (runif(1) < sensitivity) {
					t_detection <- tt
					break
				}
			} else {
				if (runif(1) < fp_rate) {
					t_detection <- tt
					break
				}
			}
		}

		# Surviving attempts: those before detection
		surviving <- attempt_times[attempt_times < t_detection]
		gi_effective <- c(gi_effective, surviving)
		n_attempts_surviving <- n_attempts_surviving + length(surviving)
	}

	all_gi_natural[[ki]] <- gi_natural
	all_gi_effective[[ki]] <- gi_effective

	# Summary statistics
	mean_nat <- mean(gi_natural)
	mean_eff <- mean(gi_effective)
	R_eff <- R0 * (n_attempts_surviving / n_attempts_total)
	frac_averted <- 1 - n_attempts_surviving / n_attempts_total

	summary_results[[ki]] <- data.frame(
		kappa = kappa,
		mean_gi_natural = mean_nat,
		sd_gi_natural = sd(gi_natural),
		mean_gi_effective = mean_eff,
		sd_gi_effective = sd(gi_effective),
		gi_shift = mean_eff - mean_nat,
		gi_shift_pct = (mean_eff - mean_nat) / mean_nat * 100,
		R_eff = R_eff,
		frac_averted = frac_averted
	)

	cat(sprintf("  kappa=%4.1f: mean GI natural=%.2f, effective=%.2f (shift=%+.2f, %+.1f%%), R_eff=%.3f\n",
	            kappa, mean_nat, mean_eff,
	            mean_eff - mean_nat,
	            (mean_eff - mean_nat) / mean_nat * 100,
	            R_eff))
}

summary_df <- do.call(rbind, summary_results)
cat("\n")
print(summary_df)

# ==============================================================================
# Section 3: Euler-Lotka growth rate analysis
# ==============================================================================

cat("\n--- Euler-Lotka growth rate analysis ---\n\n")

#' Solve Euler-Lotka equation: 1 = R_eff * integral[exp(-r*tau) * g(tau) dtau]
#' where g(tau) is specified empirically via a vector of generation intervals.
#' Uses numerical root-finding. Allows negative r (declining epidemics).
solve_euler_lotka <- function(R_eff, gi_samples, r_range = c(-0.5, 5)) {
	if (length(gi_samples) < 10) return(NA_real_)

	euler_lotka_fn <- function(r_val) {
		R_eff * mean(exp(-r_val * gi_samples)) - 1
	}

	# Check that root exists in range
	f_lo <- euler_lotka_fn(r_range[1])
	f_hi <- euler_lotka_fn(r_range[2])
	if (sign(f_lo) == sign(f_hi)) {
		# Try wider range
		r_range <- c(-1.5, 10)
		f_lo <- euler_lotka_fn(r_range[1])
		f_hi <- euler_lotka_fn(r_range[2])
		if (sign(f_lo) == sign(f_hi)) return(NA_real_)
	}

	uniroot(euler_lotka_fn, interval = r_range, tol = 1e-8)$root
}

# --- Approach: Fixed R_eff thought experiment ---
# The simulated R_eff values are all < 1 (the intervention controls the
# epidemic). To isolate the GI distortion effect on growth rate, we ask a
# thought experiment question:
#
#   "For a given R_eff > 1, what growth rate would you predict using
#    the natural GI vs the intervention-distorted GI?"
#
# The gap between these two predictions, across kappa, shows how the
# GI distortion partially undermines the growth rate reduction for
# smooth profiles.

R_eff_values <- c(1.2, 1.5, 2.0, 2.5)

# Natural growth rate (no intervention)
r_natural <- solve_euler_lotka(R0, all_gi_natural[[1]])
cat(sprintf("  Natural growth rate (R0=%g): r = %.4f\n\n", R0, r_natural))

# For each R_eff value and kappa, compare growth rate with natural vs
# distorted GI
growth_list <- list()

for (R_eff_val in R_eff_values) {
	cat(sprintf("  R_eff = %.1f:\n", R_eff_val))
	for (ki in seq_along(kappas)) {
		# Growth rate using natural GI (what you'd predict from R_eff alone)
		r_natural_gi <- solve_euler_lotka(R_eff_val, all_gi_natural[[ki]])

		# Growth rate using distorted GI (what actually happens)
		r_distorted_gi <- solve_euler_lotka(R_eff_val, all_gi_effective[[ki]])

		offset <- r_distorted_gi - r_natural_gi

		growth_list[[length(growth_list) + 1]] <- data.frame(
			R_eff = R_eff_val,
			kappa = kappas[ki],
			r_natural_gi = r_natural_gi,
			r_distorted_gi = r_distorted_gi,
			offset = offset
		)

		cat(sprintf("    kappa=%4.1f: r(natural GI)=%.4f, r(distorted GI)=%.4f, offset=%+.4f\n",
		            kappas[ki], r_natural_gi, r_distorted_gi, offset))
	}
	cat("\n")
}

growth_df <- do.call(rbind, growth_list)
growth_df$R_eff_label <- factor(
	paste0("R[eff] == ", growth_df$R_eff),
	levels = paste0("R[eff] == ", R_eff_values)
)

print(growth_df)

# ==============================================================================
# Section 4: Figures
# ==============================================================================

cat("\n--- Generating figures ---\n")

# --- Figure 1: Generation interval density comparison -------------------------
# Show natural vs effective GI distributions for selected kappa values

selected_kappas <- c(0.5, 2, 6, 9.5)
tau_grid <- seq(0, 15, length.out = 500)

gi_density_df <- list()
for (kappa in selected_kappas) {
	ki <- which(kappas == kappa)

	# Natural: theoretical Gamma density
	gi_density_df[[length(gi_density_df) + 1]] <- data.frame(
		tau = tau_grid,
		density = dgamma(tau_grid, shape = alpha_total, rate = r),
		type = "Natural",
		kappa = kappa
	)

	# Effective: kernel density estimate from simulated data
	gi_eff <- all_gi_effective[[ki]]
	if (length(gi_eff) > 100) {
		kde <- density(gi_eff, from = 0, to = 15, n = 500, adjust = 1.2)
		gi_density_df[[length(gi_density_df) + 1]] <- data.frame(
			tau = kde$x,
			density = kde$y,
			type = "Under intervention",
			kappa = kappa
		)
	}
}

gi_density_df <- do.call(rbind, gi_density_df)
gi_density_df$kappa_label <- factor(
	sprintf("kappa == %s", gi_density_df$kappa),
	levels = sprintf("kappa == %s", selected_kappas)
)

fig_gi_densities <- ggplot(gi_density_df,
                           aes(x = tau, y = density, color = type,
                               linetype = type)) +
	geom_line(linewidth = 0.8) +
	facet_wrap(~kappa_label, nrow = 1, labeller = label_parsed,
	           scales = "free_y") +
	scale_color_manual(values = c("Natural" = "black",
	                              "Under intervention" = "red")) +
	scale_linetype_manual(values = c("Natural" = "dashed",
	                                 "Under intervention" = "solid")) +
	labs(x = expression(tau ~ "(days since infection)"),
	     y = "Density",
	     color = NULL, linetype = NULL) +
	theme_classic() +
	theme(legend.position = "bottom",
	      strip.background = element_blank())

ggsave("figures/fig_gi_distortion_densities.pdf", fig_gi_densities,
       width = 14, height = 4)
cat("  Saved figures/fig_gi_distortion_densities.pdf\n")

# --- Figure 2: Mean GI shift and R_eff vs kappa (two panels) -----------------

fig_shift <- ggplot(summary_df, aes(x = kappa, y = gi_shift)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
	labs(x = expression(kappa),
	     y = "Shift in mean generation\ninterval (days)",
	     tag = "A") +
	theme_classic()

fig_reff <- ggplot(summary_df, aes(x = kappa, y = R_eff)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	labs(x = expression(kappa),
	     y = expression(R[eff] ~ "under intervention"),
	     tag = "B") +
	theme_classic()

fig_shift_reff <- cowplot::plot_grid(fig_shift, fig_reff, nrow = 1)
ggsave("figures/fig_gi_distortion_shift.pdf", fig_shift_reff,
       width = 10, height = 4.5)
cat("  Saved figures/fig_gi_distortion_shift.pdf\n")

# --- Figure 3: Growth rate comparison (faceted by R_eff) ----------------------
# For a fixed R_eff, compare the growth rate implied by the natural GI vs
# the distorted GI. The gap shows the self-defeating effect of GI shortening.

growth_long <- growth_df %>%
	pivot_longer(cols = c(r_natural_gi, r_distorted_gi),
	             names_to = "scenario", values_to = "growth_rate") %>%
	mutate(scenario = recode(scenario,
	                         "r_natural_gi" = "Natural GI",
	                         "r_distorted_gi" = "Distorted GI"))

fig_growth <- ggplot(growth_long,
                     aes(x = kappa, y = growth_rate, color = scenario,
                         linetype = scenario)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 1.5) +
	facet_wrap(~R_eff_label, nrow = 1, labeller = label_parsed,
	           scales = "free_y") +
	scale_color_manual(values = c("Natural GI" = "blue",
	                              "Distorted GI" = "red")) +
	scale_linetype_manual(values = c("Natural GI" = "dashed",
	                                 "Distorted GI" = "solid")) +
	labs(x = expression(kappa),
	     y = "Epidemic growth rate (r)",
	     color = NULL, linetype = NULL) +
	theme_classic() +
	theme(legend.position = "bottom",
	      strip.background = element_blank())

ggsave("figures/fig_gi_distortion_growth.pdf", fig_growth,
       width = 14, height = 4.5)
cat("  Saved figures/fig_gi_distortion_growth.pdf\n")

# --- Figure 4: The offset (distorted - natural) vs kappa ---------------------
# Directly shows how much the GI shortening inflates the growth rate.
# Positive = GI distortion makes growth faster than R_eff alone would predict.

fig_offset <- ggplot(growth_df,
                     aes(x = kappa, y = offset,
                         color = factor(R_eff))) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
	scale_color_brewer(palette = "Set1") +
	labs(x = expression(kappa),
	     y = expression(r[distorted] - r[natural]),
	     color = expression(R[eff])) +
	theme_classic() +
	theme(legend.position = "bottom")

ggsave("figures/fig_gi_distortion_offset.pdf", fig_offset,
       width = 7, height = 5)
cat("  Saved figures/fig_gi_distortion_offset.pdf\n")

cat("\n=== Analysis complete ===\n")
