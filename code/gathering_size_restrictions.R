library(tidyverse)
library(cowplot)
source("code/utils.R")

# ==============================================================================
# Gathering size restrictions and punctuated infectiousness
#
# Key idea: spiky profiles sample the contact process at one time point,
# inheriting the full variance of gathering sizes (including the superspreading
# tail from large gatherings). Smooth profiles average over many gathering
# periods, yielding low-variance R_i regardless of gathering size distribution.
# Capping gathering sizes removes the large-gathering tail, disproportionately
# eliminating superspreading events for spiky profiles.
#
# Model: individual contact process c_i(t) is piecewise-constant (NegBin
# gathering sizes, Exp holding times). Interaction with biological profile
# a_i(tau) via Poisson thinning.
# ==============================================================================

set.seed(42)

# ==============================================================================
# Section 1: Parameters and helpers
# ==============================================================================

cat("=== Section 1: Parameters and helpers ===\n\n")

# Standard epi parameters
mu          <- 5
R0          <- 2
alpha_total <- 10
r           <- alpha_total / mu

# Punctuation sweep
kappas <- c(0.5, 1, 2, 4, 6, 8, 9.5)

# Contact process parameters
mu_gather   <- 30    # mean gathering size
k_nb        <- 1     # NegBin size parameter (geometric: heavy tail)
rate_change <- 2     # ~2 location changes per day

# Gathering size caps
caps <- c(5, 10, 20, Inf)

# Monte Carlo parameters
n_ind     <- 50000   # individual MC reps for R_i distribution
nsim_epi  <- 200     # epidemic simulations per scenario
n_pop     <- 1000    # population size for epidemics
epi_kappas <- c(0.5, 2, 6, 9.5)  # subset for epidemic sims

#' Generate a piecewise-constant contact process
#'
#' @param t_max   Maximum time horizon
#' @param mu_gather Mean gathering size (NegBin)
#' @param k_nb    NegBin size parameter
#' @param rate_change Rate of location changes (Exp holding times)
#' @param cap     Maximum gathering size (Inf = no cap)
#' @return List with t_breaks (change-point times) and sizes (gathering size per interval)
generate_contact_process <- function(t_max, mu_gather, k_nb, rate_change, cap) {
	t_breaks <- numeric(0)
	sizes    <- numeric(0)
	t_cur    <- 0

	while (t_cur < t_max) {
		# Draw gathering size
		if (is.finite(cap)) {
			# Rejection sampling from truncated NegBin
			repeat {
				g <- rnbinom(1, size = k_nb, mu = mu_gather)
				if (g <= cap) break
			}
		} else {
			g <- rnbinom(1, size = k_nb, mu = mu_gather)
		}

		t_breaks <- c(t_breaks, t_cur)
		sizes    <- c(sizes, g)

		# Holding time until next location change
		t_cur <- t_cur + rexp(1, rate_change)
	}

	list(t_breaks = t_breaks, sizes = sizes)
}

#' Evaluate contact process at a given time
#'
#' @param t       Time(s) to evaluate
#' @param process Output of generate_contact_process
#' @return Gathering size at each time t
evaluate_contact_process <- function(t, process) {
	idx <- findInterval(t, process$t_breaks)
	idx[idx < 1L] <- 1L
	process$sizes[idx]
}

#' Factory: shifted Gamma profile with stochastic gathering size contact process
#'
#' Uses Poisson thinning: oversample by c_max / mu_gather, then accept each
#' candidate with probability c_i(tau) / c_max.
#'
#' @param mu          Mean generation time
#' @param R0          Basic reproduction number
#' @param alpha_total Population kernel shape
#' @param kappa       Profile shape (small = spiky)
#' @param mu_gather   Mean gathering size
#' @param k_nb        NegBin size parameter
#' @param rate_change Rate of location changes
#' @param cap         Gathering size cap (Inf = no cap)
#' @return A closure function(t_inf) -> sorted infection attempt times
make_profile_gamma_gathering <- function(mu, R0, alpha_total, kappa,
                                         mu_gather, k_nb, rate_change, cap) {
	stopifnot(kappa > 0, kappa < alpha_total)
	r_rate <- alpha_total / mu
	alpha_shift <- alpha_total - kappa

	# c_max for thinning
	if (is.finite(cap)) {
		c_max <- cap
	} else {
		c_max <- qnbinom(0.999, size = k_nb, mu = mu_gather)
	}

	function(t_inf) {
		# Oversample
		n_proposal <- rpois(1, R0 * c_max / mu_gather)
		if (n_proposal == 0L) return(numeric(0))

		# Onset shift (shared)
		s_i <- rgamma(1, shape = alpha_shift, rate = r_rate)

		# Independent jitter for each candidate
		eps <- rgamma(n_proposal, shape = kappa, rate = r_rate)

		# Absolute candidate times
		tau_candidates <- s_i + eps

		# Generate contact process over the span of candidate times
		t_max_proc <- max(tau_candidates) + 1
		proc <- generate_contact_process(t_max_proc, mu_gather, k_nb, rate_change, cap)

		# Evaluate contact process at each candidate time
		c_vals <- evaluate_contact_process(tau_candidates, proc)

		# Thinning: accept with probability c_i(tau) / c_max
		accept_prob <- c_vals / c_max
		keep <- runif(n_proposal) < accept_prob

		if (!any(keep)) return(numeric(0))
		sort.int(t_inf + tau_candidates[keep])
	}
}

# ==============================================================================
# Section 2: Individual-level R_i distribution
# ==============================================================================

cat("=== Section 2: Individual-level R_i distribution ===\n\n")

# Cap labels for display
cap_labels <- ifelse(is.finite(caps), paste0("Cap = ", caps), "No cap")
names(cap_labels) <- caps

ri_results <- expand_grid(kappa = kappas, cap = caps) %>%
	mutate(cap_label = cap_labels[as.character(cap)],
	       mean_ri   = NA_real_,
	       var_ri    = NA_real_,
	       cv_ri     = NA_real_,
	       vmr_ri    = NA_real_,
	       p_ge5     = NA_real_)

# Store full distributions for histogram plotting
ri_dist_list <- list()

for (idx in seq_len(nrow(ri_results))) {
	kp  <- ri_results$kappa[idx]
	cap <- ri_results$cap[idx]

	alpha_shift <- alpha_total - kp

	# c_max for thinning
	if (is.finite(cap)) {
		c_max <- cap
	} else {
		c_max <- qnbinom(0.999, size = k_nb, mu = mu_gather)
	}

	ri_vec <- numeric(n_ind)
	for (i in seq_len(n_ind)) {
		# Draw proposal count
		n_proposal <- rpois(1, R0 * c_max / mu_gather)
		if (n_proposal == 0L) { ri_vec[i] <- 0L; next }

		# Onset shift
		s_i <- rgamma(1, shape = alpha_shift, rate = r)

		# Jitter per candidate
		eps <- rgamma(n_proposal, shape = kp, rate = r)
		tau_candidates <- s_i + eps

		# Contact process
		t_max_proc <- max(tau_candidates) + 1
		proc <- generate_contact_process(t_max_proc, mu_gather, k_nb, rate_change, cap)
		c_vals <- evaluate_contact_process(tau_candidates, proc)

		# Thinning
		accept_prob <- c_vals / c_max
		ri_vec[i] <- sum(runif(n_proposal) < accept_prob)
	}

	ri_results$mean_ri[idx] <- mean(ri_vec)
	ri_results$var_ri[idx]  <- var(ri_vec)
	ri_results$cv_ri[idx]   <- sd(ri_vec) / mean(ri_vec)
	ri_results$vmr_ri[idx]  <- var(ri_vec) / mean(ri_vec)
	ri_results$p_ge5[idx]   <- mean(ri_vec >= 5)

	ri_dist_list[[paste(kp, cap, sep = "_")]] <- ri_vec

	if (idx %% 4 == 0) cat(sprintf("  %d / %d scenarios done\n", idx, nrow(ri_results)))
}

cat("\nDone computing R_i distributions.\n\n")

# Print summary
cat("=== R_i summary table ===\n\n")
ri_results %>%
	select(kappa, cap_label, mean_ri, cv_ri, vmr_ri, p_ge5) %>%
	print(n = Inf)

# Verification: mean R_i should be roughly constant across kappa for each cap
cat("\n=== Verification: Mean R_i by cap (should be ~constant across kappa) ===\n\n")
ri_results %>%
	group_by(cap_label) %>%
	summarise(
		min_mean = min(mean_ri),
		max_mean = max(mean_ri),
		range    = max(mean_ri) - min(mean_ri),
		.groups = "drop"
	) %>%
	print()

# ==============================================================================
# Section 3: Epidemic simulations
# ==============================================================================

cat("\n=== Section 3: Epidemic simulations ===\n\n")

epi_df <- tibble()

for (kp in epi_kappas) {
	for (cap in caps) {
		cap_lab <- cap_labels[as.character(cap)]
		cat(sprintf("  Running: kappa=%.1f, %s\n", kp, cap_lab))

		gen <- make_profile_gamma_gathering(
			mu = mu, R0 = R0, alpha_total = alpha_total, kappa = kp,
			mu_gather = mu_gather, k_nb = k_nb, rate_change = rate_change, cap = cap
		)

		for (sim in seq_len(nsim_epi)) {
			tinf <- sim_stochastic_fast(n = n_pop, gen_inf_attempts = gen)
			n_infected <- sum(tinf < Inf)
			epi_df <- bind_rows(epi_df, tibble(
				sim       = sim,
				kappa     = kp,
				cap       = cap,
				cap_label = cap_lab,
				final_size = n_infected
			))
		}
	}
}

epi_df <- epi_df %>%
	mutate(established = final_size >= 0.1 * n_pop)

cat("\nDone with epidemic simulations.\n\n")

# Summary
cat("=== Epidemic summary ===\n\n")
epi_summary <- epi_df %>%
	group_by(kappa, cap_label) %>%
	summarise(
		mean_fs       = mean(final_size),
		p_established = mean(established),
		.groups = "drop"
	)
epi_summary %>% print(n = Inf)

# ==============================================================================
# Section 4: Summary decomposition
# ==============================================================================

cat("\n=== Section 4: Differential effect decomposition ===\n\n")

# Join to no-cap baseline
baseline <- ri_results %>%
	filter(cap == Inf) %>%
	select(kappa, mean_ri_base = mean_ri, cv_ri_base = cv_ri)

diff_df <- ri_results %>%
	filter(is.finite(cap)) %>%
	left_join(baseline, by = "kappa") %>%
	mutate(
		frac_red_mean = 1 - mean_ri / mean_ri_base,
		frac_red_cv   = 1 - cv_ri / cv_ri_base
	)

cat("=== Fractional reduction in E[R_i] (should be ~constant across kappa) ===\n\n")
diff_df %>%
	select(kappa, cap_label, frac_red_mean) %>%
	pivot_wider(names_from = cap_label, values_from = frac_red_mean) %>%
	print(n = Inf)

cat("\n=== Fractional reduction in CV(R_i) (should slope with kappa) ===\n\n")
diff_df %>%
	select(kappa, cap_label, frac_red_cv) %>%
	pivot_wider(names_from = cap_label, values_from = frac_red_cv) %>%
	print(n = Inf)

# Epidemic differential
epi_baseline <- epi_summary %>%
	filter(cap_label == "No cap") %>%
	select(kappa, mean_fs_base = mean_fs, p_est_base = p_established)

epi_diff <- epi_summary %>%
	filter(cap_label != "No cap") %>%
	left_join(epi_baseline, by = "kappa") %>%
	mutate(
		frac_red_fs  = 1 - mean_fs / mean_fs_base,
		frac_red_est = 1 - p_established / p_est_base
	)

cat("\n=== Epidemic: fractional reduction in mean final size ===\n\n")
epi_diff %>%
	select(kappa, cap_label, frac_red_fs) %>%
	pivot_wider(names_from = cap_label, values_from = frac_red_fs) %>%
	print(n = Inf)

# ==============================================================================
# Section 5: Figures
# ==============================================================================

cat("\n=== Section 5: Generating figures ===\n\n")

# Color palette
cap_colors <- c("Cap = 5" = "red", "Cap = 10" = "orange",
                "Cap = 20" = "steelblue", "No cap" = "black")

# Order factor for cap labels
cap_level_order <- c("Cap = 5", "Cap = 10", "Cap = 20", "No cap")

ri_results <- ri_results %>%
	mutate(cap_label = factor(cap_label, levels = cap_level_order))
diff_df <- diff_df %>%
	mutate(cap_label = factor(cap_label, levels = cap_level_order))

# --- Figure 1: Faceted histograms of R_i ---

hist_kappas <- c(0.5, 2, 6, 9.5)
hist_df <- tibble()
for (kp in hist_kappas) {
	for (cap in caps) {
		key <- paste(kp, cap, sep = "_")
		cap_lab <- cap_labels[as.character(cap)]
		ri_vec <- ri_dist_list[[key]]
		hist_df <- bind_rows(hist_df, tibble(
			ri        = ri_vec,
			kappa     = kp,
			cap       = cap,
			cap_label = cap_lab
		))
	}
}

hist_df <- hist_df %>%
	mutate(kappa_label = paste0("kappa == ", kappa),
	       cap_label   = factor(cap_label, levels = cap_level_order))

fig_ri_dist <- hist_df %>%
	ggplot(aes(x = ri)) +
	geom_histogram(aes(y = after_stat(density)), binwidth = 1,
	               fill = "steelblue", col = "white", alpha = 0.7) +
	facet_grid(kappa_label ~ cap_label,
	           labeller = labeller(kappa_label = label_parsed, cap_label = label_value)) +
	theme_classic() +
	coord_cartesian(xlim = c(0, 15)) +
	labs(x = expression(R[i]),
	     y = "Density",
	     title = "Individual reproduction number distributions",
	     subtitle = "Rows = kappa (spiky to smooth), Columns = gathering size cap")

ggsave("figures/fig_gathering_ri_dist.pdf", fig_ri_dist, width = 14, height = 12)
cat("Saved figures/fig_gathering_ri_dist.pdf\n")

# --- Figure 2: R_i summary (mean and CV) ---

p_mean <- ri_results %>%
	ggplot(aes(x = kappa, y = mean_ri, col = cap_label, group = cap_label)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	scale_color_manual(values = cap_colors, name = "Restriction") +
	theme_classic() +
	theme(legend.position = "none") +
	labs(x = expression(kappa), y = expression(E(R[i])),
	     title = expression("(A) Mean " * R[i] * " vs " * kappa))

p_cv <- ri_results %>%
	ggplot(aes(x = kappa, y = cv_ri, col = cap_label, group = cap_label)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	scale_color_manual(values = cap_colors, name = "Restriction") +
	theme_classic() +
	theme(legend.position = "none") +
	labs(x = expression(kappa), y = expression(CV(R[i])),
	     title = expression("(B) CV of " * R[i] * " vs " * kappa))

# Extract shared legend
legend_plot <- ri_results %>%
	ggplot(aes(x = kappa, y = mean_ri, col = cap_label)) +
	geom_point() +
	scale_color_manual(values = cap_colors, name = "Restriction") +
	theme_classic() +
	theme(legend.position = "bottom")
shared_legend <- get_legend(legend_plot)

fig_ri_summary <- plot_grid(
	plot_grid(p_mean, p_cv, ncol = 2),
	shared_legend,
	ncol = 1, rel_heights = c(1, 0.1)
)

ggsave("figures/fig_gathering_ri_summary.pdf", fig_ri_summary, width = 12, height = 5)
cat("Saved figures/fig_gathering_ri_summary.pdf\n")

# --- Figure 3: Epidemic outcomes ---

epi_summary <- epi_summary %>%
	mutate(cap_label = factor(cap_label, levels = cap_level_order))

p_fs <- epi_summary %>%
	ggplot(aes(x = kappa, y = mean_fs, col = cap_label, group = cap_label)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	scale_color_manual(values = cap_colors, name = "Restriction") +
	theme_classic() +
	theme(legend.position = "none") +
	labs(x = expression(kappa), y = "Mean final size",
	     title = "(A) Mean final size vs kappa")

p_est <- epi_summary %>%
	ggplot(aes(x = kappa, y = p_established, col = cap_label, group = cap_label)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	scale_color_manual(values = cap_colors, name = "Restriction") +
	theme_classic() +
	theme(legend.position = "none") +
	labs(x = expression(kappa), y = "P(establishment)",
	     title = "(B) P(establishment) vs kappa")

fig_epidemic <- plot_grid(
	plot_grid(p_fs, p_est, ncol = 2),
	shared_legend,
	ncol = 1, rel_heights = c(1, 0.1)
)

ggsave("figures/fig_gathering_epidemic.pdf", fig_epidemic, width = 12, height = 5)
cat("Saved figures/fig_gathering_epidemic.pdf\n")

# --- Figure 4: Differential effect ---

diff_df <- diff_df %>%
	mutate(cap_label = factor(cap_label, levels = cap_level_order))

p_red_mean <- diff_df %>%
	ggplot(aes(x = kappa, y = frac_red_mean, col = cap_label, group = cap_label)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	geom_hline(yintercept = 0, linetype = "dashed", col = "grey50") +
	scale_color_manual(values = cap_colors, name = "Restriction") +
	theme_classic() +
	theme(legend.position = "none") +
	labs(x = expression(kappa),
	     y = expression("Fractional reduction in E[" * R[i] * "]"),
	     title = expression("(A) Reduction in mean " * R[i] * " (should be flat)"))

p_red_cv <- diff_df %>%
	ggplot(aes(x = kappa, y = frac_red_cv, col = cap_label, group = cap_label)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	geom_hline(yintercept = 0, linetype = "dashed", col = "grey50") +
	scale_color_manual(values = cap_colors, name = "Restriction") +
	theme_classic() +
	theme(legend.position = "none") +
	labs(x = expression(kappa),
	     y = expression("Fractional reduction in CV(" * R[i] * ")"),
	     title = expression("(B) Reduction in CV (differential effect)"))

# Shared legend for differential plot (only finite caps)
legend_diff <- diff_df %>%
	ggplot(aes(x = kappa, y = frac_red_mean, col = cap_label)) +
	geom_point() +
	scale_color_manual(values = cap_colors, name = "Restriction") +
	theme_classic() +
	theme(legend.position = "bottom")
shared_legend_diff <- get_legend(legend_diff)

fig_differential <- plot_grid(
	plot_grid(p_red_mean, p_red_cv, ncol = 2),
	shared_legend_diff,
	ncol = 1, rel_heights = c(1, 0.1)
)

ggsave("figures/fig_gathering_differential.pdf", fig_differential, width = 12, height = 5)
cat("Saved figures/fig_gathering_differential.pdf\n")

cat("\nDone.\n")
