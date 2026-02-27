library(tidyverse)
source("code/utils.R")

# ==============================================================================
# Why does the same test-and-isolate intervention yield different R_eff
# reductions depending on kappa?
#
# The mechanism has two components:
#
# 1. ALL-OR-NOTHING VS PARTIAL: For spike profiles, detection before the spike
#    averts 100% of transmission; detection after averts 0%. For smooth profiles,
#    detection averts only the remaining fraction — some transmission has already
#    occurred. The conditional E[frac averted | detected] is lower for smooth.
#
# 2. TAIL LEAKAGE: Smooth profiles have transmission outside the positivity
#    window (before it opens or after it closes), which is uncatchable.
#    Spike profiles concentrate everything near the peak.
#
# This script decomposes the R_eff reduction into these components.
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

sensitivity <- 0.85
fp_rate <- 0.02
window_before <- 3
window_after <- 3
screen_interval <- 3

n_ind <- 100000

get_peak_offset <- function(kappa, r) max(0, (kappa - 1) / r)

cat("=== R_eff reduction decomposition ===\n")
cat(sprintf("  mu=%g, R0=%g, alpha_total=%g, r=%g\n", mu, R0, alpha_total, r))
cat(sprintf("  screen_interval=%d, sensitivity=%.2f, window=[-%g, +%g]\n",
            screen_interval, sensitivity, window_before, window_after))
cat(sprintf("  n_ind=%d\n\n", n_ind))

# ==============================================================================
# Section 2: Per-individual simulation with detailed tracking
# ==============================================================================

decomp_results <- list()
frac_averted_distributions <- list()  # store per-individual frac_averted

for (ki in seq_along(kappas)) {
	kappa <- kappas[ki]
	alpha_shift <- alpha_total - kappa
	peak_offset <- get_peak_offset(kappa, r)

	# Per-individual tracking vectors
	n_attempts_vec <- integer(n_ind)
	n_pre_window <- integer(n_ind)     # attempts before positivity window
	n_in_window <- integer(n_ind)      # attempts within positivity window
	n_post_window <- integer(n_ind)    # attempts after positivity window
	detected <- logical(n_ind)
	detected_before_peak <- logical(n_ind)
	frac_averted <- numeric(n_ind)
	frac_averted_given_det <- numeric(n_ind)

	for (i in 1:n_ind) {
		s_i <- rgamma(1, shape = alpha_shift, rate = r)
		nattempts <- rpois(1, R0)
		n_attempts_vec[i] <- nattempts

		if (nattempts == 0L) {
			frac_averted[i] <- NA  # no attempts to avert
			next
		}

		# Attempt times (relative to infection at t=0)
		jitters <- rgamma(nattempts, shape = kappa, rate = r)
		attempt_times <- s_i + sort(jitters)
		t_peak <- s_i + peak_offset

		# Window boundaries
		win_start <- t_peak - window_before
		win_end <- t_peak + window_after

		# Classify attempts by timing relative to window
		pre <- sum(attempt_times < win_start)
		post <- sum(attempt_times > win_end)
		within <- nattempts - pre - post
		n_pre_window[i] <- pre
		n_in_window[i] <- within
		n_post_window[i] <- post

		# Simulate screening
		t_max_ind <- max(win_end + screen_interval, max(attempt_times))
		phase <- runif(1, 0, screen_interval)
		test_times <- seq(phase, t_max_ind, by = screen_interval)

		# Find detection time
		t_detection <- Inf
		for (tt in test_times) {
			in_win <- (tt >= win_start) & (tt <= win_end)
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

		if (is.finite(t_detection)) {
			detected[i] <- TRUE
			detected_before_peak[i] <- (t_detection < t_peak)
			frac_averted[i] <- mean(attempt_times > t_detection)
		} else {
			detected[i] <- FALSE
			frac_averted[i] <- 0
		}
	}

	# Store distribution of frac_averted for detected individuals with attempts
	has_attempts <- n_attempts_vec > 0
	det_with_attempts <- detected & has_attempts
	frac_averted_distributions[[ki]] <- data.frame(
		kappa = kappa,
		frac_averted = frac_averted[det_with_attempts]
	)

	# Summary statistics
	has_att <- has_attempts
	total_attempts <- sum(n_attempts_vec[has_att])
	total_pre <- sum(n_pre_window[has_att])
	total_in <- sum(n_in_window[has_att])
	total_post <- sum(n_post_window[has_att])

	decomp_results[[ki]] <- data.frame(
		kappa = kappa,
		# Detection statistics
		p_detect = mean(detected[has_att]),
		p_detect_before_peak = mean(detected_before_peak[has_att]),
		p_detect_after_peak = mean(detected[has_att] & !detected_before_peak[has_att]),
		# Conditional fractions averted
		frac_averted_unconditional = mean(frac_averted[has_att]),
		frac_averted_given_detected = mean(frac_averted[det_with_attempts]),
		frac_averted_given_det_before = ifelse(
			sum(detected_before_peak[has_att]) > 0,
			mean(frac_averted[has_att & detected_before_peak]),
			NA),
		frac_averted_given_det_after = ifelse(
			sum(detected[has_att] & !detected_before_peak[has_att]) > 0,
			mean(frac_averted[has_att & detected & !detected_before_peak]),
			NA),
		# Attempt location decomposition
		frac_attempts_pre_window = total_pre / total_attempts,
		frac_attempts_in_window = total_in / total_attempts,
		frac_attempts_post_window = total_post / total_attempts,
		# R_eff
		R_eff = R0 * (1 - mean(frac_averted[has_att]))
	)

	cat(sprintf("  kappa=%4.1f: P(det)=%.3f, P(det before peak)=%.3f, E[frac|det]=%.3f, R_eff=%.3f\n",
	            kappa,
	            decomp_results[[ki]]$p_detect,
	            decomp_results[[ki]]$p_detect_before_peak,
	            decomp_results[[ki]]$frac_averted_given_detected,
	            decomp_results[[ki]]$R_eff))
}

decomp_df <- do.call(rbind, decomp_results)
cat("\n")

# Print full decomposition table
cat("=== Full decomposition ===\n\n")
cat("Detection:\n")
print(decomp_df[, c("kappa", "p_detect", "p_detect_before_peak",
                     "p_detect_after_peak")])
cat("\nFraction averted:\n")
print(decomp_df[, c("kappa", "frac_averted_unconditional",
                     "frac_averted_given_detected",
                     "frac_averted_given_det_before",
                     "frac_averted_given_det_after")])
cat("\nAttempt location:\n")
print(decomp_df[, c("kappa", "frac_attempts_pre_window",
                     "frac_attempts_in_window",
                     "frac_attempts_post_window")])
cat("\nR_eff:\n")
print(decomp_df[, c("kappa", "R_eff")])

# ==============================================================================
# Section 3: Figures
# ==============================================================================

cat("\n--- Generating figures ---\n")

# --- Figure 1: Distribution of frac_averted given detection ------------------
# This is the key figure: spike profiles show bimodal (0 or 1) distributions,
# smooth profiles show unimodal distributions centered below 1.

frac_dist_df <- do.call(rbind, frac_averted_distributions)
selected_kappas <- c(0.5, 2, 6, 9.5)
frac_dist_selected <- frac_dist_df %>%
	filter(kappa %in% selected_kappas) %>%
	mutate(kappa_label = factor(
		sprintf("kappa == %s", kappa),
		levels = sprintf("kappa == %s", selected_kappas)
	))

fig_dist <- ggplot(frac_dist_selected,
                   aes(x = frac_averted)) +
	geom_histogram(aes(y = after_stat(density)),
	               bins = 50, fill = "steelblue", color = "white",
	               linewidth = 0.2) +
	facet_wrap(~kappa_label, nrow = 1, labeller = label_parsed,
	           scales = "free_y") +
	labs(x = "Fraction of transmission averted (given detection)",
	     y = "Density") +
	theme_classic() +
	theme(strip.background = element_blank())

ggsave("figures/fig_reff_frac_averted_dist.pdf", fig_dist,
       width = 14, height = 4)
cat("  Saved figures/fig_reff_frac_averted_dist.pdf\n")

# --- Figure 2: Decomposition of the R_eff reduction mechanism ----------------
# Three-panel: (A) P(detected) & P(detected before peak),
#              (B) E[frac averted | detected],
#              (C) R_eff

fig_2a <- decomp_df %>%
	select(kappa, p_detect, p_detect_before_peak) %>%
	pivot_longer(cols = c(p_detect, p_detect_before_peak),
	             names_to = "type", values_to = "probability") %>%
	mutate(type = recode(type,
	                     "p_detect" = "P(detected)",
	                     "p_detect_before_peak" = "P(detected before peak)")) %>%
	ggplot(aes(x = kappa, y = probability, color = type, linetype = type)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	scale_color_manual(values = c("P(detected)" = "black",
	                              "P(detected before peak)" = "red")) +
	scale_linetype_manual(values = c("P(detected)" = "solid",
	                                 "P(detected before peak)" = "dashed")) +
	labs(x = expression(kappa), y = "Probability",
	     color = NULL, linetype = NULL, tag = "A") +
	theme_classic() +
	theme(legend.position = c(0.65, 0.3),
	      legend.background = element_blank())

fig_2b <- ggplot(decomp_df,
                 aes(x = kappa, y = frac_averted_given_detected)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	labs(x = expression(kappa),
	     y = "E[frac averted | detected]",
	     tag = "B") +
	theme_classic()

fig_2c <- ggplot(decomp_df, aes(x = kappa, y = R_eff)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	geom_hline(yintercept = R0, linetype = "dotted", color = "grey50") +
	annotate("text", x = 0.5, y = R0 + 0.03,
	         label = paste0("R0 = ", R0), hjust = 0, size = 3, color = "grey50") +
	labs(x = expression(kappa),
	     y = expression(R[eff]),
	     tag = "C") +
	theme_classic()

fig_decomp <- cowplot::plot_grid(fig_2a, fig_2b, fig_2c, nrow = 1,
                                  rel_widths = c(1.2, 1, 1))
ggsave("figures/fig_reff_decomposition.pdf", fig_decomp,
       width = 14, height = 4.5)
cat("  Saved figures/fig_reff_decomposition.pdf\n")

# --- Figure 3: Attempt location decomposition --------------------------------
# Stacked bar showing fraction of attempts pre-window / in-window / post-window

attempt_loc_df <- decomp_df %>%
	select(kappa, frac_attempts_pre_window, frac_attempts_in_window,
	       frac_attempts_post_window) %>%
	pivot_longer(cols = -kappa, names_to = "location", values_to = "fraction") %>%
	mutate(location = recode(location,
	                         "frac_attempts_pre_window" = "Before window\n(uncatchable)",
	                         "frac_attempts_in_window" = "Within window\n(catchable)",
	                         "frac_attempts_post_window" = "After window\n(uncatchable)")) %>%
	mutate(location = factor(location,
	                         levels = c("After window\n(uncatchable)",
	                                    "Within window\n(catchable)",
	                                    "Before window\n(uncatchable)")))

fig_location <- ggplot(attempt_loc_df,
                       aes(x = factor(kappa), y = fraction, fill = location)) +
	geom_col(position = "stack", width = 0.7) +
	scale_fill_manual(values = c("Before window\n(uncatchable)" = "firebrick",
	                             "Within window\n(catchable)" = "steelblue",
	                             "After window\n(uncatchable)" = "darkorange")) +
	labs(x = expression(kappa),
	     y = "Fraction of transmission attempts",
	     fill = "Location relative\nto positivity window") +
	theme_classic() +
	theme(legend.position = "right")

ggsave("figures/fig_reff_attempt_location.pdf", fig_location,
       width = 8, height = 5)
cat("  Saved figures/fig_reff_attempt_location.pdf\n")

# --- Figure 4: Conditional frac averted by detection timing -------------------
# E[frac averted | detected before peak] vs E[frac averted | detected after peak]

cond_df <- decomp_df %>%
	select(kappa, frac_averted_given_det_before, frac_averted_given_det_after) %>%
	pivot_longer(cols = -kappa, names_to = "timing", values_to = "frac") %>%
	mutate(timing = recode(timing,
	                       "frac_averted_given_det_before" = "Detected before peak",
	                       "frac_averted_given_det_after" = "Detected after peak"))

fig_conditional <- ggplot(cond_df,
                          aes(x = kappa, y = frac, color = timing,
                              linetype = timing)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	scale_color_manual(values = c("Detected before peak" = "blue",
	                              "Detected after peak" = "red")) +
	scale_linetype_manual(values = c("Detected before peak" = "solid",
	                                 "Detected after peak" = "dashed")) +
	labs(x = expression(kappa),
	     y = "E[frac averted | detected, timing]",
	     color = NULL, linetype = NULL) +
	theme_classic() +
	theme(legend.position = "bottom")

ggsave("figures/fig_reff_conditional.pdf", fig_conditional,
       width = 8, height = 5)
cat("  Saved figures/fig_reff_conditional.pdf\n")

cat("\n=== Analysis complete ===\n")
