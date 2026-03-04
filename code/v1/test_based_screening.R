library(tidyverse)
source("code/utils.R")

# ==============================================================================
# Impact of punctuated infectiousness on test-based screening and isolation
#
# Corresponds to writeup section:
#   "The impact of punctuated infectiousness on test-based screening and
#    isolation countermeasures"
#
# Analyzes how the punctuatedness (kappa) of individual infectiousness
# profiles affects the effectiveness of:
#   1. Regular periodic screening
#   2. Pre-event screening
#   3. Test-to-exit-isolation
#   4. Epidemic-level screening interventions
# ==============================================================================

set.seed(2024)

# ==============================================================================
# Section 1: Parameters and helpers
# ==============================================================================

# Standard epidemiological parameters (matching codebase conventions)
mu <- 5
R0 <- 2
alpha_total <- 10
r <- alpha_total / mu

# Kappa sweep (small = punctuated, large = smooth)
kappas <- c(0.5, 1, 2, 4, 6, 8, 9.5)

# Test parameters
sensitivity <- 0.85
fp_rate <- 0.02
window_before <- 3    # days before peak where test can detect
window_after <- 4     # days after peak where test can detect

# Screening intervals (days)
screen_intervals <- c(1, 3, 7)

# Monte Carlo sizes
n_ind <- 50000     # individual-level MC reps
nsim_epi <- 200    # epidemic simulations per scenario

cat("=== Test-based screening analysis ===\n")
cat(sprintf("  mu=%g, R0=%g, alpha_total=%g, r=%g\n", mu, R0, alpha_total, r))
cat(sprintf("  sensitivity=%.2f, fp_rate=%.2f, window=[-%g, +%g]\n",
            sensitivity, fp_rate, window_before, window_after))
cat(sprintf("  n_ind=%d, nsim_epi=%d\n\n", n_ind, nsim_epi))

# --- Helper functions ---------------------------------------------------------

#' Peak offset of the individual Gamma(kappa, r) profile relative to onset s_i
get_peak_offset <- function(kappa, r) {
	max(0, (kappa - 1) / r)
}

#' Run a single diagnostic test
#'
#' Returns TRUE (positive) if t_test is within the positivity window around
#' t_peak and the test fires (with prob = sensitivity), OR if it's outside
#' the window and we get a false positive (with prob = fp_rate).
run_test <- function(t_test, t_peak, sensitivity, fp_rate,
                     window_before, window_after) {
	in_window <- (t_test >= t_peak - window_before) &
		(t_test <= t_peak + window_after)
	if (in_window) {
		return(runif(1) < sensitivity)
	} else {
		return(runif(1) < fp_rate)
	}
}

#' Vectorized test: test one time point against multiple peak times
#' Returns logical vector of length(t_peak_vector)
run_test_vec <- function(t_test_scalar, t_peak_vector, sensitivity, fp_rate,
                         window_before, window_after) {
	in_window <- (t_test_scalar >= t_peak_vector - window_before) &
		(t_test_scalar <= t_peak_vector + window_after)
	result <- logical(length(t_peak_vector))
	result[in_window] <- runif(sum(in_window)) < sensitivity
	result[!in_window] <- runif(sum(!in_window)) < fp_rate
	result
}

#' Generate one individual's infection profile details
#'
#' Returns a list with onset shift s_i, peak time t_peak (absolute),
#' number of infection attempts, and their times.
generate_individual <- function(kappa, r, alpha_total, R0, t_inf = 0) {
	alpha_shift <- alpha_total - kappa
	s_i <- rgamma(1, shape = alpha_shift, rate = r)
	nattempts <- rpois(1, R0)
	if (nattempts == 0L) {
		attempt_times <- numeric(0)
	} else {
		attempt_times <- t_inf + s_i + sort(rgamma(nattempts, shape = kappa, rate = r))
	}
	peak_offset <- get_peak_offset(kappa, r)
	t_peak <- t_inf + s_i + peak_offset
	list(
		s_i = s_i,
		t_peak = t_peak,
		nattempts = nattempts,
		attempt_times = attempt_times
	)
}

# ==============================================================================
# Section 2: Individual-level — Regular periodic screening
# ==============================================================================

cat("--- Section 2: Regular periodic screening ---\n")

regular_results <- list()
idx <- 0

for (kappa in kappas) {
	for (interval in screen_intervals) {

		detected <- logical(n_ind)
		frac_averted <- numeric(n_ind)
		is_false_pos <- logical(n_ind)

		for (i in 1:n_ind) {
			ind <- generate_individual(kappa, r, alpha_total, R0, t_inf = 0)

			# Maximum time of interest: end of positivity window or last attempt
			t_max_ind <- max(c(ind$t_peak + window_after + interval,
			                   ind$attempt_times, 0))

			# Random screening phase
			phase <- runif(1, 0, interval)
			test_times <- seq(phase, t_max_ind, by = interval)

			# Test sequentially until detected or done
			t_detection <- Inf
			detection_in_window <- FALSE
			for (tt in test_times) {
				in_win <- (tt >= ind$t_peak - window_before) &
					(tt <= ind$t_peak + window_after)
				if (in_win) {
					if (runif(1) < sensitivity) {
						t_detection <- tt
						detection_in_window <- TRUE
						break
					}
				} else {
					if (runif(1) < fp_rate) {
						t_detection <- tt
						detection_in_window <- FALSE
						break
					}
				}
			}

			if (is.finite(t_detection)) {
				detected[i] <- TRUE
				is_false_pos[i] <- !detection_in_window
				if (ind$nattempts > 0) {
					frac_averted[i] <- mean(ind$attempt_times > t_detection)
				} else {
					frac_averted[i] <- 0
				}
			} else {
				detected[i] <- FALSE
				frac_averted[i] <- 0
			}
		}

		idx <- idx + 1
		regular_results[[idx]] <- data.frame(
			kappa = kappa,
			interval = interval,
			p_detect = mean(detected),
			frac_averted_all = mean(frac_averted),
			frac_averted_given_detected = ifelse(sum(detected) > 0,
			                                     mean(frac_averted[detected]), NA),
			p_false_pos = mean(is_false_pos[detected])
		)

		cat(sprintf("  kappa=%.1f, interval=%d: p_detect=%.3f, frac_averted=%.3f\n",
		            kappa, interval,
		            regular_results[[idx]]$p_detect,
		            regular_results[[idx]]$frac_averted_all))
	}
}

regular_results <- do.call(rbind, regular_results)
regular_results$interval_label <- factor(
	paste0("Every ", regular_results$interval, " day(s)"),
	levels = paste0("Every ", screen_intervals, " day(s)")
)

cat("\n")
print(regular_results)
cat("\n")

# ==============================================================================
# Section 3: Individual-level — Pre-event screening
# ==============================================================================

cat("--- Section 3: Pre-event screening ---\n")

preevent_results <- list()

for (ki in seq_along(kappas)) {
	kappa <- kappas[ki]

	T_max <- qgamma(0.99, shape = alpha_total, rate = r)

	detected <- logical(n_ind)
	frac_averted <- numeric(n_ind)

	for (i in 1:n_ind) {
		ind <- generate_individual(kappa, r, alpha_total, R0, t_inf = 0)

		# Single test at random time
		t_test <- runif(1, 0, T_max)

		in_win <- (t_test >= ind$t_peak - window_before) &
			(t_test <= ind$t_peak + window_after)
		if (in_win) {
			pos <- runif(1) < sensitivity
		} else {
			pos <- runif(1) < fp_rate
		}

		if (pos) {
			detected[i] <- TRUE
			if (ind$nattempts > 0) {
				frac_averted[i] <- mean(ind$attempt_times > t_test)
			} else {
				frac_averted[i] <- 0
			}
		} else {
			detected[i] <- FALSE
			frac_averted[i] <- 0
		}
	}

	preevent_results[[ki]] <- data.frame(
		kappa = kappa,
		p_detect = mean(detected),
		frac_averted_given_detected = ifelse(sum(detected) > 0,
		                                     mean(frac_averted[detected]), NA),
		frac_averted_all = mean(frac_averted)
	)

	cat(sprintf("  kappa=%.1f: p_detect=%.3f, frac_averted_all=%.3f\n",
	            kappa, preevent_results[[ki]]$p_detect,
	            preevent_results[[ki]]$frac_averted_all))
}

preevent_results <- do.call(rbind, preevent_results)

cat("\n")
print(preevent_results)
cat("\n")

# ==============================================================================
# Section 4: Individual-level — Exit-isolation screening
# ==============================================================================

cat("--- Section 4: Exit-isolation screening ---\n")

exit_results <- list()

for (ki in seq_along(kappas)) {
	kappa <- kappas[ki]

	iso_days <- numeric(n_ind)
	frac_covered <- numeric(n_ind)
	premature <- logical(n_ind)

	for (i in 1:n_ind) {
		ind <- generate_individual(kappa, r, alpha_total, R0, t_inf = 0)

		# Isolation starts at onset: t_inf + s_i = s_i (since t_inf=0)
		t_iso_start <- ind$s_i

		# Must first test positive before entering test-to-exit protocol.
		# Wait until inside positivity window, then require a positive test.
		# Daily testing starting day 1 after isolation.
		day <- 1
		entered_protocol <- FALSE
		released <- FALSE
		while (!released) {
			t_test <- t_iso_start + day

			in_win <- (t_test >= ind$t_peak - window_before) &
				(t_test <= ind$t_peak + window_after)

			if (in_win) {
				pos <- runif(1) < sensitivity
			} else {
				pos <- runif(1) < fp_rate
			}

			if (!entered_protocol) {
				# Need a positive test to enter the test-to-exit protocol
				if (pos) {
					entered_protocol <- TRUE
				}
			} else {
				# Already in protocol: release on first negative
				if (!pos) {
					released <- TRUE
					iso_days[i] <- day
					t_release <- t_test

					# Premature if remaining transmission attempts after release
					if (ind$nattempts > 0) {
						premature[i] <- any(ind$attempt_times > t_release)
						frac_covered[i] <- mean(ind$attempt_times >= t_iso_start &
						                        ind$attempt_times <= t_release)
					} else {
						premature[i] <- FALSE
						frac_covered[i] <- 1
					}
				}
			}

			day <- day + 1
			# Safety valve: don't exceed 60 days
			if (day > 60) {
				iso_days[i] <- 60
				t_release <- t_iso_start + 60
				premature[i] <- FALSE
				if (ind$nattempts > 0) {
					frac_covered[i] <- mean(ind$attempt_times >= t_iso_start &
					                        ind$attempt_times <= t_release)
				} else {
					frac_covered[i] <- 1
				}
				break
			}
		}
	}

	exit_results[[ki]] <- data.frame(
		kappa = kappa,
		mean_iso_days = mean(iso_days),
		sd_iso_days = sd(iso_days),
		frac_covered = mean(frac_covered),
		p_premature = mean(premature)
	)

	cat(sprintf("  kappa=%.1f: mean_iso=%.1f days, frac_covered=%.3f, p_premature=%.3f\n",
	            kappa, exit_results[[ki]]$mean_iso_days,
	            exit_results[[ki]]$frac_covered,
	            exit_results[[ki]]$p_premature))
}

exit_results <- do.call(rbind, exit_results)

cat("\n")
print(exit_results)
cat("\n")

# ==============================================================================
# Section 5: Epidemic simulation with screening
# ==============================================================================

cat("--- Section 5: Epidemic simulation with screening ---\n")

#' Simulate a stochastic epidemic with periodic test-based screening
#'
#' Extends the sim_stochastic_fast pattern with:
#'   - Source tracking for each infection attempt (to check isolation status)
#'   - Periodic screening events that test all infected-not-isolated individuals
#'   - Isolation: detected individuals' future infection attempts are blocked
#'
#' @param n              Population size
#' @param mu             Mean generation time
#' @param R0             Basic reproduction number
#' @param alpha_total    Shape of population kernel
#' @param kappa          Profile punctuation parameter
#' @param screen_interval Screening interval in days (Inf = no screening)
#' @param sensitivity    Test sensitivity within positivity window
#' @param fp_rate        False positive rate outside positivity window
#' @param window_before  Days before peak for positivity window
#' @param window_after   Days after peak for positivity window
#' @return List with tinf, final_size, n_isolated, n_false_iso
sim_epidemic_screening <- function(n, mu, R0, alpha_total, kappa,
                                   screen_interval, sensitivity, fp_rate,
                                   window_before, window_after) {

	r_rate <- alpha_total / mu
	alpha_shift <- alpha_total - kappa

	tinf <- rep(Inf, n)
	t_peak <- rep(Inf, n)
	isolated <- logical(n)
	source_id <- integer(n)

	# Generate infection attempts for a newly infected individual
	gen_attempts <- function(t_inf_val, id) {
		s_i <- rgamma(1, shape = alpha_shift, rate = r_rate)
		nattempts <- rpois(1, R0)
		peak_off <- get_peak_offset(kappa, r_rate)
		t_peak[id] <<- t_inf_val + s_i + peak_off
		if (nattempts == 0L) return(matrix(numeric(0), ncol = 2))
		times <- t_inf_val + s_i + sort(rgamma(nattempts, shape = kappa, rate = r_rate))
		cbind(time = times, source = rep(id, nattempts))
	}

	# Seed infection
	seed <- sample.int(n, 1)
	tinf[seed] <- 0
	new_events <- gen_attempts(0, seed)
	# Queue: 2-col matrix [time, source]
	queue <- new_events
	qi <- 1L

	# Set up screening schedule
	if (is.finite(screen_interval)) {
		# Determine reasonable t_max for screening schedule
		t_max_screen <- qgamma(0.999, shape = alpha_total, rate = r_rate) * 3
		screen_times <- seq(screen_interval, t_max_screen, by = screen_interval)
		si <- 1L  # screen time index
	} else {
		screen_times <- numeric(0)
		si <- 1L
	}

	n_false_iso <- 0L

	# Process events
	while (qi <= nrow(queue)) {
		t_event <- queue[qi, 1]
		src <- as.integer(queue[qi, 2])
		qi <- qi + 1L

		# Fire any screening events that occur before this infection attempt
		while (si <= length(screen_times) && screen_times[si] <= t_event) {
			t_screen <- screen_times[si]
			si <- si + 1L

			# Test all infected, non-isolated individuals
			infected_ids <- which(tinf < Inf & !isolated)
			if (length(infected_ids) > 0) {
				test_results <- run_test_vec(t_screen, t_peak[infected_ids],
				                             sensitivity, fp_rate,
				                             window_before, window_after)
				newly_isolated <- infected_ids[test_results]
				if (length(newly_isolated) > 0) {
					isolated[newly_isolated] <- TRUE
				}
			}

			# False positives on uninfected susceptibles
			n_susceptible <- sum(tinf == Inf)
			if (n_susceptible > 0) {
				n_false_iso <- n_false_iso + rbinom(1, n_susceptible, fp_rate)
			}
		}

		# Process the infection attempt: skip if source is isolated
		if (isolated[src]) next

		target <- sample.int(n, 1)
		if (tinf[target] == Inf) {
			tinf[target] <- t_event
			new_events <- gen_attempts(t_event, target)
			if (nrow(new_events) > 0) {
				remaining <- if (qi <= nrow(queue)) queue[qi:nrow(queue), , drop = FALSE] else matrix(numeric(0), ncol = 2)
				queue <- rbind(remaining, new_events)
				ord <- order(queue[, 1])
				queue <- queue[ord, , drop = FALSE]
				qi <- 1L
			}
		}
	}

	final_size <- sum(tinf < Inf)
	list(
		tinf = tinf,
		final_size = final_size,
		n_isolated = sum(isolated),
		n_false_iso = n_false_iso
	)
}

# --- Run epidemic sweep -------------------------------------------------------

n_pop_epi <- 1000  # population size for epidemic sims

# Baseline (no screening) using sim_stochastic_fast with make_profile_gamma
cat("  Running baseline (no screening) epidemics...\n")
baseline_fs <- list()

for (kappa in kappas) {
	fs_vec <- numeric(nsim_epi)
	prof <- make_profile_gamma(mu = mu, R0 = R0, alpha_total = alpha_total, kappa = kappa)
	for (s in 1:nsim_epi) {
		tinf <- sim_stochastic_fast(n = n_pop_epi, gen_inf_attempts = prof)
		fs_vec[s] <- sum(tinf < Inf)
	}
	baseline_fs[[as.character(kappa)]] <- fs_vec
	cat(sprintf("    kappa=%.1f baseline: mean FS = %.1f\n", kappa, mean(fs_vec)))
}

# Screening epidemics
cat("  Running screening epidemics...\n")
epi_results <- list()
idx <- 0

for (kappa in kappas) {
	for (interval in screen_intervals) {
		fs_vec <- numeric(nsim_epi)
		for (s in 1:nsim_epi) {
			res <- sim_epidemic_screening(
				n = n_pop_epi, mu = mu, R0 = R0, alpha_total = alpha_total,
				kappa = kappa, screen_interval = interval,
				sensitivity = sensitivity, fp_rate = fp_rate,
				window_before = window_before, window_after = window_after
			)
			fs_vec[s] <- res$final_size
		}

		bl <- mean(baseline_fs[[as.character(kappa)]])
		sc <- mean(fs_vec)
		reduction <- ifelse(bl > 1, (bl - sc) / bl * 100, 0)

		idx <- idx + 1
		epi_results[[idx]] <- data.frame(
			kappa = kappa,
			interval = interval,
			mean_fs_baseline = bl,
			mean_fs_screened = sc,
			fs_reduction_pct = reduction
		)

		cat(sprintf("    kappa=%.1f, interval=%d: baseline=%.0f, screened=%.0f, reduction=%.1f%%\n",
		            kappa, interval, bl, sc, reduction))
	}
}

epi_results <- do.call(rbind, epi_results)
epi_results$interval_label <- factor(
	paste0("Every ", epi_results$interval, " day(s)"),
	levels = paste0("Every ", screen_intervals, " day(s)")
)

cat("\n")
print(epi_results)
cat("\n")

# ==============================================================================
# Section 6: Figures
# ==============================================================================

cat("--- Section 6: Generating figures ---\n")

interval_colors <- c("Every 1 day(s)" = "blue",
                     "Every 3 day(s)" = "darkgreen",
                     "Every 7 day(s)" = "red")

# --- Figure 1: Regular screening (detection prob + frac averted vs kappa) -----

fig_regular_A <- ggplot(regular_results,
                        aes(x = kappa, y = p_detect, color = interval_label)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	scale_color_manual(values = interval_colors) +
	labs(x = expression(kappa), y = "Detection probability",
	     color = "Screening interval", tag = "A") +
	theme_classic() +
	theme(legend.position = "bottom")

fig_regular_B <- ggplot(regular_results,
                        aes(x = kappa, y = frac_averted_all, color = interval_label)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	scale_color_manual(values = interval_colors) +
	labs(x = expression(kappa), y = "Fraction of transmission averted\n(unconditional)",
	     color = "Screening interval", tag = "B") +
	theme_classic() +
	theme(legend.position = "bottom")

fig_regular <- cowplot::plot_grid(
	fig_regular_A + theme(legend.position = "none"),
	fig_regular_B + theme(legend.position = "none"),
	nrow = 1
)
legend_regular <- cowplot::get_legend(fig_regular_A)
fig_regular <- cowplot::plot_grid(fig_regular, legend_regular,
                                  ncol = 1, rel_heights = c(1, 0.1))

ggsave("figures/fig_screening_regular.pdf", fig_regular, width = 12, height = 5)
cat("  Saved figures/fig_screening_regular.pdf\n")

# --- Figure 2: Pre-event screening -------------------------------------------

fig_pre_A <- ggplot(preevent_results, aes(x = kappa, y = p_detect)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	labs(x = expression(kappa), y = "Detection probability", tag = "A") +
	theme_classic()

fig_pre_B <- ggplot(preevent_results,
                    aes(x = kappa, y = frac_averted_given_detected)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	labs(x = expression(kappa), y = "Frac. averted | detected", tag = "B") +
	theme_classic()

fig_pre_C <- ggplot(preevent_results, aes(x = kappa, y = frac_averted_all)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	labs(x = expression(kappa), y = "Frac. averted (unconditional)", tag = "C") +
	theme_classic()

fig_preevent <- cowplot::plot_grid(fig_pre_A, fig_pre_B, fig_pre_C, nrow = 1)
ggsave("figures/fig_screening_preevent.pdf", fig_preevent, width = 14, height = 5)
cat("  Saved figures/fig_screening_preevent.pdf\n")

# --- Figure 3: Exit-isolation screening --------------------------------------

fig_exit_A <- ggplot(exit_results, aes(x = kappa, y = mean_iso_days)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	labs(x = expression(kappa), y = "Mean isolation days", tag = "A") +
	theme_classic()

fig_exit_B <- ggplot(exit_results, aes(x = kappa, y = frac_covered)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	labs(x = expression(kappa), y = "Fraction of infectiousness\ncovered by isolation", tag = "B") +
	theme_classic()

fig_exit_C <- ggplot(exit_results, aes(x = kappa, y = p_premature)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	labs(x = expression(kappa), y = "P(premature release)", tag = "C") +
	theme_classic()

fig_exit <- cowplot::plot_grid(fig_exit_A, fig_exit_B, fig_exit_C, nrow = 1)
ggsave("figures/fig_screening_exit.pdf", fig_exit, width = 14, height = 5)
cat("  Saved figures/fig_screening_exit.pdf\n")

# --- Figure 4: Epidemic-level final size reduction ----------------------------

fig_epidemic <- ggplot(epi_results,
                       aes(x = kappa, y = fs_reduction_pct,
                           color = interval_label)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2) +
	scale_color_manual(values = interval_colors) +
	labs(x = expression(kappa),
	     y = "Final size reduction (%)",
	     color = "Screening interval") +
	theme_classic() +
	theme(legend.position = "bottom")

ggsave("figures/fig_screening_epidemic.pdf", fig_epidemic, width = 8, height = 6)
cat("  Saved figures/fig_screening_epidemic.pdf\n")

cat("\n=== Analysis complete ===\n")
