library(tidyverse)
source("code/utils.R")

# ==============================================================================
# Analysis of early epidemic growth: smooth vs. spike
#
# The spike case takes off slower because of a "minimum of order statistics"
# mechanism: in the smooth case, an infector's k infection attempts race independently
# across the generation interval, so the first transmission happens earlier.
# In the spike case, all k infection attempts occur at the same time.
# ==============================================================================

set.seed(123)
n_pop <- 1000
e_dur <- 2 
i_dur <- 3
R0 <- 2

# ==============================================================================
# 1. First-transmission time from a single infector (no epidemic needed)
# ==============================================================================

nrep <- 100000
first_tx <- list()
nattempts <- list()

for(profile in c("smooth", "spike")) {
	ft <- numeric(nrep)
	nc <- numeric(nrep)
	for(i in 1:nrep) {
		k <- rpois(1, R0)
		nc[i] <- k
		if(k == 0) { ft[i] <- Inf; next }
		if(profile == "smooth") {
			ft[i] <- min(rexp(k, 1/e_dur) + rexp(k, 1/i_dur))
		} else {
			ft[i] <- rexp(1, 1/e_dur) + rexp(1, 1/i_dur)
		}
	}
	first_tx[[profile]] <- ft
	nattempts[[profile]] <- nc
}

# Overall comparison (conditional on >= 1 infection attempt)
cat("=== First transmission time from a single infector ===\n")
cat("(conditional on >= 1 infection attempt)\n\n")
for(p in c("smooth", "spike")) {
	vals <- first_tx[[p]][first_tx[[p]] < Inf]
	cat(sprintf("  %s: mean=%.3f, median=%.3f, sd=%.3f (n=%d)\n",
	            p, mean(vals), median(vals), sd(vals), length(vals)))
}
cat(sprintf("\n  Difference in means: %.3f days\n",
            mean(first_tx[["spike"]][first_tx[["spike"]] < Inf]) -
            mean(first_tx[["smooth"]][first_tx[["smooth"]] < Inf])))

# Stratified by number of infection attempts
cat("\n=== First transmission time by number of infection attempts ===\n")
cat(sprintf("  E[GenInterval] = e_dur + i_dur = %.1f days\n\n", e_dur + i_dur))
for(k in 1:6) {
	for(p in c("smooth", "spike")) {
		vals <- first_tx[[p]][nattempts[[p]] == k]
		vals <- vals[vals < Inf]
		if(length(vals) > 10)
			cat(sprintf("  k=%d, %s: %.3f\n", k, p, mean(vals)))
	}
	cat("\n")
}

# ==============================================================================
# 2. Early epidemic dynamics from full simulations
# ==============================================================================

set.seed(42)
nsim <- 5000
max_early <- 50

early_times <- list()
for(profile in c("smooth", "spike")) {
	mat <- matrix(NA, nrow = nsim, ncol = max_early)
	for(i in 1:nsim) {
		tinf <- sim_stochastic_fast(n = n_pop, e_dur = e_dur, i_dur = i_dur,
		                            R0 = R0, profiletype = profile)
		sorted <- sort(tinf[tinf < Inf])
		k <- min(length(sorted), max_early)
		if(k > 0) mat[i, 1:k] <- sorted[1:k]
	}
	early_times[[profile]] <- mat
	cat("Done:", profile, "\n")
}

# Inter-event times for the first 20 infections
cat("\n=== Mean inter-event times (days) ===\n\n")
cat(sprintf("%12s %10s %10s\n", "Interval", "Smooth", "Spike"))
for(j in 2:20) {
	gap_sm <- early_times[["smooth"]][,j] - early_times[["smooth"]][,j-1]
	gap_sp <- early_times[["spike"]][,j]  - early_times[["spike"]][,j-1]
	cat(sprintf("%8d->%d: %10.3f %10.3f\n", j-1, j,
	            mean(gap_sm, na.rm = TRUE), mean(gap_sp, na.rm = TRUE)))
}

# Exponential growth rate (log-linear regression, cases 5-50)
cat("\n=== Empirical exponential growth rate (cases 5-50) ===\n\n")
for(profile in c("smooth", "spike")) {
	mat <- early_times[[profile]]
	has50 <- !is.na(mat[, 50])
	rates <- numeric(sum(has50))
	k <- 0
	for(i in which(has50)) {
		k <- k + 1
		rates[k] <- coef(lm(log(5:50) ~ mat[i, 5:50]))[2]
	}
	rates <- rates[1:k]
	cat(sprintf("  %s: mean=%.4f, sd=%.4f (n=%d)\n",
	            profile, mean(rates), sd(rates), k))
}

# Burstiness: CV of inter-event times in first 20 infections
cat("\n=== Burstiness: CV of inter-event times (first 20 infections) ===\n\n")
for(profile in c("smooth", "spike")) {
	mat <- early_times[[profile]]
	has20 <- !is.na(mat[, 20])
	cvs <- numeric(sum(has20))
	k <- 0
	for(i in which(has20)) {
		gaps <- diff(mat[i, 1:20])
		k <- k + 1
		cvs[k] <- sd(gaps) / mean(gaps)
	}
	cvs <- cvs[1:k]
	cat(sprintf("  %s: mean CV=%.3f, sd=%.3f (n=%d)\n",
	            profile, mean(cvs), sd(cvs), k))
}

# Same-day clustering
cat("\n=== Same-day burst fraction (first 20 infections) ===\n\n")
for(profile in c("smooth", "spike")) {
	mat <- early_times[[profile]]
	has20 <- !is.na(mat[, 20])
	bf <- numeric(sum(has20))
	k <- 0
	for(i in which(has20)) {
		k <- k + 1
		bf[k] <- mean(diff(floor(mat[i, 1:20])) == 0)
	}
	bf <- bf[1:k]
	cat(sprintf("  %s: mean=%.3f, sd=%.3f\n", profile, mean(bf), sd(bf)))
}

# ==============================================================================
# 3. Time-to-threshold analysis across R0 values
#
# Sweep R0 = 1.2, 1.5, 2.0, 3.0, 5.0 and compare smooth vs spike profiles.
# For each (R0, profile) pair, run many simulations, condition on epidemics
# reaching >= 100 infections, and record time to 10 and time to 100 infections.
# ==============================================================================

set.seed(2024)

R0_values <- c(1.2, 1.5, 2.0, 3.0, 5.0)
nsim_r0   <- 3000  # simulations per (R0, profile) pair
threshold_lo  <- 10
threshold_hi  <- 100

cat("\n=== Time-to-threshold analysis across R0 ===\n")
cat(sprintf("  nsim=%d per (R0, profile), N=%d, e_dur=%g, i_dur=%g\n",
            nsim_r0, n_pop, e_dur, i_dur))
cat(sprintf("  Thresholds: %d and %d infections\n", threshold_lo, threshold_hi))
cat(sprintf("  Conditioning on epidemics reaching >= %d infections\n\n",
            threshold_hi))

threshold_results <- list()

for(R0_val in R0_values) {
	for(profile in c("smooth", "spike")) {
		t_lo <- numeric(nsim_r0)
		t_hi <- numeric(nsim_r0)
		reached <- logical(nsim_r0)

		for(i in 1:nsim_r0) {
			tinf <- sim_stochastic_fast(n = n_pop, e_dur = e_dur, i_dur = i_dur,
			                            R0 = R0_val, profiletype = profile)
			sorted <- sort(tinf[tinf < Inf])
			if(length(sorted) >= threshold_hi) {
				reached[i] <- TRUE
				t_lo[i] <- sorted[threshold_lo]
				t_hi[i] <- sorted[threshold_hi]
			}
		}

		threshold_results[[paste(R0_val, profile)]] <- list(
			R0 = R0_val, profile = profile,
			n_established = sum(reached),
			t_lo = t_lo[reached], t_hi = t_hi[reached]
		)
		cat(sprintf("  R0=%.1f, %s: %d/%d reached %d infections\n",
		            R0_val, profile, sum(reached), nsim_r0, threshold_hi))
	}
}

# -- Table 1: Mean time to 10 and 100 infections ------------------------------

cat("\n=== Table: Time to threshold (established epidemics) ===\n\n")
cat(sprintf("%5s | %12s %12s %11s | %13s %13s %12s\n",
            "R0", "t10 smooth", "t10 spike", "delay(t10)",
            "t100 smooth", "t100 spike", "delay(t100)"))
cat(paste(rep("-", 95), collapse = ""), "\n")

for(R0_val in R0_values) {
	sm <- threshold_results[[paste(R0_val, "smooth")]]
	sp <- threshold_results[[paste(R0_val, "spike")]]
	t10_sm <- mean(sm$t_lo); t10_sp <- mean(sp$t_lo)
	t100_sm <- mean(sm$t_hi); t100_sp <- mean(sp$t_hi)
	delay10  <- (t10_sp - t10_sm) / t10_sm * 100
	delay100 <- (t100_sp - t100_sm) / t100_sm * 100
	cat(sprintf("%5.1f | %12.1f %12.1f %+10.1f%% | %13.1f %13.1f %+11.1f%%\n",
	            R0_val, t10_sm, t10_sp, delay10, t100_sm, t100_sp, delay100))
}

# -- Table 2: Variability (SD) in time to 100 infections ----------------------

cat("\n=== Table: Variability in time to 100 infections ===\n\n")
cat(sprintf("%5s | %15s %15s | %18s\n",
            "R0", "SD(t100) smooth", "SD(t100) spike", "spike/smooth ratio"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for(R0_val in R0_values) {
	sm <- threshold_results[[paste(R0_val, "smooth")]]
	sp <- threshold_results[[paste(R0_val, "spike")]]
	sd_sm <- sd(sm$t_hi); sd_sp <- sd(sp$t_hi)
	cat(sprintf("%5.1f | %15.1f %15.1f | %18.2f\n",
	            R0_val, sd_sm, sd_sp, sd_sp / sd_sm))
}
