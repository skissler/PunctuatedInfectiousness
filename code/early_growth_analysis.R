library(tidyverse)
source("code/utils.R")

# ==============================================================================
# Analysis of early epidemic growth: smooth vs. spike
#
# The spike case takes off slower because of a "minimum of order statistics"
# mechanism: in the smooth case, an infector's k contacts race independently
# across the generation interval, so the first transmission happens earlier.
# In the spike case, all k contacts occur at the same time.
# ==============================================================================

set.seed(123)
n_pop <- 1000
e_dur <- 2; i_dur <- 3; R0 <- 2

# ==============================================================================
# 1. First-transmission time from a single infector (no epidemic needed)
# ==============================================================================

nrep <- 100000
first_tx <- list()
ncontacts <- list()

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
	ncontacts[[profile]] <- nc
}

# Overall comparison (conditional on >= 1 contact)
cat("=== First transmission time from a single infector ===\n")
cat("(conditional on >= 1 contact)\n\n")
for(p in c("smooth", "spike")) {
	vals <- first_tx[[p]][first_tx[[p]] < Inf]
	cat(sprintf("  %s: mean=%.3f, median=%.3f, sd=%.3f (n=%d)\n",
	            p, mean(vals), median(vals), sd(vals), length(vals)))
}
cat(sprintf("\n  Difference in means: %.3f days\n",
            mean(first_tx[["spike"]][first_tx[["spike"]] < Inf]) -
            mean(first_tx[["smooth"]][first_tx[["smooth"]] < Inf])))

# Stratified by number of contacts
cat("\n=== First transmission time by number of contacts ===\n")
cat(sprintf("  E[GenInterval] = e_dur + i_dur = %.1f days\n\n", e_dur + i_dur))
for(k in 1:6) {
	for(p in c("smooth", "spike")) {
		vals <- first_tx[[p]][ncontacts[[p]] == k]
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
