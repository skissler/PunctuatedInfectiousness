library(tidyverse)
source("code/utils.R")
source("code/parameters.R")

# ==============================================================================
# Diagnose bias in compute_growth_rate() and explore threshold choices
# ==============================================================================

set.seed(2026)

# ---- Pure NHPP simulator (rate lambda0 * exp(r*t)) -------------------------
sim_nhpp <- function(K, r, lambda0 = 1) {
	t <- numeric(K)
	s <- 0
	for (k in seq_len(K)) {
		U <- runif(1)
		x <- (1 / r) * log(1 - r * log(1 - U) / (lambda0 * exp(r * s)))
		s <- s + x
		t[k] <- s
	}
	t
}

n_reps <- 2000

# ---- (A) Verify fix on pure NHPP at r = 0.25 and r = 0.38 ------------------
cat("==== (A) Estimator verification on pure NHPP ====\n")
for (r_true in c(0.25, 0.38)) {
	rhat <- replicate(n_reps, {
		t <- sim_nhpp(K = 510, r = r_true)
		compute_growth_rate(t, 100, 500)
	})
	cat(sprintf("  r_true = %.2f: mean r_hat = %.4f (bias = %+.1f%%)\n",
	            r_true, mean(rhat), 100 * (mean(rhat)/r_true - 1)))
}

# ---- (B) Threshold sweep under saturation ----------------------------------
# For each pathogen, compute the theoretical effective r at cumulative =
# growth_threshold (i.e., the "instantaneous" r at the right edge of the fit
# window) and averaged over the window via midpoint.
cat("\n==== (B) Residual saturation bias (theoretical) for N = 10000 ====\n")
cat("    Effective r at the mid-point of the fit window, as % of true r:\n\n")
N <- 10000

thresholds <- list(
	c(min = 10,   max = 100),
	c(min = 50,   max = 200),
	c(min = 50,   max = 300),
	c(min = 100,  max = 500),
	c(min = 200,  max = 1000)
)

cat(sprintf("%-10s  %-18s  %-18s  %-18s  %-18s  %-18s\n",
            "pathogen", "10-100", "50-200", "50-300", "100-500", "200-1000"))

for (p in parslist) {
	r0 <- p$r
	vals <- sapply(thresholds, function(th) {
		mid <- (th["min"] + th["max"]) / 2
		S <- 1 - mid / N
		r_sat <- p$beta * ((p$R0 * S)^(1 / p$alpha) - 1)
		100 * (r_sat / r0 - 1)
	})
	cat(sprintf("%-10s  %+7.1f%%           %+7.1f%%           %+7.1f%%           %+7.1f%%           %+7.1f%%\n",
	            p$pathogen, vals[1], vals[2], vals[3], vals[4], vals[5]))
}

# ---- (C) Full empirical test on branching simulations ----------------------
# Run actual simulations for each pathogen with different window thresholds.
# This captures saturation AND residual estimator bias together.
cat("\n==== (C) Empirical bias: full simulation (2000 sims) ====\n\n")

n_sims_per_pathogen <- 500

set.seed(42)
for (p in parslist) {
	pathogen <- p$pathogen
	T <- p$Tgen; alpha <- p$alpha; beta <- p$beta; R0 <- p$R0
	r_true <- p$r

	# Single gamma generator (psi=1, smooth case) for speed — pure renewal
	gfun <- gen_inf_attempts_gamma(T, R0, alpha, psi = 1)

	cat(sprintf("  Running %s (r_true = %.3f) ...\n", pathogen, r_true))

	results <- matrix(NA_real_, n_sims_per_pathogen, length(thresholds))
	colnames(results) <- sapply(thresholds, function(th) sprintf("%d-%d", th["min"], th["max"]))

	for (s in seq_len(n_sims_per_pathogen)) {
		tinf <- sim_stochastic_fast(n = N, gen_inf_attempts = gfun)
		infection_times <- sort(tinf[is.finite(tinf)])
		if (length(infection_times) < 0.1 * N) next  # not established
		for (j in seq_along(thresholds)) {
			th <- thresholds[[j]]
			results[s, j] <- compute_growth_rate(infection_times, th["min"], th["max"])
		}
	}

	means <- colMeans(results, na.rm = TRUE)
	sds   <- apply(results, 2, sd, na.rm = TRUE)

	cat(sprintf("    %-10s  ", "window"))
	cat(paste(sprintf("%-14s", colnames(results)), collapse = ""), "\n")
	cat(sprintf("    %-10s  ", "mean r_hat"))
	cat(paste(sprintf("%-14.4f", means), collapse = ""), "\n")
	cat(sprintf("    %-10s  ", "bias (%)"))
	cat(paste(sprintf("%-14s", sprintf("%+.1f%%", 100*(means/r_true - 1))), collapse = ""), "\n")
	cat(sprintf("    %-10s  ", "SD"))
	cat(paste(sprintf("%-14.4f", sds), collapse = ""), "\n\n")
}
