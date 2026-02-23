source("code/utils.R")

# Statistical equivalence test: runs many simulations with each function
# and compares the distributions of key summary statistics via KS tests.
# (The two functions use different RNG internals so identical-seed comparison
# is not possible, but their output distributions should be indistinguishable.)
test_equivalence <- function(nsim = 500, n = 500, e_dur = 2, i_dur = 3, R0 = 2) {

	profiles <- c("stepwise", "smooth", "spike")
	all_passed <- TRUE

	for(profile in profiles) {
		fs_orig <- numeric(nsim)
		fs_fast <- numeric(nsim)
		for(i in 1:nsim) {
			r1 <- sim_stochastic(n=n, e_dur=e_dur, i_dur=i_dur,
			                     R0=R0, profiletype=profile)
			fs_orig[i] <- sum(r1 < Inf)

			r2 <- sim_stochastic_fast(n=n, e_dur=e_dur, i_dur=i_dur,
			                          R0=R0, profiletype=profile)
			fs_fast[i] <- sum(r2 < Inf)
		}

		ks <- ks.test(fs_orig, fs_fast)
		passed <- ks$p.value > 0.01
		status <- if(passed) "PASS" else "FAIL"
		cat(sprintf("  %s: %s — final size KS p=%.4f (mean orig=%.1f, fast=%.1f)\n",
		            status, profile, ks$p.value, mean(fs_orig), mean(fs_fast)))
		if(!passed) all_passed <- FALSE
	}

	cat(sprintf("\n%s\n", if(all_passed) "All tests passed." else "Some tests FAILED."))
	invisible(all_passed)
}

test_equivalence()
