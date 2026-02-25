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

# Test that mixture profile with w=0 matches smooth and w=1 matches spike,
# and that the gen_inf_attempts interface produces results equivalent to profiletype.
test_mixture <- function(nsim = 500, n = 500, e_dur = 2, i_dur = 3, R0 = 2) {

	all_passed <- TRUE

	# --- Factory-based profiles should match built-in profiletype ---
	cat("Factory vs built-in equivalence:\n")
	factory_tests <- list(
		list(name = "stepwise", factory = make_profile_stepwise(e_dur, i_dur, R0)),
		list(name = "smooth",   factory = make_profile_smooth(e_dur, i_dur, R0)),
		list(name = "spike",    factory = make_profile_spike(e_dur, i_dur, R0))
	)
	for (ft in factory_tests) {
		fs_builtin <- numeric(nsim)
		fs_factory <- numeric(nsim)
		for (i in 1:nsim) {
			r1 <- sim_stochastic_fast(n=n, e_dur=e_dur, i_dur=i_dur,
			                          R0=R0, profiletype=ft$name)
			fs_builtin[i] <- sum(r1 < Inf)
			r2 <- sim_stochastic_fast(n=n, gen_inf_attempts=ft$factory)
			fs_factory[i] <- sum(r2 < Inf)
		}
		ks <- ks.test(fs_builtin, fs_factory)
		passed <- ks$p.value > 0.01
		status <- if (passed) "PASS" else "FAIL"
		cat(sprintf("  %s: factory %s — KS p=%.4f (mean builtin=%.1f, factory=%.1f)\n",
		            status, ft$name, ks$p.value, mean(fs_builtin), mean(fs_factory)))
		if (!passed) all_passed <- FALSE
	}

	# --- Mixture w=0 should match smooth ---
	cat("\nMixture boundary tests:\n")
	mix_smooth <- make_profile_mixture(e_dur, i_dur, R0, w = 0)
	fs_smooth  <- numeric(nsim)
	fs_mix0    <- numeric(nsim)
	for (i in 1:nsim) {
		r1 <- sim_stochastic_fast(n=n, e_dur=e_dur, i_dur=i_dur,
		                          R0=R0, profiletype="smooth")
		fs_smooth[i] <- sum(r1 < Inf)
		r2 <- sim_stochastic_fast(n=n, gen_inf_attempts=mix_smooth)
		fs_mix0[i] <- sum(r2 < Inf)
	}
	ks <- ks.test(fs_smooth, fs_mix0)
	passed <- ks$p.value > 0.01
	status <- if (passed) "PASS" else "FAIL"
	cat(sprintf("  %s: mixture(w=0) vs smooth — KS p=%.4f (mean smooth=%.1f, mix0=%.1f)\n",
	            status, ks$p.value, mean(fs_smooth), mean(fs_mix0)))
	if (!passed) all_passed <- FALSE

	# --- Mixture w=1 should match spike ---
	mix_spike <- make_profile_mixture(e_dur, i_dur, R0, w = 1)
	fs_spike  <- numeric(nsim)
	fs_mix1   <- numeric(nsim)
	for (i in 1:nsim) {
		r1 <- sim_stochastic_fast(n=n, e_dur=e_dur, i_dur=i_dur,
		                          R0=R0, profiletype="spike")
		fs_spike[i] <- sum(r1 < Inf)
		r2 <- sim_stochastic_fast(n=n, gen_inf_attempts=mix_spike)
		fs_mix1[i] <- sum(r2 < Inf)
	}
	ks <- ks.test(fs_spike, fs_mix1)
	passed <- ks$p.value > 0.01
	status <- if (passed) "PASS" else "FAIL"
	cat(sprintf("  %s: mixture(w=1) vs spike  — KS p=%.4f (mean spike=%.1f, mix1=%.1f)\n",
	            status, ks$p.value, mean(fs_spike), mean(fs_mix1)))
	if (!passed) all_passed <- FALSE

	# --- Smoke test: w=0.5 runs without error ---
	mix_half <- make_profile_mixture(e_dur, i_dur, R0, w = 0.5)
	r <- sim_stochastic_fast(n = n, gen_inf_attempts = mix_half)
	ninfected <- sum(r < Inf)
	cat(sprintf("\n  Smoke test: mixture(w=0.5) — %d/%d infected\n", ninfected, n))

	cat(sprintf("\n%s\n", if(all_passed) "All mixture tests passed." else "Some mixture tests FAILED."))
	invisible(all_passed)
}

test_equivalence()
cat("\n")
test_mixture()
