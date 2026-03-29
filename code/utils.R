library(tidyverse) 
library(odin) 

# ==============================================================================
# Utility functions
# 
# Defines functions for simulating deterministic and stochastic epidemics, with
# helper functions 
# ==============================================================================

# ==============================================================================
# Deterministic models 
# ==============================================================================

#' SEIR ODE model
#'
#' Generates a standard SEIR epidemic trajectory, initialized with a single
#' person in the infectious (I) state. 
#'
#' @param R0    Basic reproduction number (default 2)
#' @param nu    E->I rate, 1/mean latent period (default 1/2)
#' @param gamma I->R rate, 1/mean infectious period (default 1/3)
#' @param N     Population size, used to set initial conditions (default 1000)
seir <- odin({
  ## Derivatives
  deriv(S) <- -beta*S*I
  deriv(E) <- beta*S*I - nu*E
  deriv(I) <- nu*E - gamma*I
  deriv(R) <- gamma*I
  deriv(cuminf) <- beta*S*I

  ## Initial conditions
  initial(S) <- 1-1/N
  initial(E) <- 0
  initial(I) <- 1/N
  initial(R) <- 0
  initial(cuminf) <- 1/N

  ## Parameters with default values 
  R0 <- user(2)
  nu <- user(1/2)
  gamma <- user(1/3)
  beta <- R0*gamma
  N <- user(1000)
})

#' Renewal equation epidemic model
#'
#' Generates a deterministic epidemic trajectory using a discrete approximation
#' to the renewal equation framework, with a Gamma-distributed infectiousness 
#' profile and initalized with a single infectious person. 
#'
#' @param R0     Basic reproduction number
#' @param alpha  Shape parameter of Gamma generation interval distribution
#' @param T_gen  Mean generation interval (days)
#' @param N      Population size, used to set initial conditions
#' @param dt     Time step for discretization (default 0.01 days)
#' @param tmax   Maximum simulation length (default 100 days)

renewal_epidemic <- function(R0, alpha, T_gen, N, dt = 0.01, tmax = 100) {

	# Set rate parameter of the generation interval distribution
  beta <- alpha / T_gen

  # Initalize time vectors
  times <- seq(0, tmax, by = dt)
  nt <- length(times)

  # Initialize incidence and cumulative infections
  j <- numeric(nt)
  C <- numeric(nt)
  C[1] <- 1 / N

  # Generation interval weights (with dt baked in for the convolution)
  g <- dgamma(times, shape = alpha, rate = beta) * dt
  # Index case kernel (no dt — it's a single person, not a rate)
  g_index <- dgamma(times, shape = alpha, rate = beta) / N

  for (i in 2:nt) {
    S <- 1 - C[i-1]
    # Index case + renewal convolution
    j[i] <- S * R0 * (g_index[i] + sum(g[2:i] * j[(i-1):1]))
    C[i] <- C[i-1] + j[i] * dt
  }

  tibble(t = times, j = j, cuminf = C, S = 1 - C)
}

# ==============================================================================
# Stochastic simulation model
# ==============================================================================

gen_inf_attempts_stepwise <- function(e_dur, i_dur, R0){
	beta <- R0/i_dur
	function(tinf){
		t1 <- rexp(1, 1/e_dur)
		t2 <- t1 + rexp(1, 1/i_dur)
		nattempts <- rpois(1, (t2 - t1)*beta) 
		if(nattempts == 0L) return(numeric(0))
		attempt_times <- tinf + sort(runif(nattempts, min=t1, max=t2))
		return(attempt_times)
	}
}

gen_inf_attempts_smooth <- function(e_dur, i_dur, R0){
	function(tinf){
		nattempts <- rpois(1,R0)
		if(nattempts == 0L) return(numeric(0))
		attempt_times <- tinf + sort(rexp(nattempts, 1/e_dur) + rexp(nattempts, 1/i_dur))
		return(attempt_times)
	}
}

gen_inf_attempts_spike <- function(e_dur, i_dur, R0){
	function(tinf){
		nattempts <- rpois(1,R0)
		if(nattempts == 0L) return(numeric(0))
		attempt_times <- tinf + rep(rexp(1, 1/e_dur) + rexp(1, 1/i_dur), nattempts)
		return(attempt_times)
	}
}

#' Infection attempt generator for gamma profiles
#'
#' Generates attempted infection times for a single person, themselves infected
#' at time tinf, using the Gamma-distributed individual infectiousness profile
#' with punctuation parameter psi.
#'
#' @param T Mean generation interval
#' @param R0 Basic reproduction number
#' @param alpha Shape parameter for the Gamma-distributed population
#'   infectiousness profile (A(tau) ~ Gamma(alpha, beta) where beta = alpha/T)
#' @param psi Parameter governing punctuation of the individual
#'   infectiousness profile (psi \in [0, 1], with psi -> 0 giving sharp
#'   infectiousness profiles and psi -> 1 giving smooth ones equivalent to
#'   the population generation interval distribution)
#' @return A function(t_inf) that returns sorted infection attempt times
gen_inf_attempts_gamma <- function(T, R0, alpha, psi){
	stopifnot(psi>=0, psi<=1)
	beta <- alpha/T

	if(psi < 1e-6){ # spike implementation
		shiftshape <- alpha
		function(tinf){
			nattempts <- rpois(1,R0)
			if(nattempts==0L) return(numeric(0))
			shift <- rgamma(1, shape=shiftshape, rate=beta)
			attempt_times <- rep(tinf + shift, nattempts)
			return(sort(attempt_times))
		}
	} else if(psi > 1-1e-6){ # smooth implementation
		function(tinf){
			nattempts <- rpois(1,R0)
			if(nattempts==0L) return(numeric(0))
			attempt_times <- tinf + rgamma(nattempts, shape=alpha, rate=beta)
			return(sort(attempt_times))
		}
	} else { # regular implementation
		shiftshape <- (1-psi)*alpha
		function(tinf){
			nattempts <- rpois(1,R0)
			if(nattempts==0L) return(numeric(0))
			shift <- rgamma(1, shape=shiftshape, rate=beta)
			attempt_times <- tinf + shift + rgamma(nattempts, shape=psi*alpha, rate=beta)
			return(sort(attempt_times))
		}
	}

}

gen_inf_attempts_gamma_contacts <- function(T, z, z_max, alpha, psi){

	stopifnot(psi>0, psi<1)
	stopifnot(z_max>0)

	beta <- alpha/T
	shiftshape <- (1-psi)*alpha

	function(tinf){
		nproposals <- rpois(1, z_max)
		if(nproposals == 0L) return(numeric(0))
		shift <- rgamma(1, shape=shiftshape, rate=beta)
		proposals <- tinf + shift + rgamma(nproposals, shape=psi*alpha, rate=beta)
		accept_prob <- z(proposals)/z_max
		keep <- runif(nproposals) < accept_prob
		if(!any(keep)) return(numeric(0))
		return(sort(proposals[keep]))
	}

}

# --- Helpers for detection-and-isolation protocols ---

# Determine isolation time under regular screening.
# mode_time: absolute time of the infectiousness profile mode.
# Returns Inf if undetected.
.detect_screening <- function(mode_time, d_pre, w, Delta, p_sens, lambda_act){
	window_start <- mode_time - d_pre
	window_end <- window_start + w
	U <- runif(1, 0, Delta)
	first_test <- window_start + U
	if(first_test > window_end) return(Inf)
	test_times <- seq(first_test, window_end, by=Delta)
	if(p_sens >= 1){
		tau_det <- test_times[1]
	} else {
		detected <- runif(length(test_times)) < p_sens
		if(!any(detected)) return(Inf)
		tau_det <- test_times[which(detected)[1]]
	}
	tau_det + rexp(1, lambda_act)
}

# Determine isolation time under symptom-triggered isolation.
# mode_time: absolute time of the infectiousness profile mode.
# Returns Inf if asymptomatic or (implicitly) if delay is infinite.
.detect_symptoms <- function(mode_time, mu_sym, sigma_sym, lambda_act, p_sym){
	if(runif(1) > p_sym) return(Inf)
	delta_sym <- if(sigma_sym > 0) rnorm(1, mu_sym, sigma_sym) else mu_sym
	mode_time + delta_sym + rexp(1, lambda_act)
}

# Remove (or thin with probability eta) infection attempts after isolation.
# Returns sorted surviving attempt times.
.thin_attempts <- function(attempt_times, tau_iso, eta){
	if(!is.finite(tau_iso) || length(attempt_times) == 0L) return(attempt_times)
	pre <- attempt_times <= tau_iso
	if(eta >= 1){
		attempt_times <- attempt_times[pre]
	} else {
		npost <- sum(!pre)
		if(npost > 0L){
			post_keep <- runif(npost) >= eta
			attempt_times <- c(attempt_times[pre], attempt_times[!pre][post_keep])
		}
	}
	if(length(attempt_times) == 0L) return(numeric(0))
	sort(attempt_times)
}

#' Infection attempt generator with fixed isolation time
#'
#' Like gen_inf_attempts_gamma, but isolates every individual at a fixed
#' time offset relative to the mode of their infectiousness profile.
#'
#' @param T Mean generation interval
#' @param R0 Basic reproduction number
#' @param alpha Shape parameter for population infectiousness profile
#' @param psi Punctuation parameter (0 = spike, 1 = smooth)
#' @param tau_offset Isolation time relative to peak infectiousness (days);
#'   negative = before peak, positive = after peak
#' @param p_adhere Probability that an individual adheres to isolation
#' @param eta Isolation effectiveness (1 = perfect, 0 = none)
#' @return A function(t_inf) that returns sorted infection attempt times
gen_inf_attempts_gamma_fixed_iso <- function(T, R0, alpha, psi, tau_offset,
                                             p_adhere=1, eta=1){
	stopifnot(psi>=0, psi<=1)
	beta <- alpha/T

	if(psi < 1e-6){ # spike implementation
		shiftshape <- alpha
		function(tinf){
			nattempts <- rpois(1, R0)
			if(nattempts==0L) return(numeric(0))
			shift <- rgamma(1, shape=shiftshape, rate=beta)
			attempt_times <- rep(tinf + shift, nattempts)
			tau_iso <- if(runif(1) < p_adhere) tinf + shift + tau_offset else Inf
			.thin_attempts(attempt_times, tau_iso, eta)
		}
	} else if(psi > 1-1e-6){ # smooth implementation
		mode_psi <- max(0, (alpha-1)/beta)
		function(tinf){
			nattempts <- rpois(1, R0)
			if(nattempts==0L) return(numeric(0))
			attempt_times <- tinf + rgamma(nattempts, shape=alpha, rate=beta)
			tau_iso <- if(runif(1) < p_adhere) tinf + mode_psi + tau_offset else Inf
			.thin_attempts(attempt_times, tau_iso, eta)
		}
	} else { # regular implementation
		shiftshape <- (1-psi)*alpha
		mode_psi <- max(0, (psi*alpha-1)/beta)
		function(tinf){
			nattempts <- rpois(1, R0)
			if(nattempts==0L) return(numeric(0))
			shift <- rgamma(1, shape=shiftshape, rate=beta)
			attempt_times <- tinf + shift + rgamma(nattempts, shape=psi*alpha, rate=beta)
			tau_iso <- if(runif(1) < p_adhere) tinf + shift + mode_psi + tau_offset else Inf
			.thin_attempts(attempt_times, tau_iso, eta)
		}
	}
}

#' Infection attempt generator with regular screening and isolation
#'
#' Like gen_inf_attempts_gamma, but after generating attempts, determines
#' whether the individual is detected via regular screening (tests every
#' Delta days with random phase, detection within a window around peak
#' infectiousness) and removes post-isolation attempts.
#'
#' @param T Mean generation interval
#' @param R0 Basic reproduction number
#' @param alpha Shape parameter for population infectiousness profile
#' @param psi Punctuation parameter (0 = spike, 1 = smooth)
#' @param d_pre Days before peak that detectability window begins
#' @param d_post Days after peak that detectability window ends
#' @param Delta Gap between tests (days)
#' @param lambda_act Rate of exponential delay from detection to isolation
#'   (Inf = immediate isolation, the default)
#' @param p_sens Test sensitivity (probability of detection per test)
#' @param eta Isolation effectiveness (1 = perfect, 0 = none)
#' @return A function(t_inf) that returns sorted infection attempt times
gen_inf_attempts_gamma_screening <- function(T, R0, alpha, psi,
                                             d_pre, d_post, Delta,
                                             lambda_act=Inf, p_sens=1, eta=1){
	stopifnot(psi>=0, psi<=1)
	beta <- alpha/T
	w <- d_pre + d_post

	if(psi < 1e-6){ # spike implementation
		shiftshape <- alpha
		function(tinf){
			nattempts <- rpois(1, R0)
			if(nattempts==0L) return(numeric(0))
			shift <- rgamma(1, shape=shiftshape, rate=beta)
			attempt_times <- rep(tinf + shift, nattempts)
			mode_time <- tinf + shift
			tau_iso <- .detect_screening(mode_time, d_pre, w, Delta, p_sens, lambda_act)
			.thin_attempts(attempt_times, tau_iso, eta)
		}
	} else if(psi > 1-1e-6){ # smooth implementation
		mode_psi <- max(0, (alpha-1)/beta)
		function(tinf){
			nattempts <- rpois(1, R0)
			if(nattempts==0L) return(numeric(0))
			attempt_times <- tinf + rgamma(nattempts, shape=alpha, rate=beta)
			mode_time <- tinf + mode_psi
			tau_iso <- .detect_screening(mode_time, d_pre, w, Delta, p_sens, lambda_act)
			.thin_attempts(attempt_times, tau_iso, eta)
		}
	} else { # regular implementation
		shiftshape <- (1-psi)*alpha
		mode_psi <- max(0, (psi*alpha-1)/beta)
		function(tinf){
			nattempts <- rpois(1, R0)
			if(nattempts==0L) return(numeric(0))
			shift <- rgamma(1, shape=shiftshape, rate=beta)
			attempt_times <- tinf + shift + rgamma(nattempts, shape=psi*alpha, rate=beta)
			mode_time <- tinf + shift + mode_psi
			tau_iso <- .detect_screening(mode_time, d_pre, w, Delta, p_sens, lambda_act)
			.thin_attempts(attempt_times, tau_iso, eta)
		}
	}
}

#' Infection attempt generator with symptom-triggered isolation
#'
#' Like gen_inf_attempts_gamma, but after generating attempts, determines
#' whether the individual develops symptoms (probability p_sym) at a time
#' anchored to peak infectiousness, and removes post-isolation attempts.
#'
#' @param T Mean generation interval
#' @param R0 Basic reproduction number
#' @param alpha Shape parameter for population infectiousness profile
#' @param psi Punctuation parameter (0 = spike, 1 = smooth)
#' @param mu_sym Mean symptom onset time relative to peak infectiousness
#'   (positive = symptoms lag peak; negative = symptoms precede peak)
#' @param sigma_sym SD of symptom onset time (0 = deterministic, the default)
#' @param lambda_act Rate of exponential delay from symptom onset to isolation
#'   (Inf = immediate isolation, the default)
#' @param p_sym Probability of being symptomatic
#' @param eta Isolation effectiveness (1 = perfect, 0 = none)
#' @return A function(t_inf) that returns sorted infection attempt times
gen_inf_attempts_gamma_symptoms <- function(T, R0, alpha, psi,
                                            mu_sym, sigma_sym=0,
                                            lambda_act=Inf, p_sym=1, eta=1){
	stopifnot(psi>=0, psi<=1)
	beta <- alpha/T

	if(psi < 1e-6){ # spike implementation
		shiftshape <- alpha
		function(tinf){
			nattempts <- rpois(1, R0)
			if(nattempts==0L) return(numeric(0))
			shift <- rgamma(1, shape=shiftshape, rate=beta)
			attempt_times <- rep(tinf + shift, nattempts)
			mode_time <- tinf + shift
			tau_iso <- .detect_symptoms(mode_time, mu_sym, sigma_sym, lambda_act, p_sym)
			.thin_attempts(attempt_times, tau_iso, eta)
		}
	} else if(psi > 1-1e-6){ # smooth implementation
		mode_psi <- max(0, (alpha-1)/beta)
		function(tinf){
			nattempts <- rpois(1, R0)
			if(nattempts==0L) return(numeric(0))
			attempt_times <- tinf + rgamma(nattempts, shape=alpha, rate=beta)
			mode_time <- tinf + mode_psi
			tau_iso <- .detect_symptoms(mode_time, mu_sym, sigma_sym, lambda_act, p_sym)
			.thin_attempts(attempt_times, tau_iso, eta)
		}
	} else { # regular implementation
		shiftshape <- (1-psi)*alpha
		mode_psi <- max(0, (psi*alpha-1)/beta)
		function(tinf){
			nattempts <- rpois(1, R0)
			if(nattempts==0L) return(numeric(0))
			shift <- rgamma(1, shape=shiftshape, rate=beta)
			attempt_times <- tinf + shift + rgamma(nattempts, shape=psi*alpha, rate=beta)
			mode_time <- tinf + shift + mode_psi
			tau_iso <- .detect_symptoms(mode_time, mu_sym, sigma_sym, lambda_act, p_sym)
			.thin_attempts(attempt_times, tau_iso, eta)
		}
	}
}

#' Stochastic epidemic simulation
#'
#' Generates the infection times for a population (size n) given individual
#' infectiousness profiles specified by gen_inf_attempts
#'
#' @param n Population size (default 1000)
#' @param gen_inf_attempts Infection attempt time generator function, 
#'   capturing the individual infectiousness profile, with signature 
#'   function(t_inf) -> numeric vector. 
#' @return Numeric vector of length n containing infection times (Inf if the 
#'   person remains uninfected at the end of the simulation)
sim_stochastic_fast <- function(n=1000, gen_inf_attempts){

	# Initialize infection times 
	tinf_vec <- rep(Inf, n)

	# Seed infection
	indexcase <- sample.int(n,1)
	tinf_vec[indexcase] <- 0
	queue <- gen_inf_attempts(0)
	qi <- 1L
	qn <- length(queue)

	# Process events in chronological order 
	while(qi <= qn){
		t_attempt <- queue[qi]
		qi <- qi + 1L 

		target <- sample.int(n,1)
		if(tinf_vec[target]==Inf){
			tinf_vec[target] <- t_attempt
			new_attempts <- gen_inf_attempts(t_attempt)
			if(length(new_attempts) > 0L){
				# Merge new events with unprocessed remainder of queue 
				remaining <- if(qi <= qn) queue[qi:qn] else numeric(0)
				queue <- sort.int(c(remaining, new_attempts))
				qi <- 1L
				qn <- length(queue)
			}
		}
	}
	return(tinf_vec)
}

# Helper: save a ggplot as both .png and .pdf in figures/
save_fig <- function(plot, name, width = 10, height = 5) {
	ggsave(file.path("figures", paste0(name, ".pdf")), plot, width = width, height = height)
	ggsave(file.path("figures", paste0(name, ".png")), plot, width = width, height = height, dpi = 300)
}

# Helper: construct cache file path encoding nsim and popsize
cache_path <- function(pathogen, nsim, popsize) {
	file.path("output", sprintf("cuminf_df_%s_n%d_s%d.csv", pathogen, popsize, nsim))
}

# Helper: load cached simulations, or return NULL if cache is missing/stale
#
# Returns the cached data frame if the file exists and contains all required
# psi values. Returns NULL otherwise (caller should re-run simulations).
load_cache <- function(pathogen, nsim, popsize, psivals) {
	cache_file <- cache_path(pathogen, nsim, popsize)
	if (!file.exists(cache_file)) return(NULL)

	df <- read_csv(cache_file, show_col_types = FALSE)
	cached_psi <- sort(unique(df$psi))
	if (!all(psivals %in% cached_psi)) {
		cat(sprintf("  %s: cache missing psi values %s, re-running\n",
		    pathogen, paste(setdiff(psivals, cached_psi), collapse = ", ")))
		return(NULL)
	}

	cat(sprintf("  %s: loading cached simulations from %s\n", pathogen, cache_file))
	df %>% mutate(psi = factor(psi))
}

