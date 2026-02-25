library(tidyverse)
library(odin)

# ==============================================================================
# Mean-field ODE models (proportions, not counts)
# ==============================================================================

#' SIR mean-field ODE model (odin)
#'
#' State variables are proportions: S + I + R = 1.
#' Transmission rate beta = R0 * gamma so that R0 = beta / gamma.
#' A single seed individual (1/N) starts in the I compartment.
#'
#' @param R0    Basic reproduction number (default 2)
#' @param gamma Recovery rate, 1/mean infectious period (default 1/5)
#' @param N     Population size, used only to set initial conditions (default 1000)
sir <- odin({
  ## Derivatives
  deriv(S) <- -beta*S*I
  deriv(I) <- beta*S*I - gamma*I
  deriv(R) <- gamma*I
  deriv(cuminf) <- beta*S*I

  ## Initial conditions
  initial(S) <- 1-1/N
  initial(I) <- 1/N
  initial(R) <- 0
  initial(cuminf) <- 1/N

  ## parameters
  R0 <- user(2)
  gamma <- user(1/5)
  beta <- R0*gamma
  N <- user(1000)
})

#' SEIR mean-field ODE model (odin)
#'
#' Adds an exposed (E) compartment with latency rate nu = 1/mean latent period.
#' R0 = beta / gamma (the E compartment adds delay but does not change R0).
#'
#' @param R0    Basic reproduction number (default 2)
#' @param nu    Latency rate, 1/mean latent period (default 1/2)
#' @param gamma Recovery rate, 1/mean infectious period (default 1/3)
#' @param N     Population size, used only to set initial conditions (default 1000)
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

  ## parameters
  R0 <- user(2)
  nu <- user(1/2)
  gamma <- user(1/3)
  beta <- R0*gamma
  N <- user(1000)
})

# ==============================================================================
# Stochastic epidemic simulation
# ==============================================================================

#' Simulate a stochastic epidemic under a specified individual infectiousness profile
#'
#' Runs an event-driven stochastic simulation of an SIR/SEIR-like epidemic.
#' Each infection attempt targets a uniformly random individual; if the target is
#' still susceptible it becomes infected, otherwise the infection attempt is wasted.
#' This gives frequency-dependent (mass-action) transmission consistent with
#' the ODE models above.
#'
#' Three individual-level infectiousness profiles are available, all of which
#' converge to the same population-level SEIR infectiousness kernel A(tau) in
#' expectation (see notes/overview.md):
#'
#'   "stepwise" — infectious at constant rate beta during a random window
#'     [onset, onset + duration], where onset ~ Exp(1/e_dur) and
#'     duration ~ Exp(1/i_dur). Number of infection attempts is Poisson(beta * duration),
#'     so the secondary infection count is overdispersed relative to Poisson(R0).
#'
#'   "smooth" — Poisson(R0) infection attempts, each independently timed according to
#'     the SEIR generation interval (Exp(1/e_dur) + Exp(1/i_dur)). Every
#'     individual has the same expected infectiousness profile; the only
#'     stochasticity is Poisson variation in the infection attempt count.
#'
#'   "spike" — Poisson(R0) infection attempts all occurring at a single random time
#'     drawn from the SEIR generation interval. Maximally punctuated:
#'     individual infectiousness is a (scaled) delta function.
#'
#' @param n           Population size (default 1000)
#' @param e_dur       Mean latent period; 0 gives SIR dynamics (default 0)
#' @param i_dur       Mean infectious period (default 5)
#' @param R0          Basic reproduction number (default 2.5)
#' @param profiletype One of "stepwise", "smooth", or "spike"
#'
#' @return Numeric vector of length n. Entry i is the infection time of
#'   individual i, or Inf if individual i was never infected.
sim_stochastic <- function(n=1000, e_dur=0, i_dur=5, R0=2.5, profiletype="stepwise"){

	beta <- R0/i_dur
	tinf <- rep(Inf,n)
	t <- 0

	# Seed infection
	newinf <- ceiling(runif(1)*n)
	tinf[newinf] <- t

	# Generate the seed's infection attempts
	queue <- c()
	if(profiletype=="stepwise"){
		t1 <- rexp(1,1/e_dur)
		t2 <- t1+rexp(1, 1/i_dur)
		nattempts <- rpois(1,(t2-t1)*beta)
		inftimes <- t + sort(runif(nattempts, min=t1, max=t2))
	} else if(profiletype=="smooth") {
		nattempts <- rpois(1,R0)
		inftimes <- t + sort(rexp(nattempts,1/e_dur) + rexp(nattempts,1/i_dur))
	} else if(profiletype=="spike"){
		nattempts <- rpois(1,R0)
		inftimes <- t + rep(rexp(1,1/e_dur) + rexp(1,1/i_dur), nattempts)
	} else {
		stop("Unrecognized profiletype")
	}
	queue <- sort(c(queue, inftimes))

	# Process events in chronological order
	while(length(queue)>0){

		t <- queue[1]
		queue <- queue[-1]

		# Draw a uniformly random target
		newinf <- ceiling(runif(1)*n)

		# If the target is susceptible, infect and generate their infection attempts
		if(tinf[newinf]==Inf){

			tinf[newinf] <- t
			if(profiletype=="stepwise"){
				t1 <- rexp(1,1/e_dur)
				t2 <- t1+rexp(1, 1/i_dur)
				nattempts <- rpois(1,(t2-t1)*beta)
				inftimes <- t + sort(runif(nattempts, min=t1, max=t2))
			} else if(profiletype=="smooth") {
				nattempts <- rpois(1,R0)
				inftimes <- t + sort(rexp(nattempts,1/e_dur) + rexp(nattempts,1/i_dur))
			} else if(profiletype=="spike"){
				nattempts <- rpois(1,R0)
				inftimes <- t + rep(rexp(1,1/e_dur) + rexp(1,1/i_dur), nattempts)
			} else {
				stop("Unrecognized profiletype")
			}
			queue <- sort(c(queue, inftimes))

		}

	}

	return(tinf)

}

# ==============================================================================
# Profile factory functions
# ==============================================================================
# Each factory takes epidemiological parameters and returns a closure with
# signature function(t_inf) -> numeric vector of infection attempt times.

#' Factory: stepwise infectiousness profile
#'
#' Infectious at constant rate beta during a random window
#' [onset, onset + duration], where onset ~ Exp(1/e_dur) and
#' duration ~ Exp(1/i_dur).
#'
#' @param e_dur Mean latent period
#' @param i_dur Mean infectious period
#' @param R0   Basic reproduction number
#' @return A function(t_inf) returning sorted infection attempt times
make_profile_stepwise <- function(e_dur, i_dur, R0) {
	beta <- R0 / i_dur
	function(t_inf) {
		t1 <- rexp(1, 1/e_dur)
		t2 <- t1 + rexp(1, 1/i_dur)
		nattempts <- rpois(1, (t2 - t1) * beta)
		if (nattempts == 0L) return(numeric(0))
		t_inf + sort(runif(nattempts, min = t1, max = t2))
	}
}

#' Factory: smooth infectiousness profile
#'
#' Poisson(R0) infection attempts, each independently timed from
#' Exp(1/e_dur) + Exp(1/i_dur).
#'
#' @inheritParams make_profile_stepwise
#' @return A function(t_inf) returning sorted infection attempt times
make_profile_smooth <- function(e_dur, i_dur, R0) {
	function(t_inf) {
		nattempts <- rpois(1, R0)
		if (nattempts == 0L) return(numeric(0))
		t_inf + sort(rexp(nattempts, 1/e_dur) + rexp(nattempts, 1/i_dur))
	}
}

#' Factory: spike infectiousness profile
#'
#' Poisson(R0) infection attempts all at one random time drawn from
#' Exp(1/e_dur) + Exp(1/i_dur).
#'
#' @inheritParams make_profile_stepwise
#' @return A function(t_inf) returning sorted infection attempt times
make_profile_spike <- function(e_dur, i_dur, R0) {
	function(t_inf) {
		nattempts <- rpois(1, R0)
		if (nattempts == 0L) return(numeric(0))
		t_inf + rep(rexp(1, 1/e_dur) + rexp(1, 1/i_dur), nattempts)
	}
}

#' Factory: mixture infectiousness profile (smooth <-> spike interpolation)
#'
#' Each individual draws Poisson(R0) infection attempts. With probability w, all infection attempts
#' share a single random time (spike); with probability 1-w, each infection attempt gets
#' an independent time (smooth). This interpolates between the smooth (w=0) and
#' spike (w=1) profiles.
#'
#' @inheritParams make_profile_stepwise
#' @param w Mixture weight in [0,1]. w=0 gives smooth, w=1 gives spike.
#' @return A function(t_inf) returning sorted infection attempt times
make_profile_mixture <- function(e_dur, i_dur, R0, w) {
	stopifnot(w >= 0, w <= 1)
	function(t_inf) {
		nattempts <- rpois(1, R0)
		if (nattempts == 0L) return(numeric(0))
		if (runif(1) < w) {
			# Spike: all infection attempts at one random time
			t_inf + rep(rexp(1, 1/e_dur) + rexp(1, 1/i_dur), nattempts)
		} else {
			# Smooth: each infection attempt gets independent timing
			t_inf + sort(rexp(nattempts, 1/e_dur) + rexp(nattempts, 1/i_dur))
		}
	}
}

#' Factory: shifted Gamma infectiousness profile
#'
#' Decomposes each infection attempt time into shift + jitter, exploiting the fact that
#' Gamma(alpha, r) + Gamma(kappa, r) = Gamma(alpha + kappa, r). Every
#' individual has the EXACT same profile shape Gamma(kappa, r), just shifted
#' to a different onset time s_i ~ Gamma(alpha_total - kappa, r).
#'
#' The population kernel A(tau) = R0 * Gamma(tau; alpha_total, r) is invariant
#' across kappa — changing punctuation does not change the population-level
#' generation interval.
#'
#' Properties:
#'   - All individual profiles have identical height, width, and shape
#'   - Each a_i integrates to R0 exactly
#'   - E[a_i(tau)] = A(tau) exactly, for all kappa
#'   - A(tau) is the same for all kappa (invariant population kernel)
#'
#' Limits:
#'   kappa -> 0:            a_i -> R0 * delta(tau - s_i)  (spike)
#'   kappa -> alpha_total:  a_i -> A(tau) for all i        (smooth)
#'
#' @param mu          Mean generation time (default e_dur + i_dur = 5)
#' @param R0          Basic reproduction number (default 2)
#' @param alpha_total Shape of the population kernel Gamma (default 10)
#' @param kappa       Profile shape in (0, alpha_total). Small = punctuated,
#'                    large = smooth. Must be > 0 and < alpha_total.
#' @return A function(t_inf) returning sorted infection attempt times
make_profile_gamma <- function(mu = 5, R0 = 2, alpha_total = 10, kappa) {
	stopifnot(kappa > 0, kappa < alpha_total)
	r <- alpha_total / mu
	alpha_shift <- alpha_total - kappa
	function(t_inf) {
		nattempts <- rpois(1, R0)
		if (nattempts == 0L) return(numeric(0))
		s_i <- rgamma(1, shape = alpha_shift, rate = r)
		t_inf + s_i + sort(rgamma(nattempts, shape = kappa, rate = r))
	}
}

#' Factory: shifted Gamma profile with time-varying contacts
#'
#' Extends make_profile_gamma by multiplying the biological timing density
#' b_i(tau) = f_kappa(tau - s_i) by a calendar-time contact function z(t).
#' The full individual infectiousness is:
#'   a_i(tau) = R0 * b_i(tau) * z(t_i + tau)
#'
#' Implementation uses Poisson thinning:
#' 1. Draw n_proposal ~ Poisson(R0 * z_max) candidate infection attempts
#' 2. Draw biological timing for each (onset shift s_i shared, jitter per infection attempt)
#' 3. Accept infection attempt j with probability z(t_inf + s_i + epsilon_j) / z_max
#'
#' The expected number of accepted infection attempts for individual i is:
#'   nu_i = R0 * integral[ b_i(tau) * z(t_i + tau) dtau ]
#' which varies across individuals — this is the superspreading mechanism.
#'
#' @param mu          Mean generation time (default 5)
#' @param R0          Basic reproduction number (default 2)
#' @param alpha_total Shape of the population kernel Gamma (default 10)
#' @param kappa       Profile shape in (0, alpha_total). Small = punctuated,
#'                    large = smooth.
#' @param contact_fn  Function z(t) giving contact rate multiplier at calendar
#'                    time t. Should average to ~1 over time.
#' @param z_max       Supremum of contact_fn, needed for thinning.
#' @return A function(t_inf) returning sorted infection attempt times
make_profile_gamma_contacts <- function(mu = 5, R0 = 2, alpha_total = 10, kappa,
                                        contact_fn, z_max) {
	stopifnot(kappa > 0, kappa < alpha_total)
	stopifnot(z_max > 0)
	r <- alpha_total / mu
	alpha_shift <- alpha_total - kappa
	function(t_inf) {
		# Oversample by factor z_max
		n_proposal <- rpois(1, R0 * z_max)
		if (n_proposal == 0L) return(numeric(0))
		# Onset shift shared across all infection attempts for this individual
		s_i <- rgamma(1, shape = alpha_shift, rate = r)
		# Independent jitter for each candidate infection attempt
		eps <- rgamma(n_proposal, shape = kappa, rate = r)
		# Calendar times of candidate infection attempts
		t_contacts <- t_inf + s_i + eps
		# Thinning: accept with probability z(t) / z_max
		accept_prob <- contact_fn(t_contacts) / z_max
		keep <- runif(n_proposal) < accept_prob
		if (!any(keep)) return(numeric(0))
		sort.int(t_contacts[keep])
	}
}

#' Faster implementation of sim_stochastic (identical model logic)
#'
#' Speed-ups over the reference implementation:
#' - Profile-specific infection attempt generator resolved once via switch(), avoiding
#'   repeated if/else branching on every event.
#' - Index-based queue consumption: failed infection attempts (target already infected)
#'   advance a pointer in O(1) rather than copying the queue via queue[-1].
#' - sort.int() instead of sort() to skip S3 method dispatch.
#' - Early return when rpois draws zero infection attempts.
#'
#' @inheritParams sim_stochastic
#' @param gen_inf_attempts Optional infection attempt generator function with signature
#'   function(t_inf) -> numeric vector. If provided, overrides profiletype.
#' @return Numeric vector of length n (same format as sim_stochastic).
sim_stochastic_fast <- function(n=1000, e_dur=0, i_dur=5, R0=2.5,
                                profiletype="stepwise", gen_inf_attempts=NULL){

	beta <- R0/i_dur
	tinf <- rep(Inf, n)

	# Define the infection attempt generation function once at entry
	if (is.null(gen_inf_attempts)) {
		gen_inf_attempts <- switch(profiletype,
			"stepwise" = function(t_inf) {
				t1 <- rexp(1, 1/e_dur)
				t2 <- t1 + rexp(1, 1/i_dur)
				nattempts <- rpois(1, (t2 - t1) * beta)
				if(nattempts == 0L) return(numeric(0))
				t_inf + sort(runif(nattempts, min = t1, max = t2))
			},
			"smooth" = function(t_inf) {
				nattempts <- rpois(1, R0)
				if(nattempts == 0L) return(numeric(0))
				t_inf + sort(rexp(nattempts, 1/e_dur) + rexp(nattempts, 1/i_dur))
			},
			"spike" = function(t_inf) {
				nattempts <- rpois(1, R0)
				if(nattempts == 0L) return(numeric(0))
				t_inf + rep(rexp(1, 1/e_dur) + rexp(1, 1/i_dur), nattempts)
			},
			stop("Unrecognized profiletype")
		)
	}

	# Seed infection
	seed <- sample.int(n, 1)
	tinf[seed] <- 0
	queue <- gen_inf_attempts(0)
	qi <- 1L
	qn <- length(queue)

	# Process events in chronological order
	while(qi <= qn) {
		t_event <- queue[qi]
		qi <- qi + 1L

		target <- sample.int(n, 1)
		if(tinf[target] == Inf) {
			tinf[target] <- t_event
			new_events <- gen_inf_attempts(t_event)
			if(length(new_events) > 0L) {
				# Merge new events with unprocessed remainder of queue
				remaining <- if(qi <= qn) queue[qi:qn] else numeric(0)
				queue <- sort.int(c(remaining, new_events))
				qi <- 1L
				qn <- length(queue)
			}
		}
	}

	return(tinf)

}
