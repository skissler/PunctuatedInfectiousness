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
#' @param R0        Basic reproduction number
#' @param popshape  Shape parameter of Gamma generation interval distribution
#' @param T_gen     Mean generation interval (days)
#' @param N         Population size, used to set initial conditions
#' @param dt    '   Time step for discretization (default 0.01 days)
#' @param tmax      Maximum simulation length (default 100 days)

renewal_epidemic <- function(R0, popshape, T_gen, N, dt = 0.01, tmax = 100) {

	# Set rate parameter of the generation interval distribution
  r <- popshape / T_gen

  # Initalize time vectors 
  times <- seq(0, tmax, by = dt)
  nt <- length(times)
  
  # Initialize incidence and cumulative infections 
  j <- numeric(nt)
  C <- numeric(nt)
  C[1] <- 1 / N

  # Generation interval weights (with dt baked in for the convolution)
  g <- dgamma(times, shape = popshape, rate = r) * dt
  # Index case kernel (no dt — it's a single person, not a rate)
  g_index <- dgamma(times, shape = popshape, rate = r) / N

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
#' with punctuation parameter kappa. 
#' 
#' @param T Mean generation interval 
#' @param R0 Basic reproduction number 
#' @param popshape Shape parameter for the Gamma-distributed population 
#'   infectiousness profile (A(tau) ~ Gamma(popshape, popshape/T))
#' @param kappa Parameter governing punctuation of the individual
#'   infectiousness profile (kappa \in [0, 1], with kappa -> 0 giving sharp
#'   infectiousness profiles and kappa -> 1 giving smooth ones equivalent to 
#'   the population generation interval distribution)
#' @return A function(t_inf) that returns sorted infection attempt times
gen_inf_attempts_gamma <- function(T, R0, popshape, kappa){
	stopifnot(kappa>=0, kappa<=1)
	r <- popshape/T

	if(kappa < 1e-6){ # spike implementation
		shiftshape <- popshape 
		function(tinf){
			nattempts <- rpois(1,R0)
			if(nattempts==0L) return(numeric(0))
			shift <- rgamma(1, shape=shiftshape, rate=r)
			attempt_times <- rep(tinf + shift, nattempts)
			return(sort(attempt_times))
		}
	} else if(kappa > 1-1e-6){ # smooth implementation
		function(tinf){
			nattempts <- rpois(1,R0)
			if(nattempts==0L) return(numeric(0))
			attempt_times <- tinf + rgamma(nattempts, shape=popshape, rate=r)
			return(sort(attempt_times))
		}
	} else { # regular implementation
		shiftshape <- (1-kappa)*popshape
		function(tinf){
			nattempts <- rpois(1,R0)
			if(nattempts==0L) return(numeric(0))
			shift <- rgamma(1, shape=shiftshape, rate=r)
			attempt_times <- tinf + shift + rgamma(nattempts, shape=kappa*popshape, rate=r)
			return(sort(attempt_times))
		}
	}
	
}

gen_inf_attempts_gamma_contacts <- function(T, z, z_max, popshape, kappa){

	stopifnot(kappa>0, kappa<1)
	stopifnot(z_max>0) 

	r <- popshape/T 
	shiftshape <- (1-kappa)*popshape 

	function(tinf){
		nproposals <- rpois(1, z_max)
		if(nproposals == 0L) return(numeric(0))
		shift <- rgamma(1, shape=shiftshape, rate=r)
		proposals <- tinf + shift + rgamma(nproposals, shape=kappa*popshape, rate=r)
		accept_prob <- z(proposals)/z_max
		keep <- runif(nproposals) < accept_prob
		if(!any(keep)) return(numeric(0))
		return(sort(proposals[keep]))
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

