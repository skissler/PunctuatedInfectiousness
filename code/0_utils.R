library(tidyverse) 
library(odin) 

# ==============================================================================
# Utility functions
# 
# Defines functions for simulating deterministic and stochastic epidemics, with
# helper functions 
# ==============================================================================

# ==============================================================================
# Mean-field ODE model
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
#'   infectiousness profile (kappa \in (0, popshape), with small kappa giving
#'   sharp infectiousness profiles')
#' @return A function(t_inf) that returns sorted infection attempt times
gen_inf_attempts_gamma <- function(T, R0, popshape, kappa){
	stopifnot(kappa>0, kappa<popshape)
	r <- popshape/T
	shiftshape <- popshape - kappa
	function(tinf){
		nattempts <- rpois(1,R0)
		if(nattempts==0L) return(numeric(0))
		shift <- rgamma(1, shape=shiftshape, rate=r)
		attempt_times <- tinf + shift + sort(rgamma(nattempts, shape=kappa, rate=r))
		return(attempt_times)
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

