library(tidyverse) 
library(odin) 

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

gen_inf_attempts_stepwise <- function(tinf, e_dur, i_dur, R0){
	beta <- R0/i_dur
	t1 <- rexp(1, 1/e_dur)
	t2 <- t1 + rexp(1, 1/i_dur)
	nattempts <- rpois(1, (t2 - t1)*beta) 
	if(nattempts == 0L) return(numeric(0))
	attempt_times <- tinf + sort(runif(nattempts, min=t1, max=t2))
	return(attempt_times)
}

gen_inf_attempts_smooth <- function(tinf, e_dur, i_dur, R0){
	nattempts <- rpois(1, R0)
	if(nattempts == 0L) return(numeric(0))
	attempt_times <- tinf + sort(rexp(nattempts, 1/e_dur) + rexp(nattempts, 1/i_dur))
	return(attempt_times)
}

gen_inf_attempts_spike <- function(tinf, e_dur, i_dur, R0){
	nattempts <- rpois(1, R0)
	if(nattempts == 0L) return(numeric(0))
	attempt_times <- tinf + rep(rexp(1, 1/e_dur) + rexp(1, 1/i_dur), nattempts)
	return(attempt_times)
}


sim_stochastic_fast <- function(n=1000, e_dur=2, i_dur=3, R0=2, 
	                            gen_inf_attempts){

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

