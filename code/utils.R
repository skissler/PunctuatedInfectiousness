library(tidyverse) 
library(odin)

# define the mean-field ODE: 
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


# Define a stochastic simulation function
sim_stochastic <- function(n=1000, e_dur=0, i_dur=5, R0=2.5, profiletype="stepwise"){

	# Transform inputs: 
	beta <- R0/i_dur

	# Initialize infection times
	tinf <- rep(Inf,n)

	# Initialize the time
	t <- 0 

	# Set a seed infection 
	newinf <- ceiling(runif(1)*n)
	tinf[newinf] <- t

	# Initialize the event queue: 
	queue <- c()
	if(profiletype=="stepwise"){
		t1 <- rexp(1,1/e_dur)
		t2 <- t1+rexp(1, 1/i_dur)
		nattempts <- rpois(1,(t2-t1)*beta) # the infectiousness profile assumption
		# nattempts <- rpois(1,R0) # the all-equal assumption
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

	# As long as there are more times in the queue 
	while(length(queue)>0){

		# Update the current time 
		t <- queue[1]
		queue <- queue[-1]

		# Draw a potential infectee: 
		newinf <- ceiling(runif(1)*n)

		# If the potential infectee is still susceptible: 
		if(tinf[newinf]==Inf){

			tinf[newinf] <- t
			if(profiletype=="stepwise"){
				t1 <- rexp(1,1/e_dur)
				t2 <- t1+rexp(1, 1/i_dur)
				nattempts <- rpois(1,(t2-t1)*beta) # the infectiousness profile assumption
				# nattempts <- rpois(1,R0) # the all-equal assumption
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

# Faster implementation of sim_stochastic with identical logic:
# - Contact generation defined once via switch (avoids repeated if/else branching)
# - Index-based queue consumption (avoids O(queue_size) vector copy on failed contacts)
# - sort.int instead of sort (less dispatch overhead)
# - Early return on zero contacts
sim_stochastic_fast <- function(n=1000, e_dur=0, i_dur=5, R0=2.5, profiletype="stepwise"){

	beta <- R0/i_dur
	tinf <- rep(Inf, n)

	# Define contact generation function once
	gen_contacts <- switch(profiletype,
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

	# Seed infection
	seed <- sample.int(n, 1)
	tinf[seed] <- 0
	queue <- gen_contacts(0)
	qi <- 1L
	qn <- length(queue)

	while(qi <= qn) {
		t_event <- queue[qi]
		qi <- qi + 1L

		target <- sample.int(n, 1)
		if(tinf[target] == Inf) {
			tinf[target] <- t_event
			new_events <- gen_contacts(t_event)
			if(length(new_events) > 0L) {
				remaining <- if(qi <= qn) queue[qi:qn] else numeric(0)
				queue <- sort.int(c(remaining, new_events))
				qi <- 1L
				qn <- length(queue)
			}
		}
	}

	return(tinf)

}

