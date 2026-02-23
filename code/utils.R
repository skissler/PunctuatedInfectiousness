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




# As inputs: 
# tinf (infection times for everyone) 
# params (profile parameters for each person: onset/offset, delta time, or erlang shape) 
# t (actual time)

# resample <- function(x, ...) x[sample.int(length(x), ...)]

# # getnext_stepwise(t=2.5, beta=2, tinf=tinf, t1vec=t1vec, t2vec=t2vec)
# getnext_stepwise <- function(t, beta, tinf, t1vec, t2vec){

# 	if(sum(tinf<Inf & tinf>t)>0){
# 		warning("At least one person has an infection time in the future.\nAll infection times should be in the past (<t) or Inf")
# 	}

# 	# Calculate remaining infectiousness
# 	infremaining <- pmax(pmin(beta*(tinf+t2vec-t), beta*(t2vec-t1vec)),0)
# 	infremaining[tinf==Inf] <- 0

# 	# Draw number of future attempted infections 
# 	nattempts <- unlist(lapply(infremaining, function(x){rpois(1,x)}))
# 	# print(nattempts)

# 	# Draw timing of those attempted infections 
# 	# Report the earliest time of attempted infection and number of infections at that time
# 	mintimes <- unlist(lapply(1:length(tinf), function(i){
# 		if(nattempts[i]>0){
# 			return(min(runif(nattempts[i]))*(tinf[i] + t2vec[i] - max(t, tinf[i]+t1vec[i])) + max(t, tinf[i]+t1vec[i]))
# 		} else {
# 			return(Inf)
# 		}
# 		}))

# 	# sort(runif(nattempts[i]))*(tinf[i] + t2vec[i] - max(t, tinf[i]+t1vec[i])) + max(t, tinf[i]+t1vec[i])

# 	return(
# 		c(
# 			time=min(mintimes), 
# 			n_attempts=1, 
# 			infector=resample(which(mintimes==min(mintimes)),size=1))
# 		)

# }


# getnext_spike <- function(t, beta, tinf, tstar){
# 	if(sum(tinf<Inf & tinf>t)>0){
# 			warning("At least one person has an infection time in the future.\nAll infection times should be in the past (<t) or Inf")
# 		}



# }

# getnext_smooth <- function(){}



# sim_stochastic_rightorder_2 <- function(n=1000, e_dur=0, i_dur=5, R0=2.5, profiletype="stepwise"){

# 	# Transform inputs: 
# 	beta <- R0/i_dur

# 	# Define the attempt times: 
# 	t1vec <- rexp(n, 1/e_dur)
# 	t2vec <- rexp(n, 1/i_dur)

# 	rates <- beta*(t2vec - t1vec)

# 	nattempts <- unlist(lapply(rates, function(x){rpois(1,x)}))
# 	tattempts <- lapply(1:n, function(i){
# 		sort(runif(nattempts[i]))*(t2vec[i]-t1vec[i]) + t1vec[i]
# 		})

# 	# Initialize infection times
# 	tinf <- rep(Inf,n)
# 	# Initialize infector array
# 	infectors <- rep(NA,n)

# 	# Initialize the time
# 	t <- 0 

# 	# Set a seed infection 
# 	newinf <- ceiling(runif(1)*n)
# 	tinf[newinf] <- t

# }



# # tinf <- sim_stochastic_rightorder(n=1000, e_dur=0, i_dur=5, R0=2)
# sim_stochastic_rightorder <- function(n=1000, e_dur=0, i_dur=5, R0=2.5, profiletype="stepwise"){

# 	# Transform inputs: 
# 	beta <- R0/i_dur
# 	if(profiletype=="stepwise"){
# 		t1vec <- rexp(n, 1/e_dur)
# 		t2vec <- t1vec + rexp(n, 1/i_dur)
# 	}

# 	# Initialize infection times
# 	tinf <- rep(Inf,n)
# 	# Initialize infector array
# 	infectors <- rep(NA,n)

# 	# Initialize the time
# 	t <- 0 

# 	# Set a seed infection 
# 	newinf <- ceiling(runif(1)*n)
# 	tinf[newinf] <- t
# 	infectors[newinf] <- 0

# 	# While there are still susceptible people to infect: 
# 	while(max(tinf) == Inf){

# 		# When does the next attempted infection happen? 
# 		nextevent <- getnext_stepwise(t, beta, tinf, t1vec, t2vec)

# 		if(nextevent["time"] == Inf){break()}

# 		# Who do we try to infect? 
# 		newattempts <- unique(sample(setdiff(1:n, nextevent["infector"]), 
# 			size=nextevent["n_attempts"],
# 			replace=TRUE))
# 		# If these people are susceptible, update their infection times: 
# 		for(attempt in newattempts){
# 			if(tinf[attempt] == Inf){
# 				tinf[attempt] <- nextevent["time"]
# 				infectors[attempt] <- nextevent["infector"]
# 			}
# 		}

# 		# Update the time
# 		t <- nextevent["time"]

# 	}

# 	return(tinf)

# }


# tinf <- sim_stochastic_rightorder(n=1000, e_dur=0, i_dur=5, R0=2)
# toplot <- tibble(tinf=sort(tinf[tinf<Inf])) %>% 
# 				mutate(cuminf=1:n())
# # ggplot(toplot, aes(x=tinf, y=cuminf)) + 
# # 	geom_line() 
# toplot %>% 
# 	ggplot(aes(x=tinf, y=cuminf)) + 
# 		geom_line(alpha=0.2, col="grey") + 
# 		geom_line(data=sir_out_daily, aes(x=day, y=cuminf*1000), alpha=0.6, linewidth=1, col="blue") + 
# 		theme_classic() 



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


# Define a stochastic simulation function, with explicit contact rates and probability-of-infection profiles: 
sim_stochastic <- function(n=1000, e_dur=0, i_dur=5, R0=2.5, contactfun=function(x){return(25)}, profiletype="stepwise"){

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



sim_stochastic_fast <- function(n = 1000, e_dur = 0, i_dur = 5, R0 = 2.5, 
                           contactfun = function(x) 25, profiletype = "stepwise") {
  
  # Precompute beta
  beta <- R0 / i_dur
  
  # Function to generate infectious contact times for a single individual
  gen_infectious_times <- function(t_current) {
    if (profiletype == "stepwise") {
      t1 <- rexp(1, 1 / e_dur)
      t2 <- t1 + rexp(1, 1 / i_dur)
      nattempts <- rpois(1, (t2 - t1) * beta)
      return(t_current + sort(runif(nattempts, min = t1, max = t2)))
    } else if (profiletype == "smooth") {
      nattempts <- rpois(1, R0)
      return(t_current + sort(rexp(nattempts, 1 / e_dur) + rexp(nattempts, 1 / i_dur)))
    } else if (profiletype == "spike") {
      nattempts <- rpois(1, R0)
      spike_time <- rexp(1, 1 / e_dur) + rexp(1, 1 / i_dur)
      return(rep(t_current + spike_time, nattempts))
    } else {
      stop("Unrecognized profiletype")
    }
  }
  
  # Initialize infection times
  tinf <- rep(Inf, n)
  
  # Seed infection
  first_inf <- sample.int(n, 1)
  tinf[first_inf] <- 0
  queue <- gen_infectious_times(0)
  
  # Main loop
  while (length(queue) > 0) {
    # Take next event time
    t <- queue[1]
    queue <- queue[-1]
    
    # Draw potential infectee
    candidate <- sample.int(n, 1)
    
    if (is.infinite(tinf[candidate])) {
      tinf[candidate] <- t
      new_events <- gen_infectious_times(t)
      # Efficient insertion maintaining sort order
      queue <- sort(c(queue, new_events), method = "quick")
    }
  }
  
  return(tinf)
}


sim_contact <- function(mean_contacts=10, rate_out=5){

	rate_in <- mean_contacts*rate_out
	ncontacts <- mean_contacts
	cvec <- c(ncontacts)
	tvec <- c(0)

	for(indexA in 1:10000){
		cumrate <- rate_in + ncontacts*rate_out
		nexttime <- rexp(1, cumrate);
		nextevent <- if(runif(1) < (rate_in/cumrate)){"in"} else {"out"}
		if(nextevent=="in"){ncontacts <- ncontacts+1} else {ncontacts <- ncontacts-1}
		cvec <- c(cvec, ncontacts)
		tvec <- c(tvec, nexttime)
	}
	
	return(tibble(t=cumsum(tvec), c=cvec))

}


meanfun1 <- function(t){
	out <- -4*cos(2*pi*t) + 8
	return(out)
}

tvals <- seq(from=0, to=7, by=0.01)
meancvals <- meanfun1(tvals)
plot(tvals, meancvals, type='l')

sim_contact_dynamic <- function(meanfun, rate_out){

	


}

# ggplot(temp, aes(x=t, y=c)) + geom_line() + theme_classic() 




