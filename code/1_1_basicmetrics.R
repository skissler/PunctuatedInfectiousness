# Simulate 1000 epidemics in populations of size 1000 
# Grab: time to pass 10 infected, 50 infected, 100 infected; peak time; exponential growth rate. Extinction probability. Final size. 

library(tidyverse) 
library(odin) 

source("code/0_utils.R") 

# ==============================================================================
# Stochastic epidemic simulations: uncontrolled dynamics 
# 
# Corresponds to writeup sections: 
#   - "The impact of individual infectiousness profile on uncontrolled epidemic
#      dynamics"
# 
# Runs nsim stochastic epidemics for each individual infectiousness profile 
# (stepwise, smooth, spike) and: 
#   - Plots epidemic trajectories (daily incidence and cumulative incidence)
#   - Calculates extinction probabilities and final size distributions 
#   - Plots the time to reach various epidemic milestones (100 infections,
#     peak) as Kaplan-Meier curves. 
#   - Calculates empirical epidemic growth rates 
# ==============================================================================

# ==============================================================================
# 1. Parameters
# ==============================================================================

popsize <- 1000
e_dur   <- 2 
i_dur   <- 3 
R0      <- 5
nsim    <- 1000 

# ==============================================================================
# 2. Mean-field ODE solution
# ==============================================================================

# Run the SEIR model 
ode_model <- seir$new(R0=R0, nu=1/e_dur, gamma=1/i_dur, N=popsize)
t <- seq(from=0, to=100, length.out=1000)
ode_out <- as_tibble(data.frame(ode_model$run(t)))

# Aggregate ODE output to daily and weekly new-infection counts 
ode_daily <- ode_out %>% 
	mutate(day=floor(t)) %>% 
	group_by(day) %>% 
	summarise(cuminf=max(cuminf)) %>% 
	arrange(day) %>% 
	mutate(newinf=cuminf-lag(cuminf)) %>% 
	mutate(newinf=case_when(is.na(newinf)~cuminf, TRUE~newinf))

ode_weekly <- ode_out %>% 
	mutate(day=floor(t), week=floor(day/7)) %>% 
	group_by(week) %>% 
	summarise(cuminf=max(cuminf)) %>% 
	arrange(week) %>% 
	mutate(newinf=cuminf-lag(cuminf)) %>% 
	mutate(newinf=case_when(is.na(newinf)~cuminf, TRUE~newinf))

# ==============================================================================
# 3. Stochastic simulations 
# ==============================================================================

# Define profile types to loop over
profiles <- list(
	stepwise = gen_inf_attempts_stepwise(e_dur, i_dur, R0),
	smooth   = gen_inf_attempts_smooth(e_dur, i_dur, R0),
	spike    = gen_inf_attempts_spike(e_dur, i_dur, R0)
)

# Run the simulations and store the output
cuminf_df <- tibble()
for(sim in 1:nsim){
    for(pname in names(profiles)){
            tinf <- sim_stochastic_fast(n=popsize, 
            	                        gen_inf_attempts=profiles[[pname]])
            cuminf_df <- bind_rows(cuminf_df,
                    tibble(tinf=sort(tinf[tinf<Inf])) %>%
                            mutate(cuminf=1:n(), sim=sim, profiletype=pname))
    }
    if(sim %% 100 == 0) print(sim)
}

# Flag established epidemics (>= 10% of population infected)
cuminf_df <- cuminf_df %>%
	mutate(profiletype = factor(profiletype, levels = names(profiles))) %>%
	group_by(sim, profiletype) %>%
	mutate(established = as.integer(n() >= 0.1 * popsize)) %>%
	ungroup()

# Aggregate to daily/weekly resolution
dayjoin <- expand_grid(
	profiletype=unique(cuminf_df$profiletype),
	sim=1:nsim,
	day=0:max(ceiling(cuminf_df$tinf)))

weekjoin <- expand_grid(
	profiletype=unique(cuminf_df$profiletype),
	sim=1:nsim,
	week=0:max(ceiling(cuminf_df$tinf / 7)))

dailyinf_df <- cuminf_df %>%
	mutate(day = floor(tinf)) %>%
	group_by(profiletype, sim, day) %>%
	summarise(ninf = n(), established = first(established), .groups = "drop") %>%
	right_join(dayjoin, by = c("profiletype", "sim", "day")) %>%
	replace_na(list(ninf = 0)) %>%
	group_by(profiletype, sim) %>%
	arrange(day, .by_group = TRUE) %>%
	mutate(established = max(established, na.rm = TRUE)) %>%
	mutate(cuminf=cumsum(ninf))

weeklyinf_df <- cuminf_df %>%
	mutate(day = floor(tinf), week = floor(day / 7)) %>%
	group_by(profiletype, sim, week) %>%
	summarise(ninf = n(), established = first(established), .groups = "drop") %>%
	right_join(weekjoin, by = c("profiletype", "sim", "week")) %>%
	replace_na(list(ninf = 0)) %>%
	group_by(profiletype, sim) %>%
	arrange(week, .by_group = TRUE) %>%
	mutate(established = max(established, na.rm = TRUE)) %>%
	mutate(cuminf=cumsum(ninf))

lastday <- ceiling(max(cuminf_df$tinf))

dailyinf_means <- dailyinf_df %>%
	filter(established == 1) %>%
	group_by(profiletype, day) %>%
	summarise(mean_ninf = mean(ninf), mean_cuminf = mean(cuminf), .groups = "drop")

weeklyinf_means <- weeklyinf_df %>%
	filter(established == 1) %>%
	group_by(profiletype, week) %>%
	summarise(mean_ninf = mean(ninf), mean_cuminf = mean(cuminf), .groups = "drop")

# ==============================================================================
# 4. Key metrics
# ==============================================================================

pest_table <- cuminf_df %>%
	group_by(sim, profiletype) %>%
	slice(1) %>%
	group_by(profiletype) %>%
	summarise(pestablished=sum(established)/nsim, .groups="drop")

fs_table <- cuminf_df %>% 
	group_by(sim, profiletype) %>% 
	filter(established==1) %>% 
	summarise(finalsize=max(cuminf), .groups="drop") %>% 
	group_by(profiletype) %>% 
	summarise(fs_mean=mean(finalsize), fs_sd=sd(finalsize), .groups="drop")

# ==============================================================================
# 5. Figures — epidemic trajectories
# ==============================================================================

# Cumulative stochastic epidemic curves (grey) with ODE overlay (blue) 
fig_cuminf_overlay <- cuminf_df %>%
	ggplot(aes(x=tinf, y=cuminf, group=sim)) +
		geom_line(alpha=0.2, col="grey") +
		geom_line(data=filter(ode_daily, day<=lastday),
			aes(x=day, y=cuminf*popsize),
			alpha=0.6, linewidth=1, col="blue") +
		# geom_line(data=dailyinf_means, aes(x=day, y=mean_cuminf), col="black") +
		theme_classic() +
		facet_wrap(~profiletype, nrow=1)

# Daily stochastic incidence curves (grey) with ODE overlay (blue) 
fig_inf_overlay <- dailyinf_df %>% 
	ggplot(aes(x=day, y=ninf, group=sim)) + 
		geom_line(alpha=0.2, col="grey") + 
		geom_line(data=filter(ode_daily, day<=lastday), 
			aes(x=day, y=newinf*popsize),
			alpha=0.6, linewidth=1, col="blue") + 
		# geom_line(data=dailyinf_means, aes(x=day, y=mean_ninf), col="black") +
		theme_classic() + 
		facet_wrap(~profiletype, nrow=1)

# ==============================================================================
# 6. Figure — "survival curve" for time to reach 100 cases
# ==============================================================================

# For each established epidemic, extract the time it first reaches 100 infections
threshold <- 100

time_to_threshold <- cuminf_df %>%
	filter(established == 1) %>%
	group_by(sim, profiletype) %>%
	filter(cuminf >= threshold) %>%
	slice(1) %>%
	ungroup() %>%
	select(sim, profiletype, tinf)

# Build the survival curve: for a grid of times, compute the fraction
# of epidemics that haven't yet reached the threshold
t_grid <- seq(0, max(time_to_threshold$tinf), length.out = 500)

survival_df <- time_to_threshold %>%
	group_by(profiletype) %>%
	summarise(
		t = list(t_grid),
		prop_below = list(sapply(t_grid, function(tt) mean(tinf > tt))),
		.groups = "drop") %>%
	unnest(cols = c(t, prop_below))

breakvals = if(max(cuminf_df$tinf)>120){seq(0,365,by=15)}else{seq(0,365,by=7)}
fig_survival_100 <- ggplot(survival_df, aes(x=t, y=prop_below, col=profiletype)) +
	geom_line(alpha=0.6, linewidth=0.8) +
	scale_x_continuous(breaks=breakvals) + 
	scale_color_manual(values=c("stepwise"="black", "smooth"="blue", "spike"="red")) + 
	labs(x = "Days",
	     y = paste0("Proportion not yet reaching ", threshold, " cases"),
	     col = "Profile") +
	theme_classic()

# ==============================================================================
# 7. Figure — "survival curve" for time to reach peak daily incidence
# ==============================================================================

# For each established epidemic, find the day of peak daily incidence
time_to_peak <- dailyinf_df %>%
	filter(established == 1) %>%
	group_by(sim, profiletype) %>%
	slice_max(ninf, n=1, with_ties=FALSE) %>%
	ungroup() %>%
	select(sim, profiletype, day)

# Build survival curve: fraction of epidemics that haven't yet peaked
t_grid_peak <- seq(0, max(time_to_peak$day), length.out = 500)

survival_peak_df <- time_to_peak %>%
	group_by(profiletype) %>%
	summarise(
		t = list(t_grid_peak),
		prop_below = list(sapply(t_grid_peak, function(tt) mean(day > tt))),
		.groups = "drop") %>%
	unnest(cols = c(t, prop_below))

breakvals = if(max(cuminf_df$tinf)>120){seq(0,365,by=15)}else{seq(0,365,by=7)}
fig_survival_peak <- ggplot(survival_peak_df, aes(x=t, y=prop_below, col=profiletype)) +
	geom_line(alpha=0.6, linewidth=0.8) +
	scale_x_continuous(breaks=breakvals) +
	scale_color_manual(values=c("stepwise"="black", "smooth"="blue", "spike"="red")) +
	labs(x = "Days",
	     y = "Proportion not yet reaching peak daily incidence",
	     col = "Profile") +
	theme_classic()

