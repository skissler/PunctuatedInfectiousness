library(tidyverse) 
library(odin) 
source("code/0_utils.R") 

# ==============================================================================
# Stochastic epidemic simulations: uncontrolled dynamics 
# 
# Corresponds to writeup sections: 
#   - "The impact of individual infectiousness profile on uncontrolled 
#      epidemic dynamics"
# 
# Runs nsim stochastic epidemics for spike (psi=0), smooth (psi=1), and
# intermediate (psi=0.25) individual infectiousness profiles, and: 
#   - Plots epidemic trajectories (daily incidence and cumulative incidence)
#   - Calculates extinction probabilities and final size distributions 
#   - Plots the time to reach various epidemic milestones (100 infections,
#     peak) as Kaplan-Meier curves. 
#   - Calculates empirical epidemic growth rates 
# ==============================================================================

# ==============================================================================
# 1. Parameters
# ==============================================================================

popsize  <- 1000 
T        <- 7.05 # 12.2 # 3.2 # 5
popshape <- 2.39 # 11.36 # 2.32 # 10 
poprate  <- popshape/T
R0       <- 6 # 12 # 2 # 12
nsim     <- 1000 

# ==============================================================================
# 2. Mean-field solution
# ==============================================================================

# Run the renewal equation model 
ren_out <- renewal_epidemic(R0, popshape, T, 1000)

# Aggregate renewal equation output to daily and weekly new-infection counts 
ren_daily <- ren_out %>% 
	mutate(day=floor(t)) %>% 
	group_by(day) %>% 
	summarise(cuminf=max(cuminf)) %>% 
	arrange(day) %>% 
	mutate(newinf=cuminf-lag(cuminf)) %>% 
	mutate(newinf=case_when(is.na(newinf)~cuminf, TRUE~newinf))

ren_weekly <- ren_out %>% 
	mutate(day=floor(t), week=floor(day/7)) %>% 
	group_by(week) %>% 
	summarise(cuminf=max(cuminf)) %>% 
	arrange(week) %>% 
	mutate(newinf=cuminf-lag(cuminf)) %>% 
	mutate(newinf=case_when(is.na(newinf)~cuminf, TRUE~newinf))

# ==============================================================================
# 3. Stochastic simulations 
# ==============================================================================

# Define kappa values to loop over
kappavals <- c(0, 0.5, 1) 

# Run the simulations and store the output
cuminf_df <- tibble()
for(sim in 1:nsim){
    for(kappa in kappavals){
            tinf <- sim_stochastic_fast(n=popsize, 
            	                        gen_inf_attempts=gen_inf_attempts_gamma(T, R0, popshape, kappa))
            cuminf_df <- bind_rows(cuminf_df,
                    tibble(tinf=sort(tinf[tinf<Inf])) %>%
                            mutate(cuminf=1:n(), sim=sim, kappa=kappa))
    }
    if(sim %% 100 == 0) print(sim)
}

# Flag established epidemics (>= 10% of population infected)
cuminf_df <- cuminf_df %>%
	mutate(kappa = factor(kappa)) %>%
	group_by(sim, kappa) %>%
	mutate(established = as.integer(n() >= 0.1 * popsize)) %>%
	ungroup()


# Aggregate to daily/weekly resolution
dayjoin <- expand_grid(
	kappa=unique(cuminf_df$kappa),
	sim=1:nsim,
	day=0:max(ceiling(cuminf_df$tinf)))

weekjoin <- expand_grid(
	kappa=unique(cuminf_df$kappa),
	sim=1:nsim,
	week=0:max(ceiling(cuminf_df$tinf / 7)))

dailyinf_df <- cuminf_df %>%
	mutate(day = floor(tinf)) %>%
	group_by(kappa, sim, day) %>%
	summarise(ninf = n(), established = first(established), .groups = "drop") %>%
	right_join(dayjoin, by = c("kappa", "sim", "day")) %>%
	replace_na(list(ninf = 0)) %>%
	group_by(kappa, sim) %>%
	arrange(day, .by_group = TRUE) %>%
	mutate(established = max(established, na.rm = TRUE)) %>%
	mutate(cuminf=cumsum(ninf))

weeklyinf_df <- cuminf_df %>%
	mutate(day = floor(tinf), week = floor(day / 7)) %>%
	group_by(kappa, sim, week) %>%
	summarise(ninf = n(), established = first(established), .groups = "drop") %>%
	right_join(weekjoin, by = c("kappa", "sim", "week")) %>%
	replace_na(list(ninf = 0)) %>%
	group_by(kappa, sim) %>%
	arrange(week, .by_group = TRUE) %>%
	mutate(established = max(established, na.rm = TRUE)) %>%
	mutate(cuminf=cumsum(ninf))

lastday <- ceiling(max(cuminf_df$tinf))

dailyinf_means <- dailyinf_df %>%
	filter(established == 1) %>%
	group_by(kappa, day) %>%
	summarise(mean_ninf = mean(ninf), mean_cuminf = mean(cuminf), .groups = "drop")

weeklyinf_means <- weeklyinf_df %>%
	filter(established == 1) %>%
	group_by(kappa, week) %>%
	summarise(mean_ninf = mean(ninf), mean_cuminf = mean(cuminf), .groups = "drop")

# ==============================================================================
# 4. Key metrics
# ==============================================================================

pest_table <- cuminf_df %>%
	group_by(sim, kappa) %>%
	slice(1) %>%
	group_by(kappa) %>%
	summarise(pestablished=sum(established)/nsim, .groups="drop")

fs_table <- cuminf_df %>% 
	group_by(sim, kappa) %>% 
	filter(established==1) %>% 
	summarise(finalsize=max(cuminf), .groups="drop") %>% 
	group_by(kappa) %>% 
	summarise(fs_mean=mean(finalsize), fs_sd=sd(finalsize), .groups="drop")

# ==============================================================================
# 5. Figures — epidemic trajectories
# ==============================================================================

# Cumulative stochastic epidemic curves (grey) with ODE overlay (blue) 
fig_cuminf_overlay <- cuminf_df %>%
	ggplot(aes(x=tinf, y=cuminf, group=sim)) +
		geom_line(alpha=0.2, col="grey") +
		geom_line(data=filter(ren_out, t<=lastday),
			aes(x=t, y=cuminf*popsize),
			alpha=0.6, linewidth=1, col="blue") +
		theme_classic() +
		facet_wrap(~kappa, nrow=1)

# Daily stochastic incidence curves (grey) with ODE overlay (blue) 
fig_inf_overlay <- dailyinf_df %>% 
	ggplot(aes(x=day, y=ninf, group=sim)) + 
		geom_line(alpha=0.2, col="grey") + 
		geom_line(data=filter(ren_daily, day<=lastday), 
			aes(x=day, y=newinf*popsize),
			alpha=0.6, linewidth=1, col="blue") + 
		theme_classic() + 
		facet_wrap(~kappa, nrow=1)


# ==============================================================================
# 6. Figure — KM curves for time to reach 100 cases
# ==============================================================================

# For each established epidemic, extract the time it first reaches 100 infections
threshold <- 100

time_to_threshold <- cuminf_df %>%
	filter(established == 1) %>%
	group_by(sim, kappa) %>%
	filter(cuminf >= threshold) %>%
	slice(1) %>%
	ungroup() %>%
	select(sim, kappa, tinf)

# Build the survival curve: for a grid of times, compute the fraction
# of epidemics that haven't yet reached the threshold
t_grid <- seq(0, max(time_to_threshold$tinf), length.out = 500)

survival_100_df <- time_to_threshold %>%
	group_by(kappa) %>%
	summarise(
		t = list(t_grid),
		prop_below = list(sapply(t_grid, function(tt) mean(tinf > tt))),
		.groups = "drop") %>%
	unnest(cols = c(t, prop_below))

breakvals = if(max(cuminf_df$tinf)>120){seq(0,365,by=15)}else{seq(0,365,by=7)}
fig_survival_100 <- ggplot(survival_100_df, aes(x=t, y=prop_below, col=kappa)) +
	geom_line(alpha=0.6, linewidth=0.8) +
	scale_x_continuous(breaks=breakvals) + 
	scale_color_manual(values=c("1"="black", "0.5"="blue", "0"="red")) + 
	labs(x = "Days",
	     y = paste0("Proportion not yet reaching ", threshold, " cases"),
	     col = "Kappa") +
	theme_classic()

