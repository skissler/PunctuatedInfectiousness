library(tidyverse)
library(odin)
source("code/utils.R")
source("code/parameters.R")

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
# 1. Global parameters
# ==============================================================================

popsize  <- 1000
nsim     <- 1000
psivals <- c(0, 0.5, 1)

# ==============================================================================
# 2. Loop over pathogens
# ==============================================================================

for (pars in parslist) {

pathogen <- pars$pathogen
T        <- pars$Tgen
popshape <- pars$popshape
poprate  <- pars$poprate
R0       <- pars$R0

cat(sprintf("\n===== %s: T=%.2f, popshape=%.2f, R0=%.1f =====\n",
            pathogen, T, popshape, R0))

# ==============================================================================
# 2a. Mean-field solution
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
# 2b. Stochastic simulations (cached)
# ==============================================================================

cache_file <- file.path("output", paste0("cuminf_df_", pathogen, ".csv"))

if (file.exists(cache_file)) {
	cat(sprintf("  %s: loading cached simulations from %s\n", pathogen, cache_file))
	cuminf_df <- read_csv(cache_file, show_col_types = FALSE) %>%
		mutate(psi = factor(psi))
} else {
	# Run the simulations and store the output
	cuminf_df <- tibble()
	for(sim in 1:nsim){
	    for(psi in psivals){
	            tinf <- sim_stochastic_fast(n=popsize,
	            	                        gen_inf_attempts=gen_inf_attempts_gamma(T, R0, popshape, psi))
	            cuminf_df <- bind_rows(cuminf_df,
	                    tibble(tinf=sort(tinf[tinf<Inf])) %>%
	                            mutate(cuminf=1:n(), sim=sim, psi=psi))
	    }
	    if(sim %% 100 == 0) cat(sprintf("  %s: sim %d/%d\n", pathogen, sim, nsim))
	}

	# Flag established epidemics (>= 10% of population infected)
	cuminf_df <- cuminf_df %>%
		mutate(psi = factor(psi)) %>%
		group_by(sim, psi) %>%
		mutate(established = as.integer(n() >= 0.1 * popsize)) %>%
		ungroup()

	# Save to cache
	write_csv(cuminf_df, cache_file)
	cat(sprintf("  %s: simulations saved to %s\n", pathogen, cache_file))
}


# Aggregate to daily/weekly resolution
dayjoin <- expand_grid(
	psi=unique(cuminf_df$psi),
	sim=1:nsim,
	day=0:max(ceiling(cuminf_df$tinf)))

weekjoin <- expand_grid(
	psi=unique(cuminf_df$psi),
	sim=1:nsim,
	week=0:max(ceiling(cuminf_df$tinf / 7)))

dailyinf_df <- cuminf_df %>%
	mutate(day = floor(tinf)) %>%
	group_by(psi, sim, day) %>%
	summarise(ninf = n(), established = first(established), .groups = "drop") %>%
	right_join(dayjoin, by = c("psi", "sim", "day")) %>%
	replace_na(list(ninf = 0)) %>%
	group_by(psi, sim) %>%
	arrange(day, .by_group = TRUE) %>%
	mutate(established = max(established, na.rm = TRUE)) %>%
	mutate(cuminf=cumsum(ninf))

weeklyinf_df <- cuminf_df %>%
	mutate(day = floor(tinf), week = floor(day / 7)) %>%
	group_by(psi, sim, week) %>%
	summarise(ninf = n(), established = first(established), .groups = "drop") %>%
	right_join(weekjoin, by = c("psi", "sim", "week")) %>%
	replace_na(list(ninf = 0)) %>%
	group_by(psi, sim) %>%
	arrange(week, .by_group = TRUE) %>%
	mutate(established = max(established, na.rm = TRUE)) %>%
	mutate(cuminf=cumsum(ninf))

lastday <- ceiling(max(cuminf_df$tinf))

dailyinf_means <- dailyinf_df %>%
	filter(established == 1) %>%
	group_by(psi, day) %>%
	summarise(mean_ninf = mean(ninf), mean_cuminf = mean(cuminf), .groups = "drop")

weeklyinf_means <- weeklyinf_df %>%
	filter(established == 1) %>%
	group_by(psi, week) %>%
	summarise(mean_ninf = mean(ninf), mean_cuminf = mean(cuminf), .groups = "drop")

# ==============================================================================
# 2c. Key metrics
# ==============================================================================

pest_table <- cuminf_df %>%
	group_by(sim, psi) %>%
	slice(1) %>%
	group_by(psi) %>%
	summarise(pestablished=sum(established)/nsim, .groups="drop")

fs_table <- cuminf_df %>%
	group_by(sim, psi) %>%
	filter(established==1) %>%
	summarise(finalsize=max(cuminf), .groups="drop") %>%
	group_by(psi) %>%
	summarise(fs_mean=mean(finalsize), fs_sd=sd(finalsize), .groups="drop")

cat(sprintf("  %s: P(established) by psi:\n", pathogen))
print(pest_table)
cat(sprintf("  %s: Final size summary:\n", pathogen))
print(fs_table)

# ==============================================================================
# 2d. Figures — epidemic trajectories
# ==============================================================================

# Cumulative stochastic epidemic curves (grey) with ODE overlay (blue)
fig_cuminf_overlay <- cuminf_df %>%
	ggplot(aes(x=tinf, y=cuminf, group=sim)) +
		geom_line(alpha=0.2, col="grey") +
		geom_line(data=filter(ren_out, t<=lastday),
			aes(x=t, y=cuminf*popsize),
			alpha=0.8, linewidth=1, col="black") +
		theme_classic() +
		labs(x="Time (days)", y="Cumulative infections", title = pathogen) +
		facet_wrap(~psi, nrow=1)

save_fig(fig_cuminf_overlay, paste0("fig_cuminf_overlay_", pathogen))

# Daily stochastic incidence curves (grey) with ODE overlay (blue)
fig_inf_overlay <- dailyinf_df %>%
	ggplot(aes(x=day, y=ninf, group=sim)) +
		geom_line(alpha=0.2, col="grey") +
		geom_line(data=filter(ren_daily, day<=lastday),
			aes(x=day, y=newinf*popsize),
			alpha=0.8, linewidth=1, col="black") +
		theme_classic() +
		labs(x="Time (days)", y="Daily new infections", title = pathogen) +
		facet_wrap(~psi, nrow=1)

save_fig(fig_inf_overlay, paste0("fig_inf_overlay_", pathogen))

# ==============================================================================
# 2e. Figure — KM curves for time to reach 100 cases
# ==============================================================================

# For each established epidemic, extract the time it first reaches 100 infections
threshold <- 100

time_to_threshold <- cuminf_df %>%
	filter(established == 1) %>%
	group_by(sim, psi) %>%
	filter(cuminf >= threshold) %>%
	slice(1) %>%
	ungroup() %>%
	select(sim, psi, tinf)

# Build the survival curve: for a grid of times, compute the fraction
# of epidemics that haven't yet reached the threshold
t_grid <- seq(0, max(time_to_threshold$tinf), length.out = 500)

survival_100_df <- time_to_threshold %>%
	group_by(psi) %>%
	summarise(
		t = list(t_grid),
		prop_below = list(sapply(t_grid, function(tt) mean(tinf > tt))),
		.groups = "drop") %>%
	unnest(cols = c(t, prop_below))

breakvals = if(max(cuminf_df$tinf)>120){seq(0,365,by=15)}else{seq(0,365,by=7)}
fig_survival_100 <- ggplot(survival_100_df, aes(x=t, y=prop_below, col=psi)) +
	geom_line(alpha=0.6, linewidth=0.8) +
	scale_x_continuous(breaks=breakvals) +
	scale_color_manual(values=c("1"="black", "0.5"="blue", "0"="red")) +
	labs(x = "Days",
	     y = paste0("Proportion not yet reaching ", threshold, " cases"),
	     col = expression(psi),
	     title = pathogen) +
	theme_classic()

save_fig(fig_survival_100, paste0("fig_survival_100_", pathogen))

cat(sprintf("  %s: figures saved.\n", pathogen))

} # end pathogen loop
