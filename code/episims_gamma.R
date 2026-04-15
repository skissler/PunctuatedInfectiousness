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

popsize  <- 1000 #10000
nsim     <- 1000 #2000
psivals <- c(0, 0.5, 1)

# ==============================================================================
# 2. Loop over pathogens
# ==============================================================================

for (pars in parslist) {

pathogen <- pars$pathogen
T        <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

cat(sprintf("\n===== %s: T=%.2f, alpha=%.2f, R0=%.1f =====\n",
            pathogen, T, alpha, R0))

# ==============================================================================
# 2a. Mean-field solution
# ==============================================================================

# Run the renewal equation model
ren_out <- renewal_epidemic(R0, alpha, T, 1000)

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

cuminf_df <- load_cache(pathogen, nsim, popsize, psivals)

if (is.null(cuminf_df)) {
	# Run the simulations and store the output
	cuminf_df <- tibble()
	for(sim in 1:nsim){
	    for(psi in psivals){
	            tinf <- sim_stochastic_fast(n=popsize,
	            	                        gen_inf_attempts=gen_inf_attempts_gamma(T, R0, alpha, psi))
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
	cache_file <- cache_path(pathogen, nsim, popsize)
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

# ==============================================================================
# 2f. Epidemic growth rate (log-linear regression on first 10-100 cases)
# ==============================================================================

# Theoretical growth rate: solve R0 * (beta/(beta+r))^alpha = 1  (Lotka-Euler)
r_malthusian <- uniroot(
	function(r) R0 * (beta / (beta + r))^alpha - 1,
	interval = c(0, 10))$root

# For each established epidemic, fit log(case_number) ~ time using the first
# 100 infection times. The slope is the empirical exponential growth rate.
growth_threshold <- 100
min_threshold <- 10

growthrate_df <- cuminf_df %>%
	filter(established == 1, cuminf >= min_threshold, cuminf <= growth_threshold) %>%
	group_by(sim, psi) %>%
	summarise(
		growthrate = coef(lm(log(cuminf) ~ tinf))[2],
		.groups = "drop")

growthrate_table <- growthrate_df %>%
	group_by(psi) %>%
	summarise(mean = mean(growthrate), sd = sd(growthrate), .groups = "drop")

cat(sprintf("  %s: theoretical growth rate = %.4f\n", pathogen, r_malthusian))
print(growthrate_table)

fig_growthrate_hists <- ggplot(growthrate_df, aes(x = growthrate)) +
	geom_histogram(aes(y = after_stat(density)), bins = 40,
	               fill = "white", col = "darkgrey") +
	geom_density(adjust = 2) +
	geom_vline(xintercept = r_malthusian, col = "blue", lty = "dashed", linewidth = 0.8) +
	geom_vline(data = growthrate_table, aes(xintercept = mean), col = "red", linewidth = 0.8) +
	theme_classic() +
	facet_wrap(~psi, nrow = 1) +
	labs(x = "Empirical growth rate (1/day)", y = "Density",
	     title = pathogen)

save_fig(fig_growthrate_hists, paste0("fig_growthrate_hists_", pathogen))

fig_growthrate_lines <- cuminf_df %>%
	filter(established == 1, cuminf <= growth_threshold, cuminf >= min_threshold) %>%
	ggplot(aes(x = tinf, y = log(cuminf), group = factor(sim))) +
		geom_point(alpha=0.1, size=0.2, col="grey") +
		geom_line(alpha=0.1, linewidth=0.2, col="grey") +
		geom_line(stat="smooth", method = "lm", alpha = 0.1, linewidth = 0.3) +
		geom_abline(intercept = 0, slope = r_malthusian,
		            col = "blue", linewidth = 0.8, lty = "dashed") +
		theme_classic() +
		facet_wrap(~psi, nrow = 1) +
		labs(x = "Time (days)", y = "log(cumulative infections)",
		     title = pathogen)

save_fig(fig_growthrate_lines, paste0("fig_growthrate_lines_", pathogen))

cat(sprintf("  %s: figures saved.\n", pathogen))

} # end pathogen loop
