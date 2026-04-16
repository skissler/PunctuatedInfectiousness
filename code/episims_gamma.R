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
# intermediate (psi=0.5) individual infectiousness profiles, and:
#   - Plots epidemic trajectories (daily incidence and cumulative incidence)
#   - Calculates extinction probabilities and final size distributions
#   - Plots the time to reach various epidemic milestones (100 infections,
#     peak) as Kaplan-Meier curves.
#   - Calculates empirical epidemic growth rates
# ==============================================================================

# ==============================================================================
# 1. Global parameters
# ==============================================================================

popsize       <- 10000 #10000
nsim          <- 1000 #2000
psivals       <- c(0, 0.5, 1)
max_plot_sims <- 1000  # max number of individual trajectories to draw on plots

# Thresholds for growth rate estimation and survival curves
threshold        <- 100
min_threshold    <- 10
growth_threshold <- 100

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
ren_out <- renewal_epidemic(R0, alpha, T, popsize)

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
# 2b. Stochastic simulations (cached, incremental processing)
# ==============================================================================

cache <- load_cache_v2(pathogen, nsim, popsize, psivals)

if (!is.null(cache)) {
	sim_summary_df <- cache$summary
	plot_df        <- cache$plot
} else {
	# Run simulations with incremental processing:
	#   - Compute per-sim summary immediately, discard raw infection times
	#   - Keep full trajectories only for sim <= max_plot_sims (for plotting)
	summary_list <- vector("list", nsim * length(psivals))
	plot_list    <- vector("list", max_plot_sims * length(psivals))
	si <- 0L  # summary index
	pi <- 0L  # plot index

	for(sim in 1:nsim){
	    for(psi in psivals){
	        tinf <- sim_stochastic_fast(n=popsize,
	                                    gen_inf_attempts=gen_inf_attempts_gamma(T, R0, alpha, psi))

	        # Extract sorted finite infection times
	        infected <- sort(tinf[is.finite(tinf)])
	        final_size <- length(infected)
	        established <- as.integer(final_size >= 0.1 * popsize)

	        # Compute time to threshold (100 cases)
	        time_to_100 <- if(final_size >= threshold) infected[threshold] else NA_real_

	        # Compute growth rate
	        growthrate <- compute_growth_rate(infected, min_threshold, growth_threshold)

	        # Store summary row
	        si <- si + 1L
	        summary_list[[si]] <- tibble(
	            sim = sim, psi = psi,
	            established = established,
	            final_size = final_size,
	            time_to_100 = time_to_100,
	            growthrate = growthrate
	        )

	        # Keep full trajectory for plotting subset
	        if(sim <= max_plot_sims && final_size > 0){
	            pi <- pi + 1L
	            plot_list[[pi]] <- tibble(
	                tinf = infected,
	                cuminf = seq_along(infected),
	                sim = sim,
	                psi = psi,
	                established = established
	            )
	        }
	    }
	    if(sim %% 100 == 0) cat(sprintf("  %s: sim %d/%d\n", pathogen, sim, nsim))
	}

	sim_summary_df <- bind_rows(summary_list[1:si])
	plot_df        <- bind_rows(plot_list[seq_len(pi)])

	# Save to v2 cache
	write_csv(sim_summary_df, cache_path_summary(pathogen, nsim, popsize))
	write_csv(plot_df, cache_path_plot(pathogen, nsim, popsize))
	cat(sprintf("  %s: simulations saved (v2 cache)\n", pathogen))

	sim_summary_df <- sim_summary_df %>% mutate(psi = factor(psi))
	plot_df        <- plot_df %>% mutate(psi = factor(psi))
}

# ==============================================================================
# 2c. Aggregate plot subset to daily/weekly resolution (for trajectory plots)
# ==============================================================================

lastday <- ceiling(max(plot_df$tinf))

dayjoin <- expand_grid(
	psi=unique(plot_df$psi),
	sim=unique(plot_df$sim),
	day=0:lastday)

dailyinf_df <- plot_df %>%
	mutate(day = floor(tinf)) %>%
	group_by(psi, sim, day) %>%
	summarise(ninf = n(), established = first(established), .groups = "drop") %>%
	right_join(dayjoin, by = c("psi", "sim", "day")) %>%
	replace_na(list(ninf = 0)) %>%
	group_by(psi, sim) %>%
	arrange(day, .by_group = TRUE) %>%
	mutate(established = max(established, na.rm = TRUE)) %>%
	mutate(cuminf=cumsum(ninf))

# ==============================================================================
# 2d. Key metrics (from summary)
# ==============================================================================

pest_table <- sim_summary_df %>%
	group_by(psi) %>%
	summarise(pestablished=mean(established), .groups="drop")

fs_table <- sim_summary_df %>%
	filter(established==1) %>%
	group_by(psi) %>%
	summarise(fs_mean=mean(final_size), fs_sd=sd(final_size), .groups="drop")

cat(sprintf("  %s: P(established) by psi:\n", pathogen))
print(pest_table)
cat(sprintf("  %s: Final size summary:\n", pathogen))
print(fs_table)

# ==============================================================================
# 2e. Figures — epidemic trajectories (plot subset only)
# ==============================================================================

# Cumulative stochastic epidemic curves (grey) with ODE overlay (black)
fig_cuminf_overlay <- plot_df %>%
	ggplot(aes(x=tinf, y=cuminf, group=sim)) +
		geom_line(alpha=0.2, col="grey") +
		geom_line(data=filter(ren_out, t<=lastday),
			aes(x=t, y=cuminf*popsize, group=1),
			alpha=0.8, linewidth=1, col="black") +
		theme_classic() +
		labs(x="Time (days)", y="Cumulative infections", title = pathogen) +
		facet_wrap(~psi, nrow=1)

save_fig(fig_cuminf_overlay, paste0("fig_cuminf_overlay_", pathogen))

# Cumulative curves with time-to-threshold milestone annotations
milestone_df <- plot_df %>%
	filter(established == 1, cuminf >= threshold) %>%
	group_by(sim, psi) %>%
	slice(1) %>%
	group_by(psi) %>%
	summarise(
		median_t = median(tinf),
		q05_t    = quantile(tinf, 0.05),
		q95_t    = quantile(tinf, 0.95),
		.groups  = "drop")

fig_cuminf_milestones <- plot_df %>%
	ggplot(aes(x=tinf, y=cuminf, group=sim)) +
		geom_line(alpha=0.2, col="grey") +
		geom_line(data=filter(ren_out, t<=lastday),
			aes(x=t, y=cuminf*popsize, group=1),
			alpha=0.8, linewidth=1, col="black") +
		geom_segment(data=milestone_df,
			aes(x=q05_t, xend=q05_t,
			    y=threshold*0.85, yend=threshold*1.15, group=1),
			col="red", linewidth=0.6) +
		geom_segment(data=milestone_df,
			aes(x=q95_t, xend=q95_t,
			    y=threshold*0.85, yend=threshold*1.15, group=1),
			col="red", linewidth=0.6) +
		geom_point(data=milestone_df,
			aes(x=median_t, y=threshold, group=1),
			col="red", size=2.5) +
		theme_classic() +
		labs(x="Time (days)", y="Cumulative infections", title = pathogen) +
		facet_wrap(~psi, nrow=1)

save_fig(fig_cuminf_milestones, paste0("fig_cuminf_milestones_", pathogen))

# Daily stochastic incidence curves (grey) with ODE overlay (black)
fig_inf_overlay <- dailyinf_df %>%
	ggplot(aes(x=day, y=ninf, group=sim)) +
		geom_line(alpha=0.2, col="grey") +
		geom_line(data=filter(ren_daily, day<=lastday),
			aes(x=day, y=newinf*popsize, group=1),
			alpha=0.8, linewidth=1, col="black") +
		theme_classic() +
		labs(x="Time (days)", y="Daily new infections", title = pathogen) +
		facet_wrap(~psi, nrow=1)

save_fig(fig_inf_overlay, paste0("fig_inf_overlay_", pathogen))

# ==============================================================================
# 2f. Figure — KM curves for time to reach 100 cases (from summary)
# ==============================================================================

time_to_threshold <- sim_summary_df %>%
	filter(established == 1, !is.na(time_to_100)) %>%
	select(sim, psi, time_to_100)

t_grid <- seq(0, max(time_to_threshold$time_to_100), length.out = 500)

survival_100_df <- time_to_threshold %>%
	group_by(psi) %>%
	summarise(
		t = list(t_grid),
		prop_below = list(sapply(t_grid, function(tt) mean(time_to_100 > tt))),
		.groups = "drop") %>%
	unnest(cols = c(t, prop_below))

breakvals = if(max(time_to_threshold$time_to_100)>120){seq(0,365,by=15)}else{seq(0,365,by=7)}
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
# 2g. Epidemic growth rate (from summary)
# ==============================================================================

# Theoretical growth rate: solve R0 * (beta/(beta+r))^alpha = 1  (Lotka-Euler)
r_malthusian <- uniroot(
	function(r) R0 * (beta / (beta + r))^alpha - 1,
	interval = c(0, 10))$root

growthrate_df <- sim_summary_df %>%
	filter(established == 1, !is.na(growthrate))

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

# Growth rate lines plot: uses plot subset trajectories
growth_incidence <- plot_df %>%
	filter(established == 1, cuminf >= min_threshold, cuminf <= growth_threshold) %>%
	mutate(day = floor(tinf)) %>%
	group_by(sim, psi, day) %>%
	summarise(count = n(), .groups = "drop") %>%
	group_by(sim, psi) %>%
	complete(day = seq(min(day), max(day)), fill = list(count = 0)) %>%
	mutate(day0 = day - min(day)) %>%
	ungroup()

fig_growthrate_lines <- growth_incidence %>%
	filter(count > 0) %>%
	ggplot(aes(x = day0, y = count, group = factor(sim))) +
		geom_point(alpha = 0.1, size = 0.3, col = "grey") +
		geom_abline(intercept = log10(min_threshold),
		            slope = r_malthusian / log(10),
		            col = "blue", linewidth = 0.8, lty = "dashed") +
		scale_y_log10() +
		theme_classic() +
		facet_wrap(~psi, nrow = 1) +
		labs(x = "Days since case 10", y = "Daily incidence",
		     title = pathogen)

save_fig(fig_growthrate_lines, paste0("fig_growthrate_lines_", pathogen))

cat(sprintf("  %s: figures saved.\n", pathogen))

} # end pathogen loop
