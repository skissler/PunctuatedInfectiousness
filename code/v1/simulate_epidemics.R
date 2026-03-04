library(tidyverse)
library(odin)

source("code/utils.R")

# ==============================================================================
# Stochastic epidemic simulations: uncontrolled dynamics
#
# Corresponds to writeup sections:
#   - "The impact of the individual infectiousness profile on uncontrolled
#      epidemic dynamics"
#
# Runs nsim stochastic epidemics for each individual infectiousness profile
# (stepwise, smooth, spike) and compares trajectories, final sizes, and
# timing metrics against the mean-field ODE solution.
# ==============================================================================

# ==============================================================================
# 1. Parameters
# ==============================================================================

popsize <- 1000
e_dur   <- 2      # mean latent period (days); 0 gives SIR dynamics
i_dur   <- 3      # mean infectious period (days)
R0      <- 2
nsim    <- 1000   # number of stochastic replicates per profile

# ==============================================================================
# 2. Mean-field ODE solution
# ==============================================================================

# Choose SIR or SEIR depending on whether there is a latent period
if(e_dur > 0){
	ode_model <- seir$new(R0=R0, nu=1/e_dur, gamma=1/i_dur, N=popsize)
	print("Used SEIR for mean-field")
} else {
	ode_model <- sir$new(R0=R0, gamma=1/i_dur, N=popsize)
	print("Used SIR for mean-field")
}

t <- seq(0, 100, length.out = 1000)
ode_out <- as_tibble(data.frame(ode_model$run(t)))

# Aggregate ODE output to daily and weekly new-infection counts
# (ODE tracks proportions, so multiply by popsize for counts)
ode_daily <- ode_out %>%
	mutate(day = floor(t)) %>%
	group_by(day) %>%
	summarise(cuminf = max(cuminf)) %>%
	arrange(day) %>%
	mutate(newinf = cuminf - lag(cuminf)) %>%
	mutate(newinf = case_when(is.na(newinf) ~ cuminf, TRUE ~ newinf))

ode_weekly <- ode_out %>%
	mutate(day = floor(t), week = floor(day / 7)) %>%
	group_by(week) %>%
	summarise(cuminf = max(cuminf)) %>%
	arrange(week) %>%
	mutate(newinf = cuminf - lag(cuminf)) %>%
	mutate(newinf = case_when(is.na(newinf) ~ cuminf, TRUE ~ newinf))

# ==============================================================================
# 3. Stochastic simulations
# ==============================================================================

profiles <- c("stepwise", "smooth", "spike")

cuminf_df <- tibble()
for(sim in 1:nsim){
	for(profile in profiles){
		tinf <- sim_stochastic_fast(
			n=popsize, e_dur=e_dur, i_dur=i_dur, R0=R0, profiletype=profile
		)
		cuminf_df <- bind_rows(cuminf_df,
			tibble(tinf = sort(tinf[tinf < Inf])) %>%
				mutate(cuminf = 1:n(), sim = sim, profiletype = profile)
		)
	}
	if(sim %% 100 == 0) print(sim)
}

# Flag established epidemics (>= 10% of population infected)
cuminf_df <- cuminf_df %>%
	mutate(profiletype = factor(profiletype, levels = profiles)) %>%
	group_by(sim, profiletype) %>%
	mutate(established = as.integer(n() >= 0.1 * popsize)) %>%
	ungroup()

# ==============================================================================
# 4. Aggregate to daily / weekly resolution
# ==============================================================================

# Create complete grids so days/weeks with zero infections are represented
dayjoin <- expand_grid(
	profiletype = unique(cuminf_df$profiletype),
	sim         = 1:nsim,
	day         = 0:max(ceiling(cuminf_df$tinf))
)

weekjoin <- expand_grid(
	profiletype = unique(cuminf_df$profiletype),
	sim         = 1:nsim,
	week        = 0:max(ceiling(cuminf_df$tinf / 7))
)

dailyinf_df <- cuminf_df %>%
	mutate(day = floor(tinf)) %>%
	group_by(profiletype, sim, day) %>%
	summarise(ninf = n(), established = first(established), .groups = "drop") %>%
	right_join(dayjoin, by = c("profiletype", "sim", "day")) %>%
	replace_na(list(ninf = 0)) %>%
	group_by(profiletype, sim) %>%
	mutate(established = max(established, na.rm = TRUE))

weeklyinf_df <- cuminf_df %>%
	mutate(day = floor(tinf), week = floor(day / 7)) %>%
	group_by(profiletype, sim, week) %>%
	summarise(ninf = n(), established = first(established), .groups = "drop") %>%
	right_join(weekjoin, by = c("profiletype", "sim", "week")) %>%
	replace_na(list(ninf = 0)) %>%
	group_by(profiletype, sim) %>%
	mutate(established = max(established, na.rm = TRUE))

# ==============================================================================
# 5. Figures — epidemic trajectories
# ==============================================================================

# One example established trajectory highlighted in black; ODE in blue
fig_cuminf_overlay <- cuminf_df %>%
	ggplot(aes(x = tinf, y = cuminf, group = sim)) +
		geom_line(alpha = 0.2, col = "grey") +
		geom_line(data = cuminf_df %>%
			filter(established == 1) %>%
			group_by(profiletype) %>%
			filter(sim == min(sim)),
			linewidth = 0.2, alpha = 0.8, col = "black") +
		geom_line(data = ode_daily,
			aes(x = day, y = cuminf * popsize),
			alpha = 0.6, linewidth = 1, col = "blue") +
		theme_classic() +
		facet_wrap(~profiletype, nrow = 1)

# Mean +/- 1 SD of cumulative infections (established epidemics only)
fig_cuminf_means <- dailyinf_df %>%
	filter(established == 1) %>%
	arrange(profiletype, sim, day) %>%
	mutate(cumninf = cumsum(ninf)) %>%
	group_by(profiletype, day) %>%
	summarise(cumninf_mean = mean(cumninf), cumninf_sd = sd(cumninf),
	          .groups = "drop") %>%
	ggplot(aes(x = day, y = cumninf_mean, col = profiletype)) +
		geom_line() +
		geom_ribbon(aes(ymin = cumninf_mean - cumninf_sd,
		                ymax = cumninf_mean + cumninf_sd,
		                fill = profiletype), alpha = 0.2) +
		geom_line(data = ode_daily,
			aes(x = day, y = cuminf * popsize),
			alpha = 0.6, linewidth = 1, col = "magenta", lty = "dashed") +
		theme_classic() +
		scale_color_manual(values = c("black", "blue", "red")) +
		scale_fill_manual(values = c("black", "blue", "red")) +
		facet_wrap(~profiletype)

# Daily incidence overlay
fig_inf_overlay <- dailyinf_df %>%
	ggplot(aes(x = day, y = ninf, group = sim)) +
		geom_line(alpha = 0.2, col = "grey") +
		geom_line(data = dailyinf_df %>%
			filter(established == 1) %>%
			group_by(profiletype) %>%
			filter(sim == min(sim)),
			linewidth = 0.2, alpha = 0.8, col = "black") +
		geom_line(data = ode_daily,
			aes(x = day, y = newinf * popsize),
			alpha = 0.6, linewidth = 1, col = "blue") +
		theme_classic() +
		facet_wrap(~profiletype, nrow = 1)

# Weekly incidence overlay
fig_inf_wk_overlay <- weeklyinf_df %>%
	ggplot(aes(x = week, y = ninf, group = sim)) +
		geom_line(alpha = 0.2, col = "grey") +
		geom_line(data = weeklyinf_df %>%
			filter(established == 1) %>%
			group_by(profiletype) %>%
			filter(sim == min(sim)),
			linewidth = 0.2, alpha = 0.8, col = "black") +
		geom_line(data = ode_weekly,
			aes(x = week, y = newinf * popsize),
			alpha = 0.6, linewidth = 1, col = "blue") +
		theme_classic() +
		facet_wrap(~profiletype, nrow = 1)

# Weekly incidence on log scale
weeklyinf_df %>%
	ggplot(aes(x = week, y = log(ninf + 1), group = sim)) +
		geom_line(alpha = 0.2, col = "grey") +
		geom_line(data = weeklyinf_df %>%
			filter(established == 1) %>%
			group_by(profiletype) %>%
			filter(sim == min(sim)),
			linewidth = 0.2, alpha = 0.8, col = "black") +
		geom_line(data = ode_weekly,
			aes(x = week, y = log(newinf * popsize + 1)),
			alpha = 0.6, linewidth = 1, col = "blue") +
		theme_classic() +
		facet_wrap(~profiletype, nrow = 1)

# ==============================================================================
# 6. Summary statistics
# ==============================================================================

# -- Establishment probability -------------------------------------------------
pest_table <- cuminf_df %>%
	group_by(sim, profiletype) %>%
	slice(1) %>%
	group_by(profiletype) %>%
	summarise(pestablished = sum(established) / nsim, .groups = "drop")

# -- Week-to-week variation (log-scale first differences) ----------------------
weekvar_table <- weeklyinf_df %>%
	filter(established == 1) %>%
	arrange(profiletype, sim, week) %>%
	group_by(profiletype, sim) %>%
	mutate(logninfdiff = log(ninf + 1) - lag(log(ninf + 1))) %>%
	group_by(profiletype) %>%
	summarise(logninfdiff_mean = mean(logninfdiff, na.rm = TRUE),
	          logninfdiff_sd   = sd(logninfdiff, na.rm = TRUE),
	          .groups = "drop")

# -- Final size ----------------------------------------------------------------
fs_df <- cuminf_df %>%
	group_by(profiletype, sim) %>%
	summarise(finalsize = max(cuminf), .groups = "drop")

# All epidemics (including die-outs)
fs_df %>%
	ggplot(aes(x = finalsize)) +
		geom_histogram(bins = 50, fill = "white", col = "darkgrey") +
		theme_classic() +
		facet_wrap(~profiletype, nrow = 1)

# Established epidemics only (final size > 100)
fig_fs <- fs_df %>%
	filter(finalsize > 100) %>%
	ggplot(aes(x = finalsize)) +
		geom_histogram(aes(y = after_stat(density)), bins = 50,
		               fill = "white", col = "darkgrey") +
		theme_classic() +
		facet_wrap(~profiletype, nrow = 1)

fs_table <- fs_df %>%
	filter(finalsize > 100) %>%
	group_by(profiletype) %>%
	summarise(fs_mean = mean(finalsize), fs_sd = sd(finalsize), .groups = "drop")

# ==============================================================================
# 7. Figures — timing metrics
# ==============================================================================

# Time to 100 infections
fig_t100 <- cuminf_df %>%
	filter(cuminf == 100) %>%
	ggplot(aes(x = tinf)) +
		geom_histogram(aes(y = after_stat(density)), bins = 50,
		               fill = "white", col = "darkgrey") +
		geom_density(adjust = 2) +
		theme_classic() +
		facet_wrap(~profiletype, nrow = 1) +
		labs(title = "Time to 100 infections")

t100_table <- cuminf_df %>%
	filter(cuminf == 100) %>%
	group_by(profiletype) %>%
	summarise(t100_mean = mean(tinf), t100_sd = sd(tinf), .groups = "drop")

# Time from 100 to 250 infections
fig_t100_250 <- cuminf_df %>%
	filter(cuminf == 100 | cuminf == 250) %>%
	pivot_wider(names_from = cuminf, values_from = tinf) %>%
	mutate(tdiff = `250` - `100`) %>%
	ggplot(aes(x = tdiff)) +
		geom_histogram(aes(y = after_stat(density)), bins = 50,
		               fill = "white", col = "darkgrey") +
		geom_density(adjust = 2) +
		theme_classic() +
		facet_wrap(~profiletype, nrow = 1) +
		labs(title = "Time from 100 to 250 infections")

# Time to epidemic die-out
fig_t_dieout <- cuminf_df %>%
	group_by(sim, profiletype) %>%
	summarise(tmax = max(tinf), .groups = "drop") %>%
	ggplot(aes(x = tmax)) +
		geom_histogram(aes(y = after_stat(density)), bins = 50,
		               fill = "white", col = "darkgrey") +
		theme_classic() +
		facet_wrap(~profiletype, nrow = 1) +
		labs(title = "Time to epidemic die-out", x = "Time (days)")

# Peak daily infections (established epidemics)
fig_peak <- dailyinf_df %>%
	filter(established == 1) %>%
	group_by(sim, profiletype) %>%
	summarise(peak = max(ninf), .groups = "drop") %>%
	ggplot(aes(x = peak)) +
		geom_histogram(aes(y = after_stat(density)), bins = 50,
		               fill = "white", col = "darkgrey") +
		theme_classic() +
		facet_wrap(~profiletype, nrow = 1) +
		labs(title = "Peak daily number of infections", x = "Number of infections")

# ==============================================================================
# 8. Save figures
# ==============================================================================

ggsave("figures/fig_cuminf_overlay.pdf", fig_cuminf_overlay, width = 14, height = 5)
ggsave("figures/fig_cuminf_means.pdf", fig_cuminf_means, width = 14, height = 5)
ggsave("figures/fig_inf_overlay.pdf", fig_inf_overlay, width = 14, height = 5)
ggsave("figures/fig_fs.pdf", fig_fs, width = 14, height = 5)
ggsave("figures/fig_t100.pdf", fig_t100, width = 14, height = 5)
ggsave("figures/fig_peak.pdf", fig_peak, width = 14, height = 5)

cat("\n=== Epidemic simulation figures saved ===\n")
cat("\nEstablishment probabilities:\n")
print(pest_table)
cat("\nFinal size (established):\n")
print(fs_table)
cat("\nTime to 100 infections:\n")
print(t100_table)
