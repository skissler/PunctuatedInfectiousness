library(tidyverse)
source("code/utils.R")
source("code/parameters.R")

# ==============================================================================
# Infinite-population growth rate analysis
#
# Runs infinite-population stochastic simulations (no susceptible depletion)
# to estimate exponential growth rates under different punctuation (psi) values,
# and compares empirical growth rates to the Lotka-Euler theoretical prediction.
# ==============================================================================

# ==============================================================================
# 1. Global parameters
# ==============================================================================

psivals              <- c(0, 0.5, 1)

# Growth rate estimation thresholds
min_growth_threshold <- 100
max_growth_threshold <- 500

# Infinite-population simulation parameters
max_cases_infpop     <- 1000   # cases per sim (must be > max_growth_threshold)
nsim_infpop          <- 1000   # number of sims
max_plot_sims_infpop <- 250    # trajectories to keep for plotting

# ==============================================================================
# 2. Loop over pathogens
# ==============================================================================

for (pars in parslist) {

pathogen <- pars$pathogen
Tgen     <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

cat(sprintf("\n===== %s: growth rate analysis =====\n", pathogen))

# Theoretical growth rate: solve R0 * (beta/(beta+r))^alpha = 1  (Lotka-Euler)
r_malthusian <- uniroot(
	function(r) R0 * (beta / (beta + r))^alpha - 1,
	interval = c(0, 10))$root

cat(sprintf("  %s: theoretical growth rate = %.4f\n", pathogen, r_malthusian))

# ==============================================================================
# 2a. Infinite-population stochastic simulations (cached)
# ==============================================================================

cache_infpop <- load_cache_infpop(pathogen, nsim_infpop, max_cases_infpop, psivals)

if (!is.null(cache_infpop)) {
	infpop_summary_df <- cache_infpop$summary
	infpop_plot_df    <- cache_infpop$plot
} else {
	summary_list_ip <- vector("list", nsim_infpop * length(psivals))
	plot_list_ip    <- vector("list", max_plot_sims_infpop * length(psivals))
	summindex_ip <- 0L
	plotindex_ip <- 0L

	for (sim in 1:nsim_infpop) {
	    for (psi in psivals) {
	        infection_times <- sim_infinite_pop(
	            max_cases = max_cases_infpop,
	            gen_inf_attempts = gen_inf_attempts_gamma(Tgen, R0, alpha, psi))

	        n_infected <- length(infection_times)

	        # Compute growth rate
	        growthrate <- compute_growth_rate(infection_times, min_growth_threshold, max_growth_threshold)

	        # Store summary row
	        summindex_ip <- summindex_ip + 1L
	        summary_list_ip[[summindex_ip]] <- tibble(
	            sim = sim, psi = psi,
	            n_infected = n_infected,
	            growthrate = growthrate
	        )

	        # Keep full trajectory for plotting subset
	        if (sim <= max_plot_sims_infpop && n_infected > 0) {
	            plotindex_ip <- plotindex_ip + 1L
	            plot_list_ip[[plotindex_ip]] <- tibble(
	                tinf = infection_times,
	                cuminf = seq_along(infection_times),
	                sim = sim,
	                psi = psi
	            )
	        }
	    }
	    if (sim %% 100 == 0) cat(sprintf("  %s [infpop]: sim %d/%d\n", pathogen, sim, nsim_infpop))
	}

	infpop_summary_df <- bind_rows(summary_list_ip[1:summindex_ip])
	infpop_plot_df    <- bind_rows(plot_list_ip[seq_len(plotindex_ip)])

	# Save cache
	write_csv(infpop_summary_df,
	          cache_path_infpop_summary(pathogen, nsim_infpop, max_cases_infpop))
	write_csv(infpop_plot_df,
	          cache_path_infpop_plot(pathogen, nsim_infpop, max_cases_infpop))
	cat(sprintf("  %s: infpop simulations saved\n", pathogen))

	infpop_summary_df <- infpop_summary_df %>% mutate(psi = factor(psi))
	infpop_plot_df    <- infpop_plot_df %>% mutate(psi = factor(psi))
}

# ==============================================================================
# 2b. Growth rate analysis and figures
# ==============================================================================

infpop_growthrate_df <- infpop_summary_df %>%
	filter(!is.na(growthrate))

infpop_growthrate_table <- infpop_growthrate_df %>%
	group_by(psi) %>%
	summarise(mean = mean(growthrate), sd = sd(growthrate), .groups = "drop")

cat(sprintf("  %s [infpop]: theoretical growth rate = %.4f\n", pathogen, r_malthusian))
print(infpop_growthrate_table)

# Histogram of growth rates by psi
fig_growthrate_infpop_hists <- ggplot(infpop_growthrate_df, aes(x = growthrate)) +
	geom_histogram(aes(y = after_stat(density)), bins = 40,
	               fill = "white", col = "darkgrey") +
	geom_density(adjust = 2) +
	geom_vline(xintercept = r_malthusian, col = "blue", lty = "dashed", linewidth = 0.8) +
	geom_vline(data = infpop_growthrate_table, aes(xintercept = mean),
	           col = "red", linewidth = 0.8) +
	theme_classic() +
	facet_wrap(~psi, nrow = 1) +
	labs(x = "Empirical growth rate (1/day)", y = "Density",
	     title = paste0(pathogen, " (infinite pop)"))

save_fig(fig_growthrate_infpop_hists, paste0("fig_growthrate_infpop_hists_", pathogen))

# Growth rate lines: daily incidence on log scale in the growth window
infpop_growth_incidence <- infpop_plot_df %>%
	filter(cuminf >= min_growth_threshold, cuminf <= max_growth_threshold) %>%
	mutate(day = floor(tinf)) %>%
	group_by(sim, psi, day) %>%
	summarise(count = n(), .groups = "drop") %>%
	group_by(sim, psi) %>%
	complete(day = seq(min(day), max(day)), fill = list(count = 0)) %>%
	mutate(day0 = day - min(day)) %>%
	filter(day0 > 0, day0 < max(day0)) %>%
	ungroup()

fig_growthrate_infpop_lines <- infpop_growth_incidence %>%
	filter(count > 0) %>%
	ggplot(aes(x = day0, y = count, group = factor(sim))) +
		geom_line(alpha = 0.1, linewidth = 0.3, col = "grey") +
		geom_point(alpha = 0.2, size = 0.3, col = "grey") +
		geom_abline(intercept = log10(r_malthusian * min_growth_threshold),
		            slope = r_malthusian / log(10),
		            col = "blue", linewidth = 0.8, lty = "dashed") +
		scale_y_log10() +
		theme_classic() +
		facet_wrap(~psi, nrow = 1) +
		labs(x = sprintf("Days since case %d", min_growth_threshold),
		     y = "Daily incidence",
		     title = paste0(pathogen, " (infinite pop)"))

save_fig(fig_growthrate_infpop_lines, paste0("fig_growthrate_infpop_lines_", pathogen))

cat(sprintf("  %s: figures saved.\n", pathogen))

} # end pathogen loop
