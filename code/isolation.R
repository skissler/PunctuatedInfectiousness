library(tidyverse)
library(patchwork)

source("code/utils.R")
source("code/parameters.R")

# ==============================================================================
# Isolation effects on TE and epidemic dynamics, across all pathogens
#
# Four analyses per pathogen:
#   1. Fixed isolation time — TE vs isolation timing relative to peak
#   2. Symptom-triggered isolation — TE vs mean symptom onset time
#   3. Regular screening — TE vs gap between tests (with & without action delay)
#   4. Epidemic simulations — fixed isolation vs naive R_eff = R0*(1-TE) adjust
#
# For sections 1-3, curves across psi values are overlaid; for section 4,
# stochastic epidemic trajectories and growth rates are compared between
# isolation and reduced-R0 controls.
# ==============================================================================

# ==============================================================================
# Global parameters (shared across pathogens)
# ==============================================================================

psi_vals_te <- c(0, 0.2, 0.5, 0.8, 1)  # for TE curves (sections 1-3)
psi_vals_sim <- c(0, 0.5, 1)           # for epidemic simulations (section 4)

# Symptom-triggered isolation
sigma_sym <- 0.5       # SD of symptom onset (days) relative to profile peak

# Regular screening
d_pre  <- 3            # days before peak at which test detects
d_post <- 7            # days after peak at which test stops detecting
w      <- d_pre + d_post
Delta_vals <- seq(from = 0.5, to = 14, by = 0.5)
lambda_act <- 1        # rate of exponential delay to action (mean 1 day)

# Epidemic simulations
tau_offset_sim <- -1   # isolation 1 day BEFORE peak (negative = before)
p_adhere       <- 0.5  # 50% of individuals adhere to isolation
popsize        <- 1000
nsim           <- 200
min_threshold    <- 10
growth_threshold <- 100

# Storage for composite figures
te_fixed_list         <- list()
te_symp_list          <- list()
te_testing_list       <- list()
te_testing_delay_list <- list()
pest_tables_list      <- list()

# ==============================================================================
# Loop over pathogens
# ==============================================================================

for (pars in parslist) {

pathogen <- pars$pathogen
T        <- pars$Tgen
alpha    <- pars$alpha
beta     <- pars$beta
R0       <- pars$R0

cat(sprintf("\n========== %s (R0=%g, alpha=%.2f, beta=%.3f) ==========\n",
            pathogen, R0, alpha, beta))

# Scale the tau_offset / mu_sym axis to cover ~2.5 SDs of the profile
sd_g             <- sqrt(alpha) / beta
tau_offset_range <- 2.5 * sd_g
tau_offset       <- seq(-tau_offset_range, tau_offset_range, length.out = 201)
mu_sym           <- seq(-tau_offset_range, tau_offset_range, length.out = 201)

# ==============================================================================
# Section 1. Fixed isolation TE curves
#
# TE(tau_off, psi) = 1 - F_Gamma(tau_off + mode; alpha*psi, beta)
# ==============================================================================

te_fixed_df <- expand_grid(psi = psi_vals_te, tau_offset = tau_offset) %>%
	mutate(mode = pmax(0, ((alpha*psi - 1)/beta))) %>%
	mutate(TE = 1 - pgamma(tau_offset + mode, alpha*psi, beta))

fig_te_fixed <- te_fixed_df %>%
	ggplot(aes(x = tau_offset, y = TE, col = factor(psi))) +
		geom_line(linewidth = 0.8, alpha = 0.8) +
		theme_classic() +
		labs(x = "Isolation time (days relative to peak infectiousness)",
		     y = "Test effectiveness",
		     col = expression(psi),
		     title = sprintf("%s  (R0 = %g)", pathogen, R0))

save_fig(fig_te_fixed, sprintf("fig_te_fixed_%s", pathogen))
te_fixed_list[[pathogen]] <- fig_te_fixed

# ==============================================================================
# Section 2. Symptom-triggered TE curves
#
# Symptom onset t ~ N(mu_sym, sigma_sym). Given t, isolation is immediate.
# TE(mu_sym, psi) = E_t[1 - F_Gamma(t + mode; alpha*psi, beta)]
# ==============================================================================

te_symp_df <- expand_grid(psi = psi_vals_te, mu_sym = mu_sym) %>%
	mutate(mode = pmax(0, ((alpha*psi - 1)/beta))) %>%
	rowwise() %>%
	mutate(TE = integrate(function(t)
		(1 - pgamma(t + mode, alpha*psi, beta)) *
		dnorm(t, mu_sym, sigma_sym),
		lower = mu_sym - 5*sigma_sym, upper = mu_sym + 5*sigma_sym)$value) %>%
	ungroup()

fig_te_symp <- te_symp_df %>%
	ggplot(aes(x = mu_sym, y = TE, col = factor(psi))) +
		geom_line(linewidth = 0.8, alpha = 0.8) +
		theme_classic() +
		labs(x = "Mean symptom onset time (days relative to peak infectiousness)",
		     y = "Test effectiveness",
		     col = expression(psi),
		     title = sprintf("%s  (R0 = %g, sigma_sym = %g)", pathogen, R0, sigma_sym))

save_fig(fig_te_symp, sprintf("fig_te_symp_%s", pathogen))
te_symp_list[[pathogen]] <- fig_te_symp

# ==============================================================================
# Section 3a. Regular testing TE (perfect sensitivity, no delay)
#
# First test falls uniformly in [0, Delta]; detection window length w.
# TE(Delta, psi) = (1/Delta) * int_0^min(Delta,w) S_psi(u - d_pre + mode) du
# ==============================================================================

te_testing_df <- expand_grid(psi = psi_vals_te, Delta = Delta_vals) %>%
	mutate(mode = pmax(0, ((alpha*psi - 1)/beta))) %>%
	rowwise() %>%
	mutate(TE = integrate(function(u)
		(1 - pgamma(u - d_pre + mode, alpha*psi, beta)),
		lower = 0, upper = min(Delta, w))$value / Delta) %>%
	ungroup()

fig_te_testing <- te_testing_df %>%
	ggplot(aes(x = Delta, y = TE, col = factor(psi))) +
		geom_line(linewidth = 0.8, alpha = 0.8) +
		theme_classic() +
		labs(x = "Gap between tests (days)",
		     y = "Test effectiveness",
		     col = expression(psi),
		     title = sprintf("%s  (R0 = %g, d_pre=%g, d_post=%g)",
		                     pathogen, R0, d_pre, d_post))

save_fig(fig_te_testing, sprintf("fig_te_testing_%s", pathogen))
te_testing_list[[pathogen]] <- fig_te_testing

# ==============================================================================
# Section 3b. Regular testing TE with exponential action delay
#
# Isolation at tau_det + delta_act, delta_act ~ Exp(lambda_act). Inner
# expectation E_delta[S_psi(t0 + delta)] closed form via Gamma-Exp convolution.
# ==============================================================================

te_testing_delay_df <- expand_grid(psi = psi_vals_te, Delta = Delta_vals) %>%
	mutate(mode = pmax(0, ((alpha*psi - 1)/beta))) %>%
	rowwise() %>%
	mutate(TE = {
		a  <- alpha*psi
		cr <- (beta / (beta + lambda_act))^a
		integrate(function(u) {
			t0 <- u - d_pre + mode
			(1 - pgamma(t0, a, beta)) -
			exp(lambda_act*t0) * cr * (1 - pgamma(t0, a, beta + lambda_act))
		}, lower = 0, upper = min(Delta, w))$value / Delta
	}) %>%
	ungroup()

fig_te_testing_delay <- te_testing_delay_df %>%
	ggplot(aes(x = Delta, y = TE, col = factor(psi))) +
		geom_line(linewidth = 0.8, alpha = 0.8) +
		theme_classic() +
		labs(x = "Gap between tests (days)",
		     y = "Test effectiveness",
		     col = expression(psi),
		     title = sprintf("%s  (R0 = %g, lambda_act = %g)",
		                     pathogen, R0, lambda_act))

save_fig(fig_te_testing_delay, sprintf("fig_te_testing_delay_%s", pathogen))
te_testing_delay_list[[pathogen]] <- fig_te_testing_delay

# ==============================================================================
# Section 4. Epidemic simulations: fixed isolation vs reduced R0
# ==============================================================================

# Expected TE for each psi used in the simulation
te_sim_df <- tibble(psi = psi_vals_sim) %>%
	mutate(mode = pmax(0, (alpha*psi - 1)/beta),
	       TE   = p_adhere * (1 - pgamma(tau_offset_sim + mode, alpha*psi, beta)))

cat(sprintf("  Expected TE by psi (tau_off = %g, p_adhere = %g):\n",
            tau_offset_sim, p_adhere))
print(te_sim_df)

results <- vector("list", nsim * length(psi_vals_sim) * 2)
idx <- 1L

for (sim in 1:nsim) {
	for (i in seq_along(psi_vals_sim)) {
		psi <- psi_vals_sim[i]
		te  <- te_sim_df$TE[i]

		# (1) Fixed isolation at tau_offset relative to peak
		tinf <- sim_stochastic_fast(n = popsize,
			gen_inf_attempts = gen_inf_attempts_gamma_fixed_iso(
				T, R0, alpha, psi,
				tau_offset = tau_offset_sim,
				p_adhere   = p_adhere))
		infected <- sort(tinf[is.finite(tinf)])
		results[[idx]] <- tibble(tinf = infected, cuminf = seq_along(infected),
		                         sim = sim, psi = psi, type = "isolation")
		idx <- idx + 1L

		# (2) Naive TE adjustment (same GI, reduced R0)
		tinf <- sim_stochastic_fast(n = popsize,
			gen_inf_attempts = gen_inf_attempts_gamma(
				T, R0*(1 - te), alpha, psi))
		infected <- sort(tinf[is.finite(tinf)])
		results[[idx]] <- tibble(tinf = infected, cuminf = seq_along(infected),
		                         sim = sim, psi = psi, type = "reduced_R0")
		idx <- idx + 1L
	}
	if (sim %% 50 == 0) cat(sprintf("  [%s] sim %d/%d\n", pathogen, sim, nsim))
}

cuminf_iso_df <- bind_rows(results)

# Flag established epidemics (>= 10% of population infected)
cuminf_iso_df <- cuminf_iso_df %>%
	group_by(sim, psi, type) %>%
	mutate(established = max(cuminf) >= 0.1 * popsize) %>%
	ungroup()

# Warn if too few established
est_counts <- cuminf_iso_df %>%
	group_by(sim, psi, type) %>%
	summarise(fs = max(cuminf), .groups = "drop") %>%
	mutate(established = fs >= 0.1 * popsize) %>%
	group_by(psi, type) %>%
	summarise(n_est = sum(established), .groups = "drop")

for (i in seq_len(nrow(est_counts))) {
	row <- est_counts[i, ]
	if (row$n_est == 0) {
		cat(sprintf("  WARNING: no epidemics established, psi=%s type=%s\n",
		            row$psi, row$type))
	} else if (row$n_est < 30) {
		cat(sprintf("  WARNING: only %d established, psi=%s type=%s\n",
		            row$n_est, row$psi, row$type))
	}
}

# --- Cumulative incidence overlay ---
fig_cuminf_iso <- cuminf_iso_df %>%
	filter(established) %>%
	ggplot(aes(x = tinf, y = cuminf,
	           group = interaction(sim, type), col = type)) +
		geom_line(alpha = 0.15, linewidth = 0.3) +
		facet_wrap(~psi, nrow = 1, labeller = label_both) +
		theme_classic() +
		labs(x = "Time (days)", y = "Cumulative infections", col = "Strategy",
		     title = sprintf("%s: isolation vs naive R_eff adjustment", pathogen)) +
		scale_color_manual(values = c("isolation" = "steelblue",
		                              "reduced_R0" = "tomato"))

save_fig(fig_cuminf_iso, sprintf("fig_cuminf_iso_%s", pathogen))

# --- P(establishment) and final size table ---
pest_iso_table <- cuminf_iso_df %>%
	group_by(sim, psi, type) %>%
	summarise(fs = max(cuminf), .groups = "drop") %>%
	group_by(psi, type) %>%
	summarise(p_est   = mean(fs >= 0.1 * popsize),
	          mean_fs = mean(fs[fs >= 0.1 * popsize]),
	          .groups = "drop") %>%
	mutate(pathogen = pathogen)

pest_tables_list[[pathogen]] <- pest_iso_table

cat(sprintf("  [%s] P(establishment) and mean final size:\n", pathogen))
print(pest_iso_table %>%
      pivot_wider(names_from = type, values_from = c(p_est, mean_fs)))

# --- Early growth rate comparison ---
growthrate_iso_df <- cuminf_iso_df %>%
	filter(established, cuminf >= min_threshold, cuminf <= growth_threshold) %>%
	group_by(sim, psi, type) %>%
	filter(n() >= 2) %>%
	summarise(growthrate = coef(lm(log(cuminf) ~ tinf))[2], .groups = "drop")

if (nrow(growthrate_iso_df) > 0) {
	fig_growthrate_iso <- growthrate_iso_df %>%
		ggplot(aes(x = growthrate, fill = type)) +
			geom_density(alpha = 0.4) +
			facet_wrap(~psi, nrow = 1, labeller = label_both) +
			theme_classic() +
			labs(x = "Early growth rate (per day)", y = "Density",
			     fill = "Strategy",
			     title = sprintf("%s: growth rate comparison", pathogen)) +
			scale_fill_manual(values = c("isolation" = "steelblue",
			                             "reduced_R0" = "tomato"))

	save_fig(fig_growthrate_iso, sprintf("fig_growthrate_iso_%s", pathogen))
}

cat(sprintf("  [%s] figures saved.\n", pathogen))

} # end pathogen loop

# ==============================================================================
# Composite figures across pathogens
# ==============================================================================

fig_te_fixed_composite <- wrap_plots(te_fixed_list, nrow = 1) +
	plot_annotation(title = "TE from fixed-time isolation, by pathogen")
save_fig(fig_te_fixed_composite, "fig_te_fixed", width = 14, height = 5)

fig_te_symp_composite <- wrap_plots(te_symp_list, nrow = 1) +
	plot_annotation(title = sprintf("TE from symptom-triggered isolation (sigma_sym = %g)", sigma_sym))
save_fig(fig_te_symp_composite, "fig_te_symp", width = 14, height = 5)

fig_te_testing_composite <- wrap_plots(te_testing_list, nrow = 1) +
	plot_annotation(title = sprintf("TE from regular testing (d_pre=%g, d_post=%g)", d_pre, d_post))
save_fig(fig_te_testing_composite, "fig_te_testing", width = 14, height = 5)

fig_te_testing_delay_composite <- wrap_plots(te_testing_delay_list, nrow = 1) +
	plot_annotation(title = sprintf("TE from regular testing w/ Exp(%g) delay", lambda_act))
save_fig(fig_te_testing_delay_composite, "fig_te_testing_delay", width = 14, height = 5)

# Combined establishment table
pest_combined <- bind_rows(pest_tables_list) %>%
	pivot_wider(names_from = type, values_from = c(p_est, mean_fs))

cat("\n===== Combined establishment / final size summary =====\n")
print(pest_combined)

cat("\nAll figures saved to figures/.\n")
