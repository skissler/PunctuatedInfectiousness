library(tidyverse) 

source("code/utils.R")
source("code/parameters.R") 

# ==============================================================================
# Fixed isolation time 
# ==============================================================================

# For psi \in {0, 0.2, 0.5, 0.8, 1}, plot S_\psi(\tau_{iso}), i.e., the remaining infectiousness as a function of isolation time. S_\psi(\tau_{iso}) should just be a Gamma(alpha*psi, beta) CDF, so this should be easy to compute. The lines for each psi value should be overlaid in a single plot. The goal of the plot is to demonstrate how, for a given isolation time, TE can vary widely depending on how spiky (psi -> 0) the infectiousness profile is. It has the added benefit of showing that the equivalence points for the curves pass through a common point where TE=0.5 and \tau_{iso} > 0, highlighting that for right-skewed infectiousness profils, we gain a bit of time after rpeak infectiousness where spike profiles still have higher TE (until we get to really incredibly spiked profiles when the equivalence breaks down). A rough sketch of how this might look is hashed out in code/te.R.

# taum means "tau mode", i.e., time indexed by the mode of the infectiousness profile rather than the time of infection (which is indexed by tau)

taum_iso <- seq(from=-5, to=5, by=0.1)
psi_vals <- c(0, 0.2, 0.5, 0.8, 1)
alpha <- parslist[[1]]$alpha
beta <- parslist[[1]]$beta

S_deterministic_df <- expand_grid(psi=psi_vals, taum_iso=taum_iso) %>% 
	mutate(mode=pmax(0,((alpha*psi-1)/beta))) %>% 
	mutate(S=1-pgamma(taum_iso+mode, alpha*psi, beta))

fig_S_deterministic <- S_deterministic_df %>% 
	ggplot(aes(x=taum_iso, y=S, col=factor(psi))) + 
		geom_line(linewidth=0.8, alpha=0.8) + 
		theme_classic() + 
		labs(x="Isolation time (days relative to peak infectiousness)", y="Test effectiveness", col="psi")

# ==============================================================================
# Variable isolation time (symptoms)
# ==============================================================================

# Next, we can consider what happens when the detection time is a random variable, rather than fixed, and when there are potentially delays between detection and isolation. We'll begin with considering symptom onset, which follows a normal distribution with mean \mu_{symp} (relative to the peak of the infectiousness profile, so \mu_{symp} = 0 means the timing of symptoms corresponds with peak infectiousness on average, while \mu_{symp} < 0 means symptoms precede peak infectiousness) and variance \sigma^2_{symp}. The first task is to generate TE curves like in the previous section, but for the variable detection time. We assume initially that isolation is immediate. Then, we will generate another set of curves assumin isolation follows some exponentially-distributed waiting time after detection. We will use these plots to see how variation in the onset of symptoms and the delay before taking action impact TE as the infectiousness profile width varies from smooth to spike. The TE curves will have TE on the vertical axis and the mean symptom onset time on the horizontal axis. We can maybe consider \mu_{symp} between -3 and 3, and let's set \sigma^2_{symp} = 0.5 for now. 

taum_iso <- seq(from=-5, to=5, by=0.1)
psi_vals <- c(0, 0.2, 0.5, 0.8, 1)
alpha <- parslist[[1]]$alpha
beta <- parslist[[1]]$beta
sigma_det <- 0.5

S_symp_df <- expand_grid(psi=psi_vals, taum_iso=taum_iso) %>%
	mutate(mode=pmax(0,((alpha*psi-1)/beta))) %>%
	rowwise() %>%
	mutate(S=integrate(function(t)
		(1-pgamma(t+mode, alpha*psi, beta))*
		dnorm(t, taum_iso, sigma_det),
		lower=taum_iso-5*sigma_det, upper=taum_iso+5*sigma_det)$value) %>%
	ungroup()

fig_S_symp <- S_symp_df %>% 
	ggplot(aes(x=taum_iso, y=S, col=factor(psi))) + 
		geom_line(linewidth=0.8, alpha=0.8) + 
		theme_classic() + 
		labs(x="Mean symptom onset time (days relative to peak infectiousness)", y="Test effectiveness", col="psi")

# ==============================================================================
# Variable isolation time (regular testing)
# ==============================================================================

# Next, we will create yet another series of TE curves, this time with TE on the vertical axis and gap-between-tests on the horizontal axis. The model here is that there's a window of detectability extending -d_{pre} before peak infectiousness and d_{post} after peak infectiousness. Testing is regular (i.e., there are exactly \Delta days between each test), with a uniform random offset so that the testing phase is random with respect to the infectiousness profile. Parameters to consider here include: test sensitivity, delay to action (like before, modeled as a exponential random variable), the gap between tests \Delta, and the length of d_{pre} and d_{post}. For reasonable values of d_{pre} and d_{post}, e.g., d_{pre} = 3 and d_{post} = 7, we can initially consider a perfectly sensitive test and see how TE varies with the gap between tests. Then, we should see how delays between testing and taking action affect things. We can also consider varying the test positivity window. Each plot should show multiple curves depicting how the punctuatedness of the individual infectiousness profile \psi impacts TE as a function of the gap between tests.

d_pre <- 3
d_post <- 7
w <- d_pre + d_post
Delta_vals <- seq(from=0.5, to=14, by=0.5)
psi_vals <- c(0, 0.2, 0.5, 0.8, 1)
alpha <- parslist[[1]]$alpha
beta <- parslist[[1]]$beta

# --- Perfect sensitivity, no delay to action ---
# The first test in the window falls at offset U ~ Uniform(0, Delta) from the
# window start (mode - d_pre), by the random-phase argument. With p_sens = 1,
# detection occurs at the first test in the window (if any).
# TE = (1/Delta) * int_0^min(Delta,w) S_psi(mode + u - d_pre) du
# This unifies Delta <= w (certain detection) and Delta > w (P(det) = w/Delta).

S_testing_df <- expand_grid(psi=psi_vals, Delta=Delta_vals) %>%
	mutate(mode=pmax(0, ((alpha*psi-1)/beta))) %>%
	rowwise() %>%
	mutate(S=integrate(function(u)
		(1-pgamma(u - d_pre + mode, alpha*psi, beta)),
		lower=0, upper=min(Delta, w))$value / Delta) %>%
	ungroup()

fig_S_testing <- S_testing_df %>%
	ggplot(aes(x=Delta, y=S, col=factor(psi))) +
		geom_line(linewidth=0.8, alpha=0.8) +
		theme_classic() +
		labs(x="Gap between tests (days)", y="Test effectiveness", col="psi")

# --- Perfect sensitivity, exponential delay to action ---
# Now isolation occurs at tau_det + delta_act, where delta_act ~ Exp(lambda_act).
# The inner expectation E_delta[S_psi(t0 + delta)] over delta ~ Exp(lambda) has
# a closed form (by Fubini + Gamma-exponential convolution):
#   I(t0) = [1 - pgamma(t0, a, beta)]
#          - exp(lambda*t0) * (beta/(beta+lambda))^a * [1 - pgamma(t0, a, beta+lambda)]
# This reduces the double integral to a single well-behaved integral over u.

lambda_act <- 1  # rate parameter; mean delay = 1/lambda_act = 1 day

S_testing_delay_df <- expand_grid(psi=psi_vals, Delta=Delta_vals) %>%
	mutate(mode=pmax(0, ((alpha*psi-1)/beta))) %>%
	rowwise() %>%
	mutate(S={
		a <- alpha*psi
		cr <- (beta/(beta+lambda_act))^a
		integrate(function(u) {
			t0 <- u - d_pre + mode
			(1-pgamma(t0, a, beta)) -
			exp(lambda_act*t0) * cr * (1-pgamma(t0, a, beta+lambda_act))
		}, lower=0, upper=min(Delta, w))$value / Delta
	}) %>%
	ungroup()

fig_S_testing_delay <- S_testing_delay_df %>%
	ggplot(aes(x=Delta, y=S, col=factor(psi))) +
		geom_line(linewidth=0.8, alpha=0.8) +
		theme_classic() +
		labs(x="Gap between tests (days)", y="Test effectiveness", col="psi")

# ==============================================================================
# How isolation impacts the secondary infection distribution and
# generation interval distribution
# ==============================================================================

# Based on theory we've developed previously (a sketch is at the end if notes/findings_controlled.md), detect-and-isolate protocols can impact the generation interval distribution. The distribution shrinks, and does so more for smooth infectiousness profiles than for spiky ones (in the spiky limit, \psi -> 0, there is no generation interval distortion). As a result, TE tells only part of the story: if we were to simulate an epidemic with R_{eff} = R_0 (1-TE), we would get a different-looking epidemic than the one we actually see with the detect-and-isolate program. It should be the case that the detect-and-isolate epidemic should grow faster and shrink faster than the simple ajdusted R_{eff} epidemic, since the shorter generation interval distribution should lead to faster dynamics overall. Additinoally: test-and-isolate protocols should lead to more overdispersion, especially in the spike case, where test-and-isolate zeroes out some people's infectiousness but completely misses other people's. Again, this should have an impact on epidemic dynamics: relative to an epidemic with R_{eff} = R_0 (1-TE), the controlled epidemics should be less likely to take off (more overdispersion causes this), but also more explosive when they do. I'd like to generate simulations using the simulation function employed in episims_gamma.R  to compare epidemics that undergo detect-and-isolate from analogous epidemics where only R_{eff} is changed to match what we'd expect from TE. I'd like to show that these epidemics meaningfully differ.

pars <- parslist[[2]]  # omicron
T <- pars$Tgen
alpha <- pars$alpha
beta <- pars$beta
R0 <- pars$R0

psi_vals <- c(0, 0.5, 1)
tau_offset <- -1      # isolation days after peak infectiousness
p_adhere <- 0.5      # probability of adhering to isolation
popsize <- 1000
nsim <- 200

# Compute exact TE for each psi: TE = p_adhere * S_psi(tau_offset)
te_df <- tibble(psi=psi_vals) %>%
	mutate(mode=pmax(0, (alpha*psi-1)/beta),
	       TE=p_adhere * (1-pgamma(tau_offset+mode, alpha*psi, beta)))

cat("Expected TE by psi:\n")
print(te_df)

# Run simulations: (1) fixed isolation vs (2) naive R0*(1-TE) adjustment
results <- vector("list", nsim * length(psi_vals) * 2)
idx <- 1L

for(sim in 1:nsim){
	for(i in seq_along(psi_vals)){
		psi <- psi_vals[i]
		te <- te_df$TE[i]

		# (1) Fixed isolation at tau_offset relative to peak
		tinf <- sim_stochastic_fast(n=popsize,
			gen_inf_attempts=gen_inf_attempts_gamma_fixed_iso(
				T, R0, alpha, psi,
				tau_offset=tau_offset,
				p_adhere=p_adhere))
		infected <- sort(tinf[tinf < Inf])
		results[[idx]] <- tibble(tinf=infected, cuminf=seq_along(infected),
		                         sim=sim, psi=psi, type="isolation")
		idx <- idx + 1L

		# (2) Naive TE adjustment (same GI, reduced R0)
		tinf <- sim_stochastic_fast(n=popsize,
			gen_inf_attempts=gen_inf_attempts_gamma(
				T, R0*(1-te), alpha, psi))
		infected <- sort(tinf[tinf < Inf])
		results[[idx]] <- tibble(tinf=infected, cuminf=seq_along(infected),
		                         sim=sim, psi=psi, type="reduced_R0")
		idx <- idx + 1L
	}
	if(sim %% 50 == 0) cat(sprintf("sim %d/%d\n", sim, nsim))
}

cuminf_iso_df <- bind_rows(results)

# Flag established epidemics (>=10% of population infected)
cuminf_iso_df <- cuminf_iso_df %>%
	group_by(sim, psi, type) %>%
	mutate(established = max(cuminf) >= 0.1 * popsize) %>%
	ungroup()

# Check establishment counts and warn if low
est_counts <- cuminf_iso_df %>%
	group_by(sim, psi, type) %>%
	summarise(fs = max(cuminf), .groups="drop") %>%
	mutate(established = fs >= 0.1 * popsize) %>%
	group_by(psi, type) %>%
	summarise(n_est = sum(established), .groups="drop")

for(i in seq_len(nrow(est_counts))){
	row <- est_counts[i,]
	if(row$n_est == 0){
		cat(sprintf("WARNING: No epidemics established for psi=%s, type=%s. Skipping plots for this group.\n",
		            row$psi, row$type))
	} else if(row$n_est < 30){
		cat(sprintf("WARNING: Only %d epidemics established for psi=%s, type=%s; aggregate statistics may be unreliable.\n",
		            row$n_est, row$psi, row$type))
	}
}

# --- Cumulative incidence overlay ---
fig_cuminf_iso <- cuminf_iso_df %>%
	filter(established) %>%
	ggplot(aes(x=tinf, y=cuminf, group=interaction(sim, type), col=type)) +
		geom_line(alpha=0.15, linewidth=0.3) +
		facet_wrap(~psi, nrow=1, labeller=label_both) +
		theme_classic() +
		labs(x="Time (days)", y="Cumulative infections", col="Strategy") +
		scale_color_manual(values=c("isolation"="steelblue", "reduced_R0"="tomato"))

# --- P(establishment) comparison ---
pest_iso_table <- cuminf_iso_df %>%
	group_by(sim, psi, type) %>%
	summarise(fs = max(cuminf), .groups="drop") %>%
	group_by(psi, type) %>%
	summarise(p_est = mean(fs >= 0.1 * popsize),
	          mean_fs = mean(fs[fs >= 0.1 * popsize]),
	          .groups="drop")

cat("\nP(establishment) and mean final size (conditional on establishment):\n")
print(pest_iso_table %>% pivot_wider(names_from=type, values_from=c(p_est, mean_fs)))

# --- Early growth rate comparison ---
growth_threshold <- 100
min_threshold <- 10

growthrate_iso_df <- cuminf_iso_df %>%
	filter(established, cuminf >= min_threshold, cuminf <= growth_threshold) %>%
	group_by(sim, psi, type) %>%
	filter(n() >= 2) %>%
	summarise(growthrate = coef(lm(log(cuminf) ~ tinf))[2], .groups="drop")

if(nrow(growthrate_iso_df) == 0){
	cat("WARNING: No established epidemics with enough data points for growth rate estimation.\n")
	fig_growthrate_iso <- NULL
} else {
	fig_growthrate_iso <- growthrate_iso_df %>%
		ggplot(aes(x=growthrate, fill=type)) +
			geom_density(alpha=0.4) +
			facet_wrap(~psi, nrow=1, labeller=label_both) +
			theme_classic() +
			labs(x="Early growth rate (per day)", y="Density", fill="Strategy") +
			scale_fill_manual(values=c("isolation"="steelblue", "reduced_R0"="tomato"))
}

# To add: 
# number of established epidemics in the output table 
# Mean growth rate in the output table(?) 


# This script investigates how detect-and-isolate interventions interact with the punctuatedness (ψ) of the
# individual infectiousness profile, through three increasingly realistic isolation models:

  # 1. Fixed isolation time (lines 7–27): Plots TE as a function of isolation timing relative to peak
  # infectiousness, for several ψ values. Shows that for a given isolation time, spikier profiles (ψ→0)
  #  have higher TE when isolation is before the peak (all-or-nothing: the spike hasn't fired yet) but
  # lower TE after the peak (the spike already fired, nothing left to avert). The curves cross near the
  #  mode.
  # 2. Symptom-triggered isolation (lines 29–54): Same idea but the isolation time is random
  # (Normal-distributed symptom onset). Smooths out the sharp step-function behavior from section 1,
  # but the qualitative ordering across ψ values is preserved.
  # 3. Regular screening (lines 56–119): TE as a function of the gap between tests (Δ), with a
  # detectability window around peak infectiousness. Computed both without and with an exponential
  # delay to action (the latter using a closed-form Gamma-exponential convolution to avoid nested
  # quadrature). Shows how testing frequency trades off against profile shape.
  # 4. Epidemic simulations (lines 121–end): The main punchline. Compares epidemics under fixed
  # isolation (using gen_inf_attempts_gamma_fixed_iso) against epidemics that simply reduce R0 to match
  #  the expected TE. The point is that TE only captures the mean reduction in transmission — it misses
  #  two things that isolation actually does:
  #   - GI distortion: Isolation preferentially removes late transmission attempts, shortening the
  # effective generation interval. This makes controlled epidemics grow (and decline) faster than a
  # naive R_eff adjustment would predict.
  #   - Overdispersion: Isolation creates correlated thinning — especially for spike profiles, where
  # it's all-or-nothing (either detected before the spike and fully averted, or not and fully
  # transmitted). This increases variance in individual reproductive numbers beyond what Poisson(R_eff)
  #  would give, reducing establishment probability but making established epidemics more explosive.

  # What it finds (so far):

  # The simulation results show these effects clearly at ψ = 0 vs ψ = 1. For the spike case, isolation
  # dramatically reduces establishment probability relative to the reduced-R0 baseline (0.157 vs 0.624
  # in the omicron/symptom parameterization), confirming that the overdispersion from all-or-nothing
  # isolation is epidemiologically consequential — TE alone substantially overpredicts the threat. For
  # smooth profiles, the gap is smaller because the thinning is more graded. The deterministic-offset
  # simulations are still being explored, but the framework is in place to systematically vary
  # tau_offset and p_adhere to map out where ψ-dependent effects are largest.