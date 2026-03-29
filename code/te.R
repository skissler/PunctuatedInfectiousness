library(tidyverse)

alpha <- 10 
beta <- 2 

dvals <- seq(from=-5, to=5, by=0.1)
psivals <- c(0.01, 0.1, 0.25, 0.5, 0.75, 1)

gam_df <- expand_grid(d=dvals, psi=psivals) %>% 
	mutate(TE=1-pgamma(d + pmax(0,((alpha*psi-1)/beta)), alpha*psi, beta)) %>%
	mutate(psi=factor(psi))

gam_df %>% 
	ggplot(aes(x=d, y=TE, col=psi)) + 
		geom_line() + 
		theme_classic()

# Lognormal analogue (moment-matched: same mean and CV^2 as Gamma)
ln_df <- expand_grid(d = dvals, psi = psivals) %>%
	mutate(
		ln_sig2 = log(1 + 1 / (alpha * psi)),
		ln_mu = log(psi * alpha / beta) - ln_sig2 / 2,
		ln_mode = exp(ln_mu - ln_sig2),
		TE = 1 - plnorm(d + ln_mode, meanlog = ln_mu, sdlog = sqrt(ln_sig2))
	) %>%
	select(-ln_sig2, -ln_mu, -ln_mode) %>%
	mutate(psi = factor(psi))

ln_df %>%
	ggplot(aes(x = d, y = TE, col = psi)) +
		geom_line() +
		theme_classic()

# Random detection time: D ~ N(d, sigma_det^2)
sigma_det <- 1

gam_rand_df <- expand_grid(d = dvals, psi = psivals) %>%
	rowwise() %>%
	mutate(TE = integrate(function(t)
		(1 - pgamma(t + pmax(0, (alpha * psi - 1) / beta), alpha * psi, beta)) *
			dnorm(t, d, sigma_det),
		lower = d - 5 * sigma_det, upper = d + 5 * sigma_det)$value
	) %>%
	ungroup() %>%
	mutate(psi = factor(psi))

gam_rand_df %>%
	ggplot(aes(x = d, y = TE, col = psi)) +
		geom_line() +
		theme_classic()

# Lognormal with random detection
ln_rand_df <- expand_grid(d = dvals, psi = psivals) %>%
	mutate(
		ln_sig2 = log(1 + 1 / (alpha * psi)),
		ln_mu = log(psi * alpha / beta) - ln_sig2 / 2,
		ln_mode = exp(ln_mu - ln_sig2)
	) %>%
	rowwise() %>%
	mutate(TE = integrate(function(t)
		(1 - plnorm(t + ln_mode, meanlog = ln_mu, sdlog = sqrt(ln_sig2))) *
			dnorm(t, d, sigma_det),
		lower = d - 5 * sigma_det, upper = d + 5 * sigma_det)$value
	) %>%
	ungroup() %>%
	select(-ln_sig2, -ln_mu, -ln_mode) %>%
	mutate(psi = factor(psi))

ln_rand_df %>%
	ggplot(aes(x = d, y = TE, col = psi)) +
		geom_line() +
		theme_classic()