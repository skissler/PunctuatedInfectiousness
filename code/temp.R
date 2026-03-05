# ==============================================================================
# Full epidemic simulations
# ==============================================================================

popsize <- 1000
nsim <- 200 

scenarios <- list(
	list(kappa=0.95, z_mean=5, z_amp=0, label="Smooth + constant"),
	list(kappa=0.05, z_mean=5, z_amp=0, label="Punctuated + constant"),
	list(kappa=0.95, z_mean=5, z_amp=4, label="Smooth + periodic"),
	list(kappa=0.05, z_mean=5, z_amp=4, label="Punctuated + periodic")
	)

epi_df <- tibble() 

for(sc in scenarios){

	z <- make_contact_fn(z_mean=sc$z_mean, z_amp=sc$z_amp, z_per=z_per)
	z_max <- sc$z_mean + sc$z_amp

	gfun <- gen_inf_attempts_gamma_contacts(T=T, z=z, z_max=z_max,
	                                        popshape=popshape, kappa=sc$kappa)

	for(sim in seq_len(nsim)){
		tinf <- sim_stochastic_fast(n=popsize, gen_inf_attempts=gfun)
		n_infected <- sum(tinf<Inf)
		epi_df <- bind_rows(epi_df, tibble(
			sim=sim,
			kappa=sc$kappa,
			z_mean=sc$z_mean,
			z_amp=sc$z_amp,
			label=sc$label,
			final_size=n_infected,
			tinf_sorted=list(sort(tinf[tinf<Inf]))
		))
		if(sim %% 20 == 0) cat(sprintf("  %d / %d\n", sim, nsim))
	}

}

epi_df <- epi_df %>% 
	mutate(established=as.integer(final_size>=0.1*popsize))

fig_final_size <- epi_df %>%
	filter(established==1) %>% 
	ggplot(aes(x = final_size)) +
	geom_histogram(binwidth = 1, fill = "steelblue", col = "white", alpha = 0.7) +
	facet_wrap(~label, ncol = 2) +
	theme_classic() +
	labs(x = "Final epidemic size",
	     y = "Count",
	     title = "Final size distributions across scenarios",
	     subtitle = sprintf("N = %d, z_mean = %g, %d simulations per scenario",
	                        popsize, z_mean, nsim))

fs_table <- epi_df %>%
	group_by(label) %>%
	summarise(
		mean_fs = mean(final_size),
		sd_fs   = sd(final_size),
		median_fs = median(final_size),
		p_established = mean(established),
		.groups = "drop"
	)


fs_table_established <- epi_df %>%
	filter(established==1) %>%
	group_by(label) %>%
	summarise(
		n = n(),
		mean_fs = mean(final_size),
		sd_fs   = sd(final_size),
		.groups = "drop"
	) %>%
	print()


curve_df <- tibble()
for (idx in seq_len(nrow(epi_df))) {
	row <- epi_df[idx, ]
	if (!row$established) next
	ts <- row$tinf_sorted[[1]]
	curve_df <- bind_rows(curve_df, tibble(
		sim = row$sim,
		label = row$label,
		tinf = ts,
		cuminf = seq_along(ts)
	))
}

# Keep first 10 established epidemics per scenario for plotting
if (nrow(curve_df) > 0) {
	plot_sims <- curve_df %>%
		group_by(label) %>%
		distinct(sim) #%>%
		#slice_head(n = 10)

	fig_curves <- curve_df %>%
		semi_join(plot_sims, by = c("label", "sim")) %>%
		ggplot(aes(x = tinf, y = cuminf, group = interaction(label, sim))) +
		geom_line(alpha = 0.5) +
		facet_wrap(~label, ncol = 2) +
		theme_classic() +
		labs(x = "Time (days)", y = "Cumulative infections",
		     title = "Example epidemic trajectories (established epidemics)")

	# ggsave("figures/fig_epidemic_curves.pdf", fig_curves, width = 10, height = 7)
	# cat("Saved figures/fig_epidemic_curves.pdf\n")
}

