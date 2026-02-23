library(tidyverse) 
library(odin)

source("code/utils.R")

# set simulation parameters: 
popsize <- 1000
e_dur <- 2  # 0
i_dur <- 3  # 5 
R0 <- 2 # 1.5 

# run the mean-field ode: 
if(e_dur > 0){
	sir_model <- seir$new(R0=R0, nu=1/e_dur, gamma=1/i_dur, N=popsize)
	t <- seq(0, 100, length.out = 1000)
	sir_out <- as_tibble(data.frame(sir_model$run(t)))
	print("Used SEIR for mean-field")
} else {
	sir_model <- sir$new(R0=R0, gamma=1/i_dur, N=popsize)
	t <- seq(0, 100, length.out = 1000)
	sir_out <- as_tibble(data.frame(sir_model$run(t)))
	print("Used SIR for mean-field")
}

# Aggregate mean-field output into daily counts: 
sir_out_daily <- sir_out %>% 
	mutate(day=floor(t)) %>% 
	group_by(day) %>% 
	summarise(cuminf=max(cuminf)) %>% 
	arrange(day) %>% 
	mutate(newinf=cuminf-lag(cuminf)) %>% 
	mutate(newinf=case_when(is.na(newinf)~cuminf, TRUE~newinf))

sir_out_weekly <- sir_out %>% 
	mutate(day=floor(t)) %>% 
	mutate(week=floor(day/7)) %>% 
	group_by(week) %>% 
	summarise(cuminf=max(cuminf)) %>% 
	arrange(week) %>% 
	mutate(newinf=cuminf-lag(cuminf)) %>% 
	mutate(newinf=case_when(is.na(newinf)~cuminf, TRUE~newinf))

# run the stochastic simulations: 
cuminf_df <- tibble()
for(sim in 1:1000){
	tinf <- sim_stochastic(n=popsize, e_dur=e_dur, i_dur=i_dur, R0=R0, profiletype="stepwise")	
	cuminf_df <- bind_rows(cuminf_df, 
		(tibble(tinf=sort(tinf[tinf<Inf])) %>% 
				mutate(cuminf=1:n()) %>% 
				mutate(sim=sim) %>% 
				mutate(profiletype="stepwise")))

	tinf <- sim_stochastic(n=popsize, e_dur=e_dur, i_dur=i_dur, R0=R0, profiletype="smooth")	
	cuminf_df <- bind_rows(cuminf_df, 
		(tibble(tinf=sort(tinf[tinf<Inf])) %>% 
				mutate(cuminf=1:n()) %>% 
				mutate(sim=sim) %>% 
				mutate(profiletype="smooth")))

	tinf <- sim_stochastic(n=popsize, e_dur=e_dur, i_dur=i_dur, R0=R0, profiletype="spike")	
	cuminf_df <- bind_rows(cuminf_df, 
		(tibble(tinf=sort(tinf[tinf<Inf])) %>% 
				mutate(cuminf=1:n()) %>% 
				mutate(sim=sim) %>% 
				mutate(profiletype="spike")))
	if(sim%%100 == 0){
		print(sim)
	}
}
cuminf_df <- cuminf_df %>% 
	mutate(profiletype=factor(profiletype, levels=c("stepwise","smooth","spike"))) %>% 
	group_by(sim, profiletype) %>% 
	mutate(ninf=n()) %>% 
	mutate(established=case_when(ninf>=0.1*popsize~1, TRUE~0)) %>% 
	select(-ninf) %>% 
	ungroup() 

# Overlay the daily infection curves: 
dayjoin <- expand_grid(
	profiletype=unique(cuminf_df$profiletype),
	sim=(1:max(cuminf_df$sim)),
	day=(0:max(ceiling(cuminf_df$tinf))))

weekjoin <- expand_grid(
	profiletype=unique(cuminf_df$profiletype),
	sim=(1:max(cuminf_df$sim)),
	week=(0:max(ceiling(cuminf_df$tinf/7))))

dailyinf_df <- cuminf_df %>% 
	mutate(day=floor(tinf)) %>% 
	group_by(profiletype, sim, day) %>% 
	summarise(ninf=n(), established=first(established)) %>% 
	right_join(dayjoin, by=c("profiletype","sim","day")) %>% 
	replace_na(list(ninf=0)) %>% 
	group_by(profiletype, sim) %>% 
	mutate(established=max(established, na.rm=TRUE))

weeklyinf_df <- cuminf_df %>% 
	mutate(day=floor(tinf)) %>% 
	mutate(week=floor(day/7)) %>% 
	group_by(profiletype, sim, week) %>% 
	summarise(ninf=n(), established=first(established)) %>% 
	right_join(weekjoin, by=c("profiletype","sim","week")) %>% 
	replace_na(list(ninf=0)) %>% 
	group_by(profiletype, sim) %>% 
	mutate(established=max(established, na.rm=TRUE))

# Overlay the cumulative infection curves: 
fig_cuminf_overlay <- cuminf_df %>% 
	ggplot(aes(x=tinf, y=cuminf, group=sim)) + 
		geom_line(alpha=0.2, col="grey") + 
		geom_line(data=(cuminf_df %>% 
			filter(established==1) %>% 
			group_by(profiletype) %>% 
			filter(sim==min(sim))), linewidth=0.2, alpha=0.8, col="black") + 
		geom_line(data=sir_out_daily, aes(x=day, y=cuminf*1000), alpha=0.6, linewidth=1, col="blue") + 
		theme_classic() + 
		facet_wrap(~profiletype, nrow=1)

fig_cuminf_means <- dailyinf_df %>% 
	filter(established==1) %>% 
	# filter(profiletype!="stepwise") %>% 
	arrange(profiletype, sim, day) %>% 
	mutate(cumninf=cumsum(ninf)) %>% 
	group_by(profiletype, day) %>% 
	summarise(cumninf_mean=mean(cumninf), cumninf_sd=sd(cumninf)) %>% 
	ggplot(aes(x=day, y=cumninf_mean, col=profiletype)) + 
		geom_line() + 
		geom_ribbon(aes(x=day, ymin=cumninf_mean-cumninf_sd, ymax=cumninf_mean+cumninf_sd, fill=profiletype), alpha=0.2) + 
		geom_line(data=sir_out_daily, aes(x=day, y=cuminf*1000), alpha=0.6, linewidth=1, col="magenta", lty="dashed") + 
		theme_classic() + 
		scale_color_manual(values=c("black","blue","red")) + 
		scale_fill_manual(values=c("black","blue","red")) + 
		facet_wrap(~profiletype)

fig_inf_overlay <- dailyinf_df %>% 
	ggplot(aes(x=day, y=ninf, group=sim)) + 
		geom_line(alpha=0.2, col="grey") + 
		geom_line(data=(dailyinf_df %>% 
			filter(established==1) %>% 
			group_by(profiletype) %>% 
			filter(sim==min(sim))), linewidth=0.2, alpha=0.8, col="black") + 
		geom_line(data=sir_out_daily, aes(x=day, y=newinf*1000), alpha=0.6, linewidth=1, col="blue") + 
		theme_classic()  + 
		facet_wrap(~profiletype, nrow=1)

fig_inf_wk_overlay <- weeklyinf_df %>% 
	ggplot(aes(x=week, y=ninf, group=sim)) + 
		geom_line(alpha=0.2, col="grey") + 
		geom_line(data=(weeklyinf_df %>% 
			filter(established==1) %>% 
			group_by(profiletype) %>% 
			filter(sim==min(sim))), linewidth=0.2, alpha=0.8, col="black") + 
		geom_line(data=sir_out_weekly, aes(x=week, y=newinf*1000), alpha=0.6, linewidth=1, col="blue") + 
		theme_classic()  + 
		facet_wrap(~profiletype, nrow=1)

weeklyinf_df %>% 
	ggplot(aes(x=week, y=log(ninf+1), group=sim)) + 
		geom_line(alpha=0.2, col="grey") + 
		geom_line(data=(weeklyinf_df %>% 
			filter(established==1) %>% 
			group_by(profiletype) %>% 
			filter(sim==min(sim))), linewidth=0.2, alpha=0.8, col="black") + 
		geom_line(data=sir_out_weekly, aes(x=week, y=log(newinf*1000+1)), alpha=0.6, linewidth=1, col="blue") + 
		theme_classic()  + 
		facet_wrap(~profiletype, nrow=1)


# Die-out probability? 
pest_table <- cuminf_df %>% 
	group_by(sim, profiletype) %>% 
	slice(1) %>% 
	group_by(profiletype) %>% 
	summarise(pestablished=sum(established)/1000)

# # Week-to-week variation? 
weekvar_table <- weeklyinf_df %>% 
	filter(established==1) %>% 
	arrange(profiletype, sim, week) %>% 
	group_by(profiletype, sim) %>% 
	# mutate(logninfdiff=(ninf)-lag((ninf))) %>% 
	mutate(logninfdiff=log(ninf+1)-lag(log(ninf+1))) %>% 
	group_by(profiletype) %>% 
	summarise(logninfdiff_mean=mean(logninfdiff, na.rm=TRUE), loginfdiff_sd = sd(logninfdiff, na.rm=TRUE))

fs_df <- cuminf_df %>% 
	group_by(profiletype, sim) %>% 
	summarise(finalsize=max(cuminf))

fs_df %>% 
	ggplot(aes(x=finalsize)) + 
		geom_histogram(bins=50, fill="white", col="darkgrey") + 
		theme_classic() + 
		facet_wrap(~profiletype, nrow=1)

fig_fs <- fs_df %>% 
	filter(finalsize>100) %>% 
	ggplot(aes(x=finalsize)) + 
		geom_histogram(aes(y=after_stat(density)), bins=50, fill="white", col="darkgrey") + 
		# geom_density(adjust=2) + 
		theme_classic() + 
		facet_wrap(~profiletype, nrow=1)

fs_table <- fs_df %>% 
	filter(finalsize>100) %>% 
	group_by(profiletype) %>% 
	summarise(fs_mean=mean(finalsize), fs_sd=sd(finalsize))


fig_t100 <- cuminf_df %>% 
	filter(cuminf==100) %>% 
	ggplot(aes(x=tinf)) + 
		geom_histogram(aes(y=after_stat(density)), bins=50, fill="white", col="darkgrey") + 
		geom_density(adjust=2) + 
		theme_classic() + 
		facet_wrap(~profiletype, nrow=1) + 
		labs(title="Time to 100 infections")

t100_table <- cuminf_df %>% 
	filter(cuminf==100) %>% 
	group_by(profiletype) %>% 
	summarise(t100_mean=mean(tinf), t100_sd=sd(tinf))

fig_t100_250 <- cuminf_df %>% 
	filter(cuminf==100 | cuminf==250) %>% 
	pivot_wider(names_from=cuminf, values_from=tinf) %>% 
	mutate(tdiff=`250`-`100`) %>% 
	ggplot(aes(x=tdiff)) + 
		geom_histogram(aes(y=after_stat(density)), bins=50, fill="white", col="darkgrey") + 
		geom_density(adjust=2) + 
		theme_classic() + 
		facet_wrap(~profiletype, nrow=1) + 
		labs(title="Time from 100 to 250 infections")

# fig_t250_750 <- cuminf_df %>% 
# 	filter(cuminf==250 | cuminf==750) %>% 
# 	pivot_wider(names_from=cuminf, values_from=tinf) %>% 
# 	mutate(tdiff=`750`-`250`) %>% 
# 	ggplot(aes(x=tdiff)) + 
# 		geom_histogram(aes(y=after_stat(density)), bins=50, fill="white", col="darkgrey") + 
# 		geom_density(adjust=2) + 
# 		theme_classic() + 
# 		facet_wrap(~profiletype, nrow=1) + 
# 		labs(title="Time from 250 to 750 infections")

# fig_t750_end <- cuminf_df %>% 
# 	filter(cuminf>=750) %>% 
# 	group_by(profiletype, sim) %>% 
# 	filter(cuminf==750 | cuminf==max(cuminf)) %>% 
# 	arrange(profiletype, sim, cuminf) %>% 
# 	group_by(profiletype, sim) %>% 
# 	mutate(tdiff=tinf-lag(tinf)) %>% 
# 	filter(!is.na(tdiff)) %>% 
# 	ggplot(aes(x=tdiff)) + 
# 		geom_histogram(aes(y=after_stat(density)), bins=50, fill="white", col="darkgrey") + 
# 		geom_density(adjust=2) + 
# 		theme_classic() + 
# 		facet_wrap(~profiletype, nrow=1) + 
# 		labs(title="Time from 750 infections to epidemic conclusion")


fig_t_dieout <- cuminf_df %>% 
	group_by(sim, profiletype) %>% 
	summarise(tmax=max(tinf)) %>% 
	ungroup() %>% 
	ggplot(aes(x=tmax)) + 
		geom_histogram(aes(y=after_stat(density)), bins=50, fill="white", col="darkgrey") + 
		# geom_density(adjust=1) + 
		theme_classic() + 
		facet_wrap(~profiletype, nrow=1) + 
		labs(title="Time to epidemic die-out", x="Time (days)")


fig_peak <- dailyinf_df %>% 
	filter(established==1) %>% 
	group_by(sim, profiletype) %>% 
	summarise(peak=max(ninf)) %>% 
	ungroup() %>% 
	ggplot(aes(x=peak)) + 
		geom_histogram(aes(y=after_stat(density)), bins=50, fill="white", col="darkgrey") + 
		# geom_density(adjust=1) + 
		theme_classic() + 
		facet_wrap(~profiletype, nrow=1) + 
		labs(title="Peak daily number of infections", x="Number of infections")


# Additional things to measure: 
# Variability in daily cases 
# can we quantify how forecastable? https://journals.ametsoc.org/view/journals/atsc/61/20/1520-0469_2004_061_2425_paitpi_2.0.co_2.xml

# How much time does it take for the infection to die out? 
# And an idea - how do these dynamics play out under various control scenarios, or various exogeneous forcing scenarios (seasonality, other types of driving)? 

# ------------------------------------------------------------------------------
# I'm wondering if the variation we're seeing is actually because the secondary infection distribution is more variable for the stepwise scenario than for the other two. Let's look at this: 

# beta <- R0/i_dur
# tend <- rexp(10000, 1/i_dur)
# secondaryinfs <- unlist(lapply(beta*tend, function(x){rpois(1, x)}))
# secondaryinfs_poisson <- rpois(10000, R0)






















