parslist <- list(
	list(pathogen="influenza",
	     Tgen=3.2,
	     Tvar=2.1^2,
	     R0=2,
	     doi="10.1016/j.epidem.2025.100815"),
	list(pathogen="omicron",
	     Tgen=2.39*2.95,
	     Tvar=2.39*2.95^2,
	     R0=6,
	     doi="10.1016/j.lanepe.2022.100446"),
	list(pathogen="measles",
	     Tgen=12.2,
	     Tvar=3.62^2,
	     R0=12,
	     doi="10.1016/j.sciencedirect.com/S0022519311003146")
)

# Compute moment-matched Gamma shape and rate for each pathogen
# Gamma(shape, rate): mean = shape/rate, var = shape/rate^2
#   => shape = mean^2 / var, rate = mean / var
for (i in seq_along(parslist)) {
	parslist[[i]]$popshape <- parslist[[i]]$Tgen^2 / parslist[[i]]$Tvar
	parslist[[i]]$poprate  <- parslist[[i]]$Tgen   / parslist[[i]]$Tvar
}

# Print summary
for (p in parslist) {
	cat(sprintf("%-12s  Tgen=%.2f  Tvar=%.2f  popshape=%.2f  poprate=%.3f\n",
	    p$pathogen, p$Tgen, p$Tvar, p$popshape, p$poprate))
}

