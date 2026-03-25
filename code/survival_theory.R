library(tidyverse)

# ==============================================================================
# Theoretical survival curves via deterministic ODE + time shift
#
# Assumes code/1_1_basicmetrics.R has been run first. Requires from environment:
#   survival_100_df, survival_peak_df, ode_out, ode_daily,
#   alpha_theoretical, breakvals, popsize, e_dur, i_dur, R0
#
# Overlays theoretical curves for spike (kappa->0) and smooth (kappa->inf)
# limits using the Morris et al. (2024) time-shift approximation:
#   run the deterministic ODE once, shift in time by tau = log(W)/alpha,
#   where W is drawn from an analytically-derived distribution.
# ==============================================================================

# ==============================================================================
# 1. Helper function
# ==============================================================================

extinction_prob <- function(R0) {
  if (R0 <= 1) return(1)
  q <- 0.5
  for (i in 1:1000) {
    q_new <- exp(R0 * (q - 1))
    if (abs(q_new - q) < 1e-12) break
    q <- q_new
  }
  q
}

# ==============================================================================
# 2. Analytical W moments for spike and smooth limits
# ==============================================================================

alpha <- alpha_theoretical  # equals 1 for these parameters

# Discount factors for the hypoexponential generation interval
rho2 <- 1 / ((1 + 2 * alpha * e_dur) * (1 + 2 * alpha * i_dur))
rho3 <- 1 / ((1 + 3 * alpha * e_dur) * (1 + 3 * alpha * i_dur))

m2 <- R0 * rho2
m3 <- R0 * rho3

q <- extinction_prob(R0)

# Smooth limit (kappa -> inf): all generation times independent
var_S_smooth  <- m2

# Spike limit (kappa -> 0): all offspring share same generation time
var_S_spike   <- m2 + R0^2 * rho2 - 1

# Unconditional W moments (E[W] = 1 always)
var_W_smooth <- var_S_smooth / (1 - m2)
var_W_spike  <- var_S_spike  / (1 - m2)

# Conditional moments (W | survival, i.e. W > 0)
E_cond <- 1 / (1 - q)
var_cond_smooth <- (1 + var_W_smooth) / (1 - q) - E_cond^2
var_cond_spike  <- (1 + var_W_spike)  / (1 - q) - E_cond^2

cat("alpha =", alpha, "\n")
cat("q     =", q, "\n")
cat("smooth: var_W =", var_W_smooth, ", cond_var =", var_cond_smooth, "\n")
cat("spike:  var_W =", var_W_spike,  ", cond_var =", var_cond_spike, "\n")

# ==============================================================================
# 3. Fit 2-moment Gamma to W|surv and sample W
# ==============================================================================

n_draws <- 100000
set.seed(42)

gamma_params <- list(
  smooth = list(
    shape = E_cond^2 / var_cond_smooth,
    rate  = E_cond   / var_cond_smooth
  ),
  spike = list(
    shape = E_cond^2 / var_cond_spike,
    rate  = E_cond   / var_cond_spike
  )
)

theory_draws <- list()
for (nm in names(gamma_params)) {
  gp <- gamma_params[[nm]]
  cat(nm, ": Gamma shape =", gp$shape, ", rate =", gp$rate, "\n")

  W   <- rgamma(n_draws, shape = gp$shape, rate = gp$rate)
  tau <- log(W) / alpha

  theory_draws[[nm]] <- tibble(profiletype = nm, W = W, tau = tau)
}

theory_draws_df <- bind_rows(theory_draws)

# ==============================================================================
# 4. Deterministic milestone times from ODE
# ==============================================================================

# t_det_100: continuous time when cumulative infections cross 100
t_det_100 <- approx(
  x = ode_out$cuminf * popsize,
  y = ode_out$t,
  xout = 100
)$y

cat("Deterministic time to 100 cases:", t_det_100, "\n")

# t_det_peak: day with maximum daily new infections
t_det_peak <- ode_daily$day[which.max(ode_daily$newinf)]
cat("Deterministic time to peak:", t_det_peak, "\n")

# ==============================================================================
# 5. Build theoretical survival curves
# ==============================================================================

# Stochastic milestone times = deterministic time - tau
theory_draws_df <- theory_draws_df %>%
  mutate(
    t_100  = t_det_100  - tau,
    t_peak = t_det_peak - tau
  )

# Survival curve for time to 100 cases
t_grid_100 <- seq(0, max(survival_100_df$t), length.out = 500)

theory_surv_100 <- theory_draws_df %>%
  group_by(profiletype) %>%
  summarise(
    t = list(t_grid_100),
    prop_below = list(sapply(t_grid_100, function(tt) mean(t_100 > tt))),
    .groups = "drop"
  ) %>%
  unnest(cols = c(t, prop_below))

# Survival curve for time to peak
t_grid_peak <- seq(0, max(survival_peak_df$t), length.out = 500)

theory_surv_peak <- theory_draws_df %>%
  group_by(profiletype) %>%
  summarise(
    t = list(t_grid_peak),
    prop_below = list(sapply(t_grid_peak, function(tt) mean(t_peak > tt))),
    .groups = "drop"
  ) %>%
  unnest(cols = c(t, prop_below))

# ==============================================================================
# 6. Overlay figures
# ==============================================================================

profile_colors <- c("stepwise" = "black", "smooth" = "blue", "spike" = "red")

# --- Time to 100 cases ---
fig_survival_100_theory <- ggplot() +
  geom_line(data = survival_100_df,
            aes(x = t, y = prop_below, col = profiletype),
            alpha = 0.6, linewidth = 0.8) +
  geom_line(data = theory_surv_100,
            aes(x = t, y = prop_below, col = profiletype),
            linetype = "dashed", linewidth = 0.8) +
  scale_x_continuous(breaks = breakvals) +
  scale_color_manual(values = profile_colors) +
  labs(x = "Days",
       y = paste0("Proportion not yet reaching 100 cases"),
       col = "Profile") +
  theme_classic()

# --- Time to peak ---
fig_survival_peak_theory <- ggplot() +
  geom_line(data = survival_peak_df,
            aes(x = t, y = prop_below, col = profiletype),
            alpha = 0.6, linewidth = 0.8) +
  geom_line(data = theory_surv_peak,
            aes(x = t, y = prop_below, col = profiletype),
            linetype = "dashed", linewidth = 0.8) +
  scale_x_continuous(breaks = breakvals) +
  scale_color_manual(values = profile_colors) +
  labs(x = "Days",
       y = "Proportion not yet reaching peak daily incidence",
       col = "Profile") +
  theme_classic()

print(fig_survival_100_theory)
print(fig_survival_peak_theory)

ggsave("figures/fig_survival_100_theory.pdf",
       fig_survival_100_theory, width = 8, height = 5)
ggsave("figures/fig_survival_peak_theory.pdf",
       fig_survival_peak_theory, width = 8, height = 5)

cat("Saved figures/fig_survival_100_theory.pdf\n")
cat("Saved figures/fig_survival_peak_theory.pdf\n")
