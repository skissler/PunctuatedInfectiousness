library(tidyverse)

# ==============================================================================
# Distribution of W as a function of kappa for different R0 values
#
# Two approximations to W | survival:
#   (a) 2-moment Gamma: matches mean and variance (closed-form, no fitting)
#   (b) 3-moment Generalized Gamma: matches mean, variance, and skewness
#       (requires numerical optimisation for shape parameters d, p)
#
# Plots mean, median, mode, and middle 50%/90% intervals for each.
# ==============================================================================

# ==============================================================================
# 1. Parameters
# ==============================================================================

T_gen    <- 5          # Mean generation interval
popshape <- 10         # Shape parameter (alpha_total)
r        <- popshape / T_gen
R0_vals  <- c(2, 5, 12)
kappa_vals <- seq(0.01, 0.99, length.out = 300)

# ==============================================================================
# 2. Helper functions
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

#' Fit generalized gamma parameters (a, d, p) to three moments.
#' GenGamma moments: E[X^k] = a^k * Gamma((d+k)/p) / Gamma(d/p)
#' Returns list(a, d, p).
fit_gengamma <- function(EW, EW2, EW3, start_d = 1, start_p = 1) {
  target_lr2 <- log(EW2 / EW^2)
  target_lr3 <- log(EW3 / EW^3)

  obj <- function(par) {
    d <- exp(par[1]); p <- exp(par[2])
    if (d < 1e-6 || p < 1e-6 || d > 100 || p > 100) return(1e10)
    lr2 <- lgamma((d + 2) / p) + lgamma(d / p) - 2 * lgamma((d + 1) / p)
    lr3 <- lgamma((d + 3) / p) + 2 * lgamma(d / p) - 3 * lgamma((d + 1) / p)
    if (is.nan(lr2) || is.nan(lr3) ||
        is.infinite(lr2) || is.infinite(lr3)) return(1e10)
    (lr2 - target_lr2)^2 + (lr3 - target_lr3)^2
  }

  # Try warm start first, then grid if needed
  best <- tryCatch(
    optim(c(log(start_d), log(start_p)), obj, method = "Nelder-Mead",
          control = list(maxit = 5000, reltol = 1e-12)),
    error = function(e) list(value = Inf)
  )

  if (best$value > 1e-10) {
    for (d0 in c(0.5, 1, 2, 5, 10)) {
      for (p0 in c(0.5, 1, 2, 3, 5)) {
        fit <- tryCatch(
          optim(c(log(d0), log(p0)), obj, method = "Nelder-Mead",
                control = list(maxit = 5000, reltol = 1e-12)),
          error = function(e) list(value = Inf)
        )
        if (fit$value < best$value) best <- fit
      }
    }
  }

  d <- exp(best$par[1]); p <- exp(best$par[2])
  a <- EW * exp(lgamma(d / p) - lgamma((d + 1) / p))
  list(a = a, d = d, p = p)
}

#' Quantile function for generalized gamma.
#' CDF: F(x) = pgamma((x/a)^p, d/p, 1)
qgengamma <- function(prob, a, d, p) {
  a * qgamma(prob, shape = d / p, rate = 1)^(1 / p)
}

# ==============================================================================
# 3. Compute exact moments of W for each (kappa, R0)
# ==============================================================================

df <- expand_grid(kappa = kappa_vals, R0 = R0_vals) %>%
  mutate(
    # Malthusian growth rate
    alpha = r * (R0^(1/popshape) - 1),

    # Discount ratios
    rho  = r / (r + alpha),
    rho2 = r / (r + 2 * alpha),
    rho3 = r / (r + 3 * alpha),

    # m2, m3 (kappa-independent)
    m2 = R0 * rho2^popshape,
    m3 = R0 * rho3^popshape,

    # Shape parameters
    kappa_shape = kappa * popshape,
    shift_shape = (1 - kappa) * popshape,

    # --- First two moments (exact, closed-form) ---
    var_S = m2 + R0^2 * rho^(2 * kappa_shape) *
      (rho2^shift_shape - rho^(2 * shift_shape)),
    var_W = var_S / (1 - m2),
    EW2 = 1 + var_W,  # since E[W] = 1

    # --- Third moment (exact, closed-form) ---
    # Poisson factorial moments: E[N(N-1)] = R0^2, E[N(N-1)(N-2)] = R0^3
    m21 = R0^2 * rho2^kappa_shape * rho^kappa_shape * rho3^shift_shape,
    m111 = R0^3 * rho^(3 * kappa_shape) * rho3^shift_shape,
    EW3 = (3 * m21 * EW2 + m111) / (1 - m3),

    # Extinction probability
    q = map_dbl(R0, extinction_prob),

    # --- Conditional moments (W | survival) ---
    E_cond  = 1 / (1 - q),
    E2_cond = EW2 / (1 - q),
    E3_cond = EW3 / (1 - q),
    var_cond = E2_cond - E_cond^2,

    # =========================================================================
    # (a) 2-moment Gamma approximation
    # =========================================================================
    g2_rate  = E_cond / var_cond,
    g2_shape = E_cond * g2_rate,
    g2_q05  = qgamma(0.05, shape = g2_shape, rate = g2_rate),
    g2_q25  = qgamma(0.25, shape = g2_shape, rate = g2_rate),
    g2_q50  = qgamma(0.50, shape = g2_shape, rate = g2_rate),
    g2_q75  = qgamma(0.75, shape = g2_shape, rate = g2_rate),
    g2_q95  = qgamma(0.95, shape = g2_shape, rate = g2_rate),
    g2_mode = ifelse(g2_shape >= 1, (g2_shape - 1) / g2_rate, 0),

    R0_label = paste0("R[0] == ", R0)
  )

# ==============================================================================
# 4. Fit 3-moment generalized gamma (requires loop for warm-starting)
# ==============================================================================

cat("Fitting generalized gamma (3-moment) ...\n")

# Pre-allocate columns
df$g3_a <- NA_real_; df$g3_d <- NA_real_; df$g3_p <- NA_real_

for (R0_val in R0_vals) {
  idx <- which(df$R0 == R0_val)
  prev_d <- 1; prev_p <- 1

  for (i in idx) {
    fit <- fit_gengamma(df$E_cond[i], df$E2_cond[i], df$E3_cond[i],
                        start_d = prev_d, start_p = prev_p)
    df$g3_a[i] <- fit$a; df$g3_d[i] <- fit$d; df$g3_p[i] <- fit$p
    prev_d <- fit$d; prev_p <- fit$p
  }
}

df <- df %>%
  mutate(
    g3_q05  = qgengamma(0.05, g3_a, g3_d, g3_p),
    g3_q25  = qgengamma(0.25, g3_a, g3_d, g3_p),
    g3_q50  = qgengamma(0.50, g3_a, g3_d, g3_p),
    g3_q75  = qgengamma(0.75, g3_a, g3_d, g3_p),
    g3_q95  = qgengamma(0.95, g3_a, g3_d, g3_p),
    g3_mode = ifelse(g3_d > 1, g3_a * ((g3_d - 1) / g3_p)^(1 / g3_p), 0)
  )

cat("Done.\n")

# ==============================================================================
# 5. Plot: 2-moment Gamma
# ==============================================================================

p_gamma <- ggplot(df, aes(x = kappa)) +
  geom_ribbon(aes(ymin = g2_q05, ymax = g2_q95),
              fill = "steelblue", alpha = 0.15) +
  geom_ribbon(aes(ymin = g2_q25, ymax = g2_q75),
              fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = g2_q50, colour = "Median"), linewidth = 0.7) +
  geom_line(aes(y = E_cond, colour = "Mean"), linewidth = 0.7) +
  geom_line(aes(y = g2_mode, colour = "Mode"), linewidth = 0.7, linetype = "dashed") +
  scale_colour_manual(
    values = c("Mean" = "black", "Median" = "firebrick", "Mode" = "grey40"),
    name = NULL
  ) +
  facet_wrap(~ R0_label, scales = "free_y", labeller = label_parsed) +
  labs(
    x = expression(kappa),
    y = expression(italic(W) ~ "|" ~ survival),
    title = expression("Distribution of " * italic(W) *
                       " conditional on survival (2-moment Gamma)"),
    subtitle = expression(paste(
      "Middle 50% (dark) and 90% (light); ",
      alpha[total], " = 10, T = 5"))
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(size = 13),
    plot.title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# ==============================================================================
# 6. Plot: 3-moment Generalized Gamma
# ==============================================================================

p_gengamma <- ggplot(df, aes(x = kappa)) +
  geom_ribbon(aes(ymin = g3_q05, ymax = g3_q95),
              fill = "steelblue", alpha = 0.15) +
  geom_ribbon(aes(ymin = g3_q25, ymax = g3_q75),
              fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = g3_q50, colour = "Median"), linewidth = 0.7) +
  geom_line(aes(y = E_cond, colour = "Mean"), linewidth = 0.7) +
  geom_line(aes(y = g3_mode, colour = "Mode"), linewidth = 0.7, linetype = "dashed") +
  scale_colour_manual(
    values = c("Mean" = "black", "Median" = "firebrick", "Mode" = "grey40"),
    name = NULL
  ) +
  facet_wrap(~ R0_label, scales = "free_y", labeller = label_parsed) +
  labs(
    x = expression(kappa),
    y = expression(italic(W) ~ "|" ~ survival),
    title = expression("Distribution of " * italic(W) *
                       " conditional on survival (3-moment Gen. Gamma)"),
    subtitle = expression(paste(
      "Middle 50% (dark) and 90% (light); ",
      alpha[total], " = 10, T = 5"))
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(size = 13),
    plot.title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(p_gamma)
print(p_gengamma)

ggsave("figures/fig_W_variance_bands_gamma.pdf", p_gamma, width = 10, height = 4.5)
ggsave("figures/fig_W_variance_bands_gengamma.pdf", p_gengamma, width = 10, height = 4.5)
