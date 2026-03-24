library(tidyverse)
library(patchwork)

# ==============================================================================
# Generation interval distortion under isolation — analytical derivation
#
# Isolation truncates individual infectiousness profiles at detection time D,
# preferentially removing late transmissions from smooth profiles. This
# shortens the effective GI and inflates the growth rate relative to what
# R_eff alone would predict.
#
# The tilted Euler-Lotka equation:
#   R₀ · (β/(β+r))^α · [1 − η · pgamma(D, ψα, β+r, lower=FALSE)] = 1
#
# vs the naive equation:
#   R_eff · (β/(β+r))^α = 1
#
# The difference: β+r replaces β inside the survival function (exponential
# tilting). For r > 0, this makes the survival function smaller → fewer
# transmissions "effectively averted" → actual growth rate > naive prediction.
#
# See notes/math_output.html for the full derivation.
# ==============================================================================

# ==============================================================================
# 1. Parameters
# ==============================================================================

T_gen    <- 5
alpha    <- 10           # population-level shape
beta     <- alpha / T_gen  # rate = 2

R0_vals  <- c(2, 5)
eta      <- 1            # full coverage for main analysis
D0_vals  <- c(0.5, 1.0, 1.5, 2.0)

psi_fine <- seq(0.05, 1.0, length.out = 200)
psi_panel_A <- c(0.1, 0.5, 1.0)

# Color palette: red (spike) → blue (smooth), consistent with overlap rule
psi_colors_3 <- setNames(
  c("#C62828", "#E65100", "#0D47A1"),
  as.character(psi_panel_A)
)

D0_colors <- setNames(
  c("#2E7D32", "#E65100", "#C62828", "#7B1FA2"),
  as.character(D0_vals)
)

R0_colors <- setNames(
  c("#1565C0", "#C62828"),
  as.character(R0_vals)
)

# Mode of Gamma(shape, rate) = max(0, (shape - 1)/rate)
gamma_mode <- function(shape, rate) {
  ifelse(shape >= 1, (shape - 1) / rate, 0)
}

# ==============================================================================
# 2. Key functions
# ==============================================================================

# Detection threshold D = m_psi + D0
D_threshold <- function(psi, D0) {
  a <- psi * alpha
  gamma_mode(a, beta) + D0
}

# Tilted survival: pgamma(D, a, beta+r, lower=FALSE)
S_tilted <- function(psi, D, r_val) {
  a <- psi * alpha
  pgamma(D, shape = a, rate = beta + r_val, lower.tail = FALSE)
}

# Untilted survival: pgamma(D, a, beta, lower=FALSE)
S_untilted <- function(psi, D) {
  a <- psi * alpha
  pgamma(D, shape = a, rate = beta, lower.tail = FALSE)
}

# R_eff from the overlap rule
R_eff_fn <- function(R0, psi, D, eta_val) {
  R0 * (1 - eta_val * S_untilted(psi, D))
}

# Naive Euler-Lotka (closed-form): r = beta * (R_eff^{1/alpha} - 1)
r_naive_fn <- function(R_eff_val) {
  if (R_eff_val <= 1) return(0)
  beta * (R_eff_val^(1/alpha) - 1)
}

# Actual Euler-Lotka (root-finding):
# R0 * (beta/(beta+r))^alpha * [1 - eta * S_tilted(psi, D, r)] = 1
r_actual_fn <- function(R0, psi, D, eta_val, r_range = c(-0.5, 10)) {
  a <- psi * alpha

  euler_lotka <- function(r_val) {
    mgf <- (beta / (beta + r_val))^alpha
    surv <- pgamma(D, shape = a, rate = beta + r_val, lower.tail = FALSE)
    R0 * mgf * (1 - eta_val * surv) - 1
  }

  # Check bracketing
  f_lo <- euler_lotka(r_range[1])
  f_hi <- euler_lotka(r_range[2])

  if (is.na(f_lo) || is.na(f_hi)) return(NA_real_)
  if (sign(f_lo) == sign(f_hi)) {
    # Try wider range
    r_range <- c(-1, 20)
    f_lo <- euler_lotka(r_range[1])
    f_hi <- euler_lotka(r_range[2])
    if (is.na(f_lo) || is.na(f_hi)) return(NA_real_)
    if (sign(f_lo) == sign(f_hi)) return(NA_real_)
  }

  uniroot(euler_lotka, interval = r_range, tol = 1e-10)$root
}

# Distorted GI density via numerical convolution
# w*(tau) = [(1-eta)*g(tau) + eta * conv(f_l, f_eps_trunc)(tau)] / (1 - TE)
w_star_fn <- function(tau_grid, psi, D0_val, eta_val = 1) {
  a <- psi * alpha
  b <- (1 - psi) * alpha
  D <- D_threshold(psi, D0_val)

  # Natural GI density
  g_nat <- dgamma(tau_grid, shape = alpha, rate = beta)

  # TE = eta * P(eps > D)
  TE <- eta_val * pgamma(D, shape = a, rate = beta, lower.tail = FALSE)

  if (TE < 1e-12) return(g_nat)  # no intervention effect

  dtau <- tau_grid[2] - tau_grid[1]

  # When b ≈ 0 (ψ ≈ 1), f_l is a point mass at 0, so tau = eps directly.
  # The truncated density is just f_eps(tau) * I(tau < D).
  if (b < 1e-6) {
    density_isolated <- dgamma(tau_grid, shape = a, rate = beta) * (tau_grid < D)
  } else {
    # Convolution: integrate f_l(tau - eps) * f_eps(eps) over eps in [0, min(tau, D)]
    density_isolated <- numeric(length(tau_grid))

    eps_grid <- seq(0, max(D, 0.01), length.out = 500)
    deps <- eps_grid[2] - eps_grid[1]

    f_eps_vals <- dgamma(eps_grid, shape = a, rate = beta)
    f_eps_vals[eps_grid >= D] <- 0

    for (i in seq_along(tau_grid)) {
      tau_i <- tau_grid[i]
      if (tau_i <= 0) next
      eps_use <- eps_grid[eps_grid < tau_i]
      if (length(eps_use) == 0) next
      f_l_vals <- dgamma(tau_i - eps_use, shape = b, rate = beta)
      f_e_vals <- f_eps_vals[seq_along(eps_use)]
      density_isolated[i] <- sum(f_l_vals * f_e_vals) * deps
    }
  }

  # Combined density (unnormalised by R0 contributions):
  # (1 - eta) * g(tau) + eta * density_isolated
  w_combined <- (1 - eta_val) * g_nat + eta_val * density_isolated
  # Normalise
  w_combined <- w_combined / (sum(w_combined) * dtau)
  w_combined
}

# Mean GI under isolation (analytical formula)
mean_gi_star <- function(psi, D0_val, eta_val = 1) {
  a <- psi * alpha
  b <- (1 - psi) * alpha
  D <- D_threshold(psi, D0_val)

  F_D_a   <- pgamma(D, shape = a, rate = beta)
  Fbar_D_a <- 1 - F_D_a
  F_D_a1  <- pgamma(D, shape = a + 1, rate = beta)

  TE <- eta_val * Fbar_D_a

  mean_natural <- alpha / beta

  # Mean of surviving transmissions from isolated individuals
  # E[tau* | isolated, surviving] = b/beta + (a/beta) * F(D; a+1, beta) / F(D; a, beta)
  if (F_D_a < 1e-15) {
    # All transmissions averted — shouldn't happen if R_eff > 0
    return(mean_natural)
  }
  mean_eps_trunc <- (a / beta) * F_D_a1 / F_D_a
  mean_isolated_surviving <- b / beta + mean_eps_trunc

  # Weighted average
  numerator <- (1 - eta_val) * mean_natural + eta_val * F_D_a * mean_isolated_surviving
  denominator <- 1 - TE
  numerator / denominator
}

# ==============================================================================
# 3. Panel A — Distorted GI densities
# ==============================================================================

cat("Computing Panel A: distorted GI densities...\n")

D0_panelA <- 1.0
tau_grid <- seq(0.01, 12, length.out = 600)

df_A <- map_dfr(psi_panel_A, function(psi) {
  w_nat  <- dgamma(tau_grid, shape = alpha, rate = beta)
  w_dist <- w_star_fn(tau_grid, psi, D0_panelA, eta)

  bind_rows(
    tibble(tau = tau_grid, density = w_nat,  type = "Natural", psi = factor(psi)),
    tibble(tau = tau_grid, density = w_dist, type = "Under isolation", psi = factor(psi))
  )
})

panel_A <- ggplot(df_A, aes(x = tau, y = density, colour = type, linetype = type)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~psi, nrow = 1, scales = "free_y",
             labeller = as_labeller(
               setNames(paste0("psi == ", psi_panel_A), as.character(psi_panel_A)),
               label_parsed
             )) +
  scale_colour_manual(values = c("Natural" = "grey40",
                                 "Under isolation" = "#C62828"),
                      name = NULL) +
  scale_linetype_manual(values = c("Natural" = "dashed",
                                   "Under isolation" = "solid"),
                        name = NULL) +
  labs(x = expression(tau ~~ "(generation interval, days)"),
       y = "Density") +
  coord_cartesian(xlim = c(0, 12)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 11))

cat("  Panel A done.\n")

# ==============================================================================
# 4. Panel B — Mean GI shift vs ψ for various D₀
# ==============================================================================

cat("Computing Panel B: mean GI shift...\n")

mean_natural <- alpha / beta

df_B <- map_dfr(D0_vals, function(d0) {
  map_dfr(psi_fine, function(psi) {
    m_star <- mean_gi_star(psi, d0, eta)
    tibble(
      psi = psi,
      delta_mean_pct = (m_star - mean_natural) / mean_natural * 100,
      D0 = factor(d0)
    )
  })
})

label_B <- df_B %>%
  filter(psi == max(psi_fine)) %>%
  mutate(label = paste0("D[0] == ", as.character(D0)))

panel_B <- ggplot(df_B, aes(x = psi, y = delta_mean_pct, colour = D0)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_text(data = label_B,
            aes(x = psi + 0.02, y = delta_mean_pct, label = label),
            parse = TRUE, hjust = 0, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = D0_colors, guide = "none") +
  labs(x = expression(psi),
       y = expression(Delta * "E[" * tau * "] / E[" * tau * "]  (%)")) +
  coord_cartesian(xlim = c(0, 1.2)) +
  theme_classic(base_size = 11)

cat("  Panel B done.\n")

# ==============================================================================
# 5. Panel C — Growth rate bias: Δr = r_actual − r_naive vs ψ
# ==============================================================================

cat("Computing Panel C: growth rate bias...\n")

df_C <- map_dfr(R0_vals, function(R0) {
  map_dfr(psi_fine, function(psi) {
    D <- D_threshold(psi, 1.0)  # D0 = 1
    R_eff <- R_eff_fn(R0, psi, D, eta)

    r_n <- r_naive_fn(R_eff)
    r_a <- r_actual_fn(R0, psi, D, eta)

    if (is.na(r_a)) r_a <- r_n  # fallback

    tibble(
      psi = psi,
      delta_r = r_a - r_n,
      R0 = factor(R0)
    )
  })
})

label_C <- df_C %>%
  filter(psi == max(psi_fine)) %>%
  mutate(label = paste0("R[0] == ", as.character(R0)))

panel_C <- ggplot(df_C, aes(x = psi, y = delta_r, colour = R0)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_text(data = label_C,
            aes(x = psi + 0.02, y = delta_r, label = label),
            parse = TRUE, hjust = 0, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = R0_colors, guide = "none") +
  labs(x = expression(psi),
       y = expression(Delta * r == r[actual] - r[naive])) +
  coord_cartesian(xlim = c(0, 1.2)) +
  theme_classic(base_size = 11)

cat("  Panel C done.\n")

# ==============================================================================
# 6. Panel D — Relative bias: r_actual / r_naive (% overestimate)
# ==============================================================================

cat("Computing Panel D: relative bias...\n")

df_D <- map_dfr(R0_vals, function(R0) {
  map_dfr(psi_fine, function(psi) {
    D <- D_threshold(psi, 1.0)
    R_eff <- R_eff_fn(R0, psi, D, eta)

    r_n <- r_naive_fn(R_eff)
    r_a <- r_actual_fn(R0, psi, D, eta)

    if (is.na(r_a)) r_a <- r_n

    rel_bias <- if (abs(r_n) > 1e-8) (r_a - r_n) / r_n * 100 else 0

    tibble(
      psi = psi,
      rel_bias_pct = rel_bias,
      R0 = factor(R0)
    )
  })
})

label_D <- df_D %>%
  filter(psi == max(psi_fine)) %>%
  mutate(label = paste0("R[0] == ", as.character(R0)))

panel_D <- ggplot(df_D, aes(x = psi, y = rel_bias_pct, colour = R0)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_text(data = label_D,
            aes(x = psi + 0.02, y = rel_bias_pct, label = label),
            parse = TRUE, hjust = 0, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = R0_colors, guide = "none") +
  labs(x = expression(psi),
       y = expression("Relative bias" ~~ (r[actual]/r[naive] - 1) %*% 100 ~ "%")) +
  coord_cartesian(xlim = c(0, 1.2)) +
  theme_classic(base_size = 11)

cat("  Panel D done.\n")

# ==============================================================================
# 7. Combine and save
# ==============================================================================

cat("Assembling figure...\n")

fig <- (panel_A) / (panel_B | panel_C | panel_D) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

ggsave("figures/fig_gi_distortion_analytical.pdf", fig,
       width = 13, height = 8)

cat("Saved figures/fig_gi_distortion_analytical.pdf\n")

# ==============================================================================
# 8. Verification checks
# ==============================================================================

cat("\n--- Verification checks ---\n")

# Check 1: Δr → 0 as ψ → 0
D_check <- D_threshold(0.05, 1.0)
R_eff_check <- R_eff_fn(2, 0.05, D_check, 1)
r_n_check <- r_naive_fn(R_eff_check)
r_a_check <- r_actual_fn(2, 0.05, D_check, 1)
cat(sprintf("  ψ=0.05: Δr = %.6f (should be ~0)\n", r_a_check - r_n_check))

# Check 2: Δr > 0 for ψ = 1 with R0 > 1
D_check2 <- D_threshold(1.0, 1.0)
R_eff_check2 <- R_eff_fn(2, 1.0, D_check2, 1)
r_n_check2 <- r_naive_fn(R_eff_check2)
r_a_check2 <- r_actual_fn(2, 1.0, D_check2, 1)
cat(sprintf("  ψ=1.0, R0=2: Δr = %.6f (should be > 0)\n", r_a_check2 - r_n_check2))

# Check 3: naive = actual when η = 0
r_n_check3 <- r_naive_fn(2)
r_a_check3 <- r_actual_fn(2, 0.5, D_threshold(0.5, 1.0), 0)
cat(sprintf("  η=0: r_naive=%.6f, r_actual=%.6f, diff=%.2e (should be ~0)\n",
            r_n_check3, r_a_check3, abs(r_a_check3 - r_n_check3)))

cat("\n=== Analysis complete ===\n")
