library(tidyverse)
source("code/utils.R")

# ==============================================================================
# Analytical computation of Var(W) and generalized gamma moment matching
#
# Derives closed-form expressions for the variance of the random variable W
# from the branching-process limit theorem, as a function of kappa in the
# Gamma convolutional family. Validates against branching-process simulation
# and fits generalized gamma distributions via moment matching.
#
# Produces:
#   figures/fig_W_variance.pdf      — CV²(W), Var(S) decomposition, time delay
#   figures/fig_W_gengamma.pdf      — GenGamma fits overlaid on simulation
#   figures/fig_W_timeshift.pdf     — Time-shift distributions by kappa
# ==============================================================================

set.seed(2024)

# ==============================================================================
# Section 1: Parameters and analytical formulas
# ==============================================================================

cat("=== Section 1: Analytical formulas for W moments ===\n\n")

# Standard parameters
mu <- 5            # mean generation time
R0 <- 2            # basic reproduction number
alpha_total <- 10  # shape of population kernel Gamma
r <- alpha_total / mu  # rate parameter (= 2)
kappas <- c(0.5, 1, 2, 4, 6, 8, 9.5)

# Closed-form Malthusian parameter
# Euler-Lotka: R0 * (r/(r+alpha))^alpha_total = 1
# => alpha = r * (R0^(1/alpha_total) - 1)
compute_alpha <- function(R0, alpha_total, r) {
  r * (R0^(1/alpha_total) - 1)
}

# m2 = E[sum_j e^{-2*alpha*tau_j}] = R0 * rho2^alpha_total (kappa-independent)
compute_m2 <- function(R0, alpha_total, r, alpha) {
  rho2 <- r / (r + 2 * alpha)
  R0 * rho2^alpha_total
}

# Var(S) as function of kappa
# S = sum_j e^{-alpha*tau_j}, decompose tau_j = s + epsilon_j
# Var(S) = m2 + R0^2 * rho^{2*kappa} * [rho2^{alpha_total - kappa} - rho^{2*(alpha_total - kappa)}]
compute_var_S <- function(R0, alpha_total, kappa, r, alpha) {
  rho  <- r / (r + alpha)
  rho2 <- r / (r + 2 * alpha)

  # Poisson+jitter term (= m2, constant across kappa)
  poisson_jitter <- R0 * rho2^alpha_total

  # Shared-shift synchronisation penalty
  sync_penalty <- R0^2 * rho^(2 * kappa) * (rho2^(alpha_total - kappa) - rho^(2 * (alpha_total - kappa)))

  poisson_jitter + sync_penalty
}

# CV²(W) = Var(S) / (1 - m2)
compute_CV2_W <- function(R0, alpha_total, kappa, r, alpha) {
  var_S <- compute_var_S(R0, alpha_total, kappa, r, alpha)
  m2 <- compute_m2(R0, alpha_total, r, alpha)
  var_S / (1 - m2)
}

# Third moment of W via m3, m21, m111
compute_third_moment_W <- function(R0, alpha_total, kappa, r, alpha) {
  rho  <- r / (r + alpha)
  rho2 <- r / (r + 2 * alpha)
  rho3 <- r / (r + 3 * alpha)

  m2 <- R0 * rho2^alpha_total
  m3 <- R0 * rho3^alpha_total

  # m21 = E[sum_{j!=k} e^{-2*alpha*tau_j} * e^{-alpha*tau_k}]
  # For Poisson(R0) offspring with shared shift s ~ Gamma(alpha_total-kappa, r):
  # E[N(N-1)] = R0^2 (Poisson factorial moment)
  # m21 = R0^2 * rho2^kappa * rho^kappa * rho3^{alpha_total - kappa}
  m21 <- R0^2 * rho2^kappa * rho^kappa * rho3^(alpha_total - kappa)

  # m111 = E[sum_{j!=k!=l} e^{-alpha*(tau_j+tau_k+tau_l)}]
  # E[N(N-1)(N-2)] = R0^3 (Poisson factorial moment)
  # = R0^3 * rho^{3*kappa} * rho3^{alpha_total - kappa}
  m111 <- R0^3 * rho^(3 * kappa) * rho3^(alpha_total - kappa)

  # E[W] = 1 (by normalisation), E[W^2] = 1 + CV2
  CV2 <- compute_CV2_W(R0, alpha_total, kappa, r, alpha)
  EW <- 1
  EW2 <- 1 + CV2

  # E[W^3] = (3*m21*EW2*EW + m111*EW^3 + m3) / (1 - m3)
  # More carefully: from the recursive moment equation
  # E[W^3] = [m3 + 3*m21*EW2 + m111] / (1 - m3)
  # But need to account for the mean contributions properly.
  # The third-moment recursion for W (with E[W]=1) is:
  # E[W^3](1 - m3) = E[S^3] where E[S^3] includes all cross-terms
  # E[S^3] = m3 + 3*m21*EW + m111*EW  ... no, let me be more careful.
  #
  # W =d sum_j D_j W_j where D_j = e^{-alpha*tau_j}
  # E[W^3] = E[(sum D_j W_j)^3]
  #         = sum_j E[D_j^3] E[W_j^3]                           (j=k=l terms)
  #           + 3 sum_{j!=k} E[D_j^2 D_k] E[W^2] E[W]          (j=k!=l terms)
  #           + sum_{j!=k!=l} E[D_j D_k D_l] (E[W])^3           (all distinct)
  # E[W^3] = m3 * E[W^3] + 3 * m21 * EW2 * EW + m111 * EW^3
  # => E[W^3] = (3 * m21 * EW2 + m111) / (1 - m3)

  EW3 <- (3 * m21 * EW2 * EW + m111 * EW^3) / (1 - m3)

  EW3
}

# Extinction probability: solve q = exp(R0*(q-1)) numerically
compute_extinction_prob <- function(R0) {
  if (R0 <= 1) return(1)
  # Fixed-point iteration
  q <- 0.5
  for (i in 1:1000) {
    q_new <- exp(R0 * (q - 1))
    if (abs(q_new - q) < 1e-12) break
    q <- q_new
  }
  q
}

# Time delay: DeltaE[T_n] approx (1-q) * [Var(S)_kappa - Var(S)_smooth] / [2*alpha*(1-m2)]
compute_delay <- function(R0, alpha_total, kappa, r, alpha, q) {
  var_S_kappa  <- compute_var_S(R0, alpha_total, kappa, r, alpha)
  var_S_smooth <- compute_var_S(R0, alpha_total, alpha_total * 0.9999, r, alpha)
  m2 <- compute_m2(R0, alpha_total, r, alpha)

  # The delay formula from Jensen's inequality on log(W):
  # E[log W]_kappa - E[log W]_smooth ≈ -[Var(W)_kappa - Var(W)_smooth] / (2 * (E[W])^2)
  # Since E[W] = 1/(1-q) for surviving lineages (conditional), and
  # tau = log(W)/alpha, we get:
  # Delta E[T_n] ≈ [CV2_kappa - CV2_smooth] / (2*alpha)
  # where CV2 = Var(W)/(E[W])^2 = Var(S)/(1-m2)

  CV2_kappa  <- var_S_kappa / (1 - m2)
  CV2_smooth <- var_S_smooth / (1 - m2)

  (CV2_kappa - CV2_smooth) / (2 * alpha)
}

# Compute all analytical quantities
alpha <- compute_alpha(R0, alpha_total, r)
m2 <- compute_m2(R0, alpha_total, r, alpha)
q <- compute_extinction_prob(R0)

cat(sprintf("  Malthusian parameter alpha = %.6f\n", alpha))
cat(sprintf("  m2 = E[sum e^{-2*alpha*tau_j}] = %.6f\n", m2))
cat(sprintf("  Extinction probability q = %.6f\n", q))
cat(sprintf("  rho = r/(r+alpha) = %.6f\n", r / (r + alpha)))
cat(sprintf("  rho2 = r/(r+2*alpha) = %.6f\n", r / (r + 2 * alpha)))
cat("\n")

# Fine grid for curves
kappa_grid <- seq(0.1, alpha_total - 0.1, length.out = 200)
CV2_curve <- sapply(kappa_grid, function(k) compute_CV2_W(R0, alpha_total, k, r, alpha))
delay_curve <- sapply(kappa_grid, function(k) compute_delay(R0, alpha_total, k, r, alpha, q))

# Decomposition of Var(S)
rho  <- r / (r + alpha)
rho2 <- r / (r + 2 * alpha)
poisson_jitter_grid <- rep(m2, length(kappa_grid))
sync_penalty_grid <- sapply(kappa_grid, function(k) {
  R0^2 * rho^(2 * k) * (rho2^(alpha_total - k) - rho^(2 * (alpha_total - k)))
})

# Summary table for standard kappa values
cat("=== Analytical quantities by kappa ===\n\n")
cat(sprintf("  %-6s  %-10s  %-10s  %-10s  %-10s\n", "kappa", "Var(S)", "CV2(W)", "E[W^3]", "Delay"))
cat(sprintf("  %-6s  %-10s  %-10s  %-10s  %-10s\n", "-----", "------", "------", "------", "-----"))
for (k in kappas) {
  vs <- compute_var_S(R0, alpha_total, k, r, alpha)
  cv2 <- compute_CV2_W(R0, alpha_total, k, r, alpha)
  ew3 <- compute_third_moment_W(R0, alpha_total, k, r, alpha)
  del <- compute_delay(R0, alpha_total, k, r, alpha, q)
  cat(sprintf("  %-6.1f  %-10.6f  %-10.6f  %-10.4f  %-10.4f\n", k, vs, cv2, ew3, del))
}
cat("\n")

# ==============================================================================
# Section 2: Generalized gamma moment matching
# ==============================================================================

cat("=== Section 2: Generalized gamma moment matching ===\n\n")

# Generalized gamma with parameters (a, d, p):
# f(x) = (p/a) * (x/a)^{d-1} * exp(-(x/a)^p) / Gamma(d/p)
# E[X^k] = a^k * Gamma((d+k)/p) / Gamma(d/p)
#
# We match the first three moments of W|W>0 (conditional on survival).
# Since E[W] = 1 (unconditional), E[W|survive] = 1/(1-q).

fit_gengamma <- function(EW, EW2, EW3) {
  # Given E[X], E[X^2], E[X^3], find (a, d, p) for GenGamma
  # Use moment ratios to reduce to 2D optimisation over (d, p).
  #
  # GenGamma moments: mu_k = a^k * Gamma((d+k)/p) / Gamma(d/p)
  # Ratios (eliminating a):
  #   r2 = mu2/mu1^2 = Gamma((d+2)/p)*Gamma(d/p) / Gamma((d+1)/p)^2
  #   r3 = mu3/mu1^3 = Gamma((d+3)/p)*Gamma(d/p)^2 / Gamma((d+1)/p)^3
  #
  # Use lgamma for numerical stability.

  target_lr2 <- log(EW2 / EW^2)
  target_lr3 <- log(EW3 / EW^3)

  obj <- function(par) {
    d <- exp(par[1])
    p <- exp(par[2])

    if (d < 1e-4 || p < 1e-4 || d > 50 || p > 50) return(1e10)

    lr2 <- lgamma((d + 2)/p) + lgamma(d/p) - 2 * lgamma((d + 1)/p)
    lr3 <- lgamma((d + 3)/p) + 2 * lgamma(d/p) - 3 * lgamma((d + 1)/p)

    if (is.nan(lr2) || is.nan(lr3) ||
        is.infinite(lr2) || is.infinite(lr3)) return(1e10)

    (lr2 - target_lr2)^2 + (lr3 - target_lr3)^2
  }

  # Try multiple starting points
  best <- list(value = Inf)
  for (d0 in c(0.3, 0.7, 1, 2, 5, 10)) {
    for (p0 in c(0.3, 0.5, 1, 1.5, 2, 3, 5)) {
      fit <- tryCatch(
        optim(c(log(d0), log(p0)), obj, method = "Nelder-Mead",
              control = list(maxit = 10000, reltol = 1e-14)),
        error = function(e) list(value = Inf)
      )
      if (fit$value < best$value) best <- fit
    }
  }

  d <- exp(best$par[1])
  p <- exp(best$par[2])
  # Recover a from E[X] = a * Gamma((d+1)/p) / Gamma(d/p)
  a <- EW * exp(lgamma(d/p) - lgamma((d + 1)/p))

  list(a = a, d = d, p = p, convergence = best$convergence, obj_value = best$value)
}

# Fit GenGamma for each kappa
cat(sprintf("  %-6s  %-10s  %-10s  %-10s  %-10s\n", "kappa", "a", "d", "p", "obj"))
cat(sprintf("  %-6s  %-10s  %-10s  %-10s  %-10s\n", "-----", "---", "---", "---", "---"))

gengamma_fits <- list()
for (k in kappas) {
  CV2 <- compute_CV2_W(R0, alpha_total, k, r, alpha)
  EW3 <- compute_third_moment_W(R0, alpha_total, k, r, alpha)

  # Unconditional moments: E[W]=1, E[W^2]=1+CV2, E[W^3]
  # W = 0 with prob q (extinction), W > 0 with prob (1-q)
  # E[W^k] = (1-q) * E[W^k | survive] + q * 0^k (for k>=1, 0^k=0)
  # So E[W^k | survive] = E[W^k] / (1-q)
  EW2_uncond <- 1 + CV2  # since E[W]=1
  EW_cond  <- 1 / (1 - q)
  EW2_cond <- EW2_uncond / (1 - q)
  EW3_cond <- EW3 / (1 - q)

  fit <- fit_gengamma(EW_cond, EW2_cond, EW3_cond)
  gengamma_fits[[as.character(k)]] <- fit

  cat(sprintf("  %-6.1f  %-10.6f  %-10.6f  %-10.6f  %-10.2e\n",
              k, fit$a, fit$d, fit$p, fit$obj_value))
}
cat("\n")

# ==============================================================================
# Section 3: Branching-process simulation for validation
# ==============================================================================

cat("=== Section 3: Branching-process simulation validation ===\n\n")

N_sim <- 1000
n_pop <- 5000   # depletion < 4% at threshold
threshold <- 200  # deep into exponential phase

sim_results <- list()

for (k in kappas) {
  cat(sprintf("  Simulating kappa = %.1f ...\n", k))

  profile_fn <- make_profile_gamma(mu = mu, R0 = R0, alpha_total = alpha_total, kappa = k)

  W_vals <- numeric(N_sim)
  T_threshold <- numeric(N_sim)
  survived <- logical(N_sim)

  for (i in 1:N_sim) {
    tinf <- sim_stochastic_fast(n = n_pop, gen_inf_attempts = profile_fn)

    # Sort infection times (finite only)
    finite_times <- sort(tinf[is.finite(tinf)])
    n_infected <- length(finite_times)

    # Check survival: at least threshold infections
    if (n_infected >= threshold) {
      survived[i] <- TRUE
      T_threshold[i] <- finite_times[threshold]

      # Estimate W = Z(T)/e^{alpha*T} at the threshold time
      # Z(T) = number infected by time T = threshold (by definition)
      W_vals[i] <- threshold * exp(-alpha * T_threshold[i])
    } else {
      survived[i] <- FALSE
      W_vals[i] <- NA
      T_threshold[i] <- NA
    }
  }

  surv_idx <- which(survived)
  W_surv <- W_vals[surv_idx]
  T_surv <- T_threshold[surv_idx]

  sim_results[[as.character(k)]] <- list(
    kappa = k,
    W = W_surv,
    T_threshold = T_surv,
    n_survived = length(surv_idx),
    survival_rate = mean(survived),
    EW = mean(W_surv),
    VarW = var(W_surv),
    CV2W = var(W_surv) / mean(W_surv)^2,
    EW2 = mean(W_surv^2),
    EW3 = mean(W_surv^3),
    mean_T = mean(T_surv),
    sd_T = sd(T_surv)
  )
}

# Conditional CV²(W|surv) from unconditional:
# CV²(W|surv) = (1-q)*CV²(W) - q
compute_CV2_W_cond <- function(R0, alpha_total, kappa, r, alpha, q) {
  (1 - q) * compute_CV2_W(R0, alpha_total, kappa, r, alpha) - q
}

# Compute simulated delays relative to smoothest kappa
ref_k <- "9.5"
ref_mean_T <- sim_results[[ref_k]]$mean_T

cat("\n=== Validation: analytical vs simulated ===\n")
cat("  (comparing conditional CV2(W|survival), which is what simulations estimate)\n\n")
cat(sprintf("  %-6s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n",
            "kappa", "CV2c_anal", "CV2c_sim", "relErr%", "delay_anal", "delay_sim", "relErr%"))
cat(sprintf("  %-6s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n",
            "-----", "---------", "--------", "-------", "----------", "---------", "-------"))

for (k in kappas) {
  ks <- as.character(k)
  cv2c_anal <- compute_CV2_W_cond(R0, alpha_total, k, r, alpha, q)
  cv2c_sim  <- sim_results[[ks]]$CV2W
  rel_cv2   <- 100 * (cv2c_sim - cv2c_anal) / cv2c_anal

  del_anal <- compute_delay(R0, alpha_total, k, r, alpha, q)
  del_sim  <- sim_results[[ks]]$mean_T - ref_mean_T
  rel_del  <- if (abs(del_anal) > 0.001) 100 * (del_sim - del_anal) / del_anal else NA

  cat(sprintf("  %-6.1f  %-10.6f  %-10.6f  %-10.1f  %-10.4f  %-10.4f  %-10.1f\n",
              k, cv2c_anal, cv2c_sim, rel_cv2, del_anal, del_sim,
              if (is.na(rel_del)) NA else rel_del))
}
cat("\n")

# ==============================================================================
# Section 4: Figures
# ==============================================================================

cat("=== Section 4: Generating figures ===\n\n")

# --- Figure 1: W variance (3 panels) ---

# For figures: use conditional CV²(W|surv) since simulations only observe survivors
CV2_cond_curve <- sapply(kappa_grid, function(k) compute_CV2_W_cond(R0, alpha_total, k, r, alpha, q))

sim_cv2 <- sapply(kappas, function(k) sim_results[[as.character(k)]]$CV2W)
sim_cv2_se <- sapply(kappas, function(k) {
  W <- sim_results[[as.character(k)]]$W
  n <- length(W)
  # Bootstrap SE of CV^2
  boot_cv2 <- replicate(500, {
    Wb <- sample(W, n, replace = TRUE)
    var(Wb) / mean(Wb)^2
  })
  sd(boot_cv2)
})

sim_delay <- sapply(kappas, function(k) sim_results[[as.character(k)]]$mean_T - ref_mean_T)
sim_delay_se <- sapply(kappas, function(k) {
  sim_results[[as.character(k)]]$sd_T / sqrt(sim_results[[as.character(k)]]$n_survived)
})

pdf("figures/fig_W_variance.pdf", width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2, 1), oma = c(0, 0, 2, 0))

# Panel A: CV²(W|surv) vs kappa
plot(kappa_grid, CV2_cond_curve, type = "l", lwd = 2, col = "steelblue",
     xlab = expression(kappa), ylab = expression(CV^2*(W~"|"~surv)),
     main = expression("(A)"~~CV^2*(W~"|"~survival)))
points(kappas, sim_cv2, pch = 16, col = "firebrick", cex = 1.2)
arrows(kappas, sim_cv2 - 1.96 * sim_cv2_se, kappas, sim_cv2 + 1.96 * sim_cv2_se,
       angle = 90, code = 3, length = 0.05, col = "firebrick")
legend("topright", legend = c("Analytical", "Simulated"),
       col = c("steelblue", "firebrick"), lty = c(1, NA), pch = c(NA, 16),
       lwd = c(2, NA), bty = "n")

# Panel B: Var(S) decomposition
plot(kappa_grid, poisson_jitter_grid + sync_penalty_grid, type = "l", lwd = 2,
     col = "steelblue", xlab = expression(kappa), ylab = "Var(S)",
     main = "(B) Var(S) decomposition", ylim = c(0, max(poisson_jitter_grid + sync_penalty_grid) * 1.1))
# Shaded regions
polygon(c(kappa_grid, rev(kappa_grid)),
        c(rep(0, length(kappa_grid)), rev(poisson_jitter_grid)),
        col = adjustcolor("grey60", 0.4), border = NA)
polygon(c(kappa_grid, rev(kappa_grid)),
        c(poisson_jitter_grid, rev(poisson_jitter_grid + sync_penalty_grid)),
        col = adjustcolor("coral", 0.4), border = NA)
lines(kappa_grid, poisson_jitter_grid + sync_penalty_grid, lwd = 2, col = "steelblue")
lines(kappa_grid, poisson_jitter_grid, lwd = 1, lty = 2, col = "grey40")
legend("topright", legend = c("Poisson + jitter (constant)", "Synchronisation penalty"),
       fill = c(adjustcolor("grey60", 0.4), adjustcolor("coral", 0.4)),
       bty = "n", cex = 0.9)

# Panel C: Expected time delay (analytical only — simulation noise exceeds signal)
plot(kappa_grid, delay_curve, type = "l", lwd = 2, col = "steelblue",
     xlab = expression(kappa), ylab = expression(Delta*E*"["*T[n]*"]"~~"(days)"),
     main = expression("(C) Expected time delay"~~Delta*E*"["*T[n]*"]"))
abline(h = 0, lty = 3, col = "grey50")
# Annotate max delay
max_del <- max(delay_curve)
text(0.5, max_del * 0.85, sprintf("Max delay: %.2f days", max_del),
     adj = c(0, 0), cex = 0.85, col = "steelblue")

mtext(expression("Variance of W and growth delay vs"~kappa), outer = TRUE, cex = 1.2)
dev.off()
cat("  Saved figures/fig_W_variance.pdf\n")

# --- Figure 2: GenGamma fits (4 panels) ---

# Generalized gamma density (log-scale computation for stability)
dgengamma <- function(x, a, d, p) {
  log_dens <- log(p) - log(a) + (d - 1) * log(x / a) - (x / a)^p - lgamma(d / p)
  exp(log_dens)
}

show_kappas <- c(0.5, 2, 6, 9.5)

# For the density overlay figure, fit GenGamma to *empirical* moments from simulation.
# The simulated W_est = n * e^{-alpha*T_n} estimates W up to an unknown scaling constant c,
# so the analytical GenGamma (fit to theoretical moments with E[W]=1) has a different scale.
# Fitting to empirical moments gives a properly scaled comparison.
gengamma_fits_emp <- list()
for (k in show_kappas) {
  ks <- as.character(k)
  W_surv <- sim_results[[ks]]$W
  gengamma_fits_emp[[ks]] <- fit_gengamma(mean(W_surv), mean(W_surv^2), mean(W_surv^3))
}

pdf("figures/fig_W_gengamma.pdf", width = 12, height = 4)
par(mfrow = c(1, 4), mar = c(4, 4, 3, 1))

for (k in show_kappas) {
  ks <- as.character(k)
  W_surv <- sim_results[[ks]]$W
  fit_emp <- gengamma_fits_emp[[ks]]

  # Use kernel density for smoother empirical estimate
  d_emp <- density(W_surv, from = 0, adjust = 1.2)
  x_upper <- quantile(W_surv, 0.98)

  # GenGamma density from empirical moment fit
  x_seq <- seq(0.01, x_upper, length.out = 500)
  gg_dens <- dgengamma(x_seq, fit_emp$a, fit_emp$d, fit_emp$p)

  # Determine y range from both
  ymax <- max(c(d_emp$y[d_emp$x <= x_upper], gg_dens), na.rm = TRUE) * 1.05

  # Plot empirical density as shaded area
  plot(d_emp, col = "steelblue", lwd = 2, xlim = c(0, x_upper), ylim = c(0, ymax),
       main = bquote(kappa == .(k)), xlab = "W | survival", ylab = "Density")
  polygon(c(d_emp$x, rev(d_emp$x)),
          c(d_emp$y, rep(0, length(d_emp$y))),
          col = adjustcolor("steelblue", 0.2), border = NA)

  # Overlay GenGamma density
  lines(x_seq, gg_dens, col = "firebrick", lwd = 2, lty = 2)

  # Annotation: show analytical fit params (the ones of interest)
  fit_anal <- gengamma_fits[[ks]]
  legend("topright",
         legend = c("Empirical", "GenGamma fit",
                    sprintf("a=%.3f", fit_anal$a),
                    sprintf("d=%.3f", fit_anal$d),
                    sprintf("p=%.3f", fit_anal$p)),
         col = c("steelblue", "firebrick", NA, NA, NA),
         lty = c(1, 2, NA, NA, NA), lwd = c(2, 2, NA, NA, NA),
         bty = "n", cex = 0.7, text.col = c("steelblue", "firebrick",
                                              "firebrick", "firebrick", "firebrick"))
}

dev.off()
cat("  Saved figures/fig_W_gengamma.pdf\n")

# --- Figure 3: Time-shift distributions ---

pdf("figures/fig_W_timeshift.pdf", width = 8, height = 5)
par(mar = c(4.5, 4.5, 2, 1))

cols <- colorRampPalette(c("firebrick", "steelblue"))(length(show_kappas))

# First pass: compute all densities to set ylim
all_dens <- list()
for (i in seq_along(show_kappas)) {
  k <- show_kappas[i]
  ks <- as.character(k)
  W_surv <- sim_results[[ks]]$W

  # tau = (log(W) - log(E[W])) / alpha
  tau_emp <- (log(W_surv) - log(mean(W_surv))) / alpha
  tau_emp <- tau_emp[is.finite(tau_emp)]
  all_dens[[i]] <- density(tau_emp, adjust = 1.2)
}

ymax <- max(sapply(all_dens, function(d) max(d$y))) * 1.1

plot(NULL, xlim = c(-20, 15), ylim = c(0, ymax),
     xlab = expression("Time shift"~~tau~~"(days)"),
     ylab = "Density",
     main = expression("Time-shift distribution"~~tau == (log*W - log*E*"["*W*"]") / alpha))

for (i in seq_along(show_kappas)) {
  k <- show_kappas[i]
  ks <- as.character(k)

  # Empirical
  lines(all_dens[[i]], col = cols[i], lwd = 2)

  # Analytical from GenGamma fit
  fit <- gengamma_fits[[ks]]
  EW_cond <- 1 / (1 - q)

  # tau = (log(W) - log(EW_cond)) / alpha
  # W ~ GenGamma(a, d, p), so log(W) has a known distribution
  # Sample from GenGamma and transform
  n_draw <- 50000
  # Sample GenGamma: X = a * Y^{1/p} where Y ~ Gamma(d/p, 1)
  Y <- rgamma(n_draw, shape = fit$d / fit$p, rate = 1)
  W_draw <- fit$a * Y^(1 / fit$p)
  tau_gg <- (log(W_draw) - log(EW_cond)) / alpha
  d_gg <- density(tau_gg, adjust = 1.2)
  lines(d_gg, col = cols[i], lwd = 1.5, lty = 2)
}

legend("topleft",
       legend = sapply(show_kappas, function(k) bquote(kappa == .(k))),
       col = cols, lwd = 2, bty = "n")
legend("topright",
       legend = c("Empirical", "GenGamma fit"),
       lty = c(1, 2), lwd = c(2, 1.5), bty = "n")

dev.off()
cat("  Saved figures/fig_W_timeshift.pdf\n")

# ==============================================================================
# Section 5: Console output summary
# ==============================================================================

cat("\n=== Final summary ===\n\n")

cat(sprintf("  Parameters: mu=%.0f, R0=%.0f, alpha_total=%.0f, r=%.1f\n", mu, R0, alpha_total, r))
cat(sprintf("  Malthusian parameter: alpha = %.6f\n", alpha))
cat(sprintf("  m2 (kappa-independent): %.6f\n", m2))
cat(sprintf("  Extinction probability: q = %.6f\n\n", q))

cat("  Analytical quantities:\n")
cat(sprintf("  %-6s  %-10s  %-10s  %-10s\n", "kappa", "CV2(W)", "Delay(d)", "E[W^3]"))
cat(sprintf("  %-6s  %-10s  %-10s  %-10s\n", "-----", "------", "--------", "------"))
for (k in kappas) {
  cv2 <- compute_CV2_W(R0, alpha_total, k, r, alpha)
  del <- compute_delay(R0, alpha_total, k, r, alpha, q)
  ew3 <- compute_third_moment_W(R0, alpha_total, k, r, alpha)
  cat(sprintf("  %-6.1f  %-10.6f  %-10.4f  %-10.4f\n", k, cv2, del, ew3))
}

cat("\n  Generalized gamma fits:\n")
cat(sprintf("  %-6s  %-10s  %-10s  %-10s\n", "kappa", "a", "d", "p"))
cat(sprintf("  %-6s  %-10s  %-10s  %-10s\n", "-----", "---", "---", "---"))
for (k in kappas) {
  fit <- gengamma_fits[[as.character(k)]]
  cat(sprintf("  %-6.1f  %-10.4f  %-10.4f  %-10.4f\n", k, fit$a, fit$d, fit$p))
}

cat("\n  Simulation validation (conditional CV2(W|surv)):\n")
cat(sprintf("  %-6s  %-10s  %-10s  %-10s\n", "kappa", "CV2c_anal", "CV2c_sim", "relErr%"))
cat(sprintf("  %-6s  %-10s  %-10s  %-10s\n", "-----", "---------", "--------", "-------"))
for (k in kappas) {
  ks <- as.character(k)
  cv2c_anal <- compute_CV2_W_cond(R0, alpha_total, k, r, alpha, q)
  cv2c_sim  <- sim_results[[ks]]$CV2W
  rel_err   <- 100 * (cv2c_sim - cv2c_anal) / cv2c_anal
  cat(sprintf("  %-6.1f  %-10.6f  %-10.6f  %-10.1f\n", k, cv2c_anal, cv2c_sim, rel_err))
}

cat("\n  CV2(W) monotonicity check (unconditional, should be decreasing):\n")
cv2_vals <- sapply(kappas, function(k) compute_CV2_W(R0, alpha_total, k, r, alpha))
cat(sprintf("  %s\n", paste(sprintf("%.4f", cv2_vals), collapse = " > ")))
cat(sprintf("  Monotonically decreasing: %s\n", all(diff(cv2_vals) < 0)))

cat("\n  Done.\n")
