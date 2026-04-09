library(ggplot2)
library(dplyr)
library(tidyr)

# --- Core function ---
var_W <- function(z, alpha, psi) {
  c_val   <- (1 + z) / (1 + 2*z)
  sig_val <- (1 + z)^2 / (1 + 2*z)
  
  numer <- c_val^alpha + sig_val^((1 - psi)*alpha) - 1
  denom <- 1 - c_val^alpha
  
  # Return NA outside validity region (denom <= 0 or numer/denom <= 0)
  ifelse(denom > 0 & numer > 0, numer / denom, NA_real_)
}

# --- Grid ---
alpha_seq <- seq(0.5, 10,  length.out = 300)
psi_seq   <- seq(0,   1,   length.out = 300)
z_vals    <- c(0.1, 0.5, 2.0)
z_labels  <- c("z = 0.1  (slow)", "z = 0.5  (moderate)", "z = 2.0  (fast)")

grid <- expand.grid(alpha = alpha_seq, psi = psi_seq, z = z_vals) |>
  mutate(
    VarW  = var_W(z, alpha, psi),
    z_lab = factor(z, levels = z_vals, labels = z_labels)
  )

# --- Hyperbola overlay: (1 - psi)*alpha = k ---
k_vals <- c(0.5, 1, 2, 4, 7)

hyperbolas <- expand.grid(k = k_vals, alpha = alpha_seq) |>
  mutate(psi = 1 - k / alpha) |>
  filter(psi >= 0, psi <= 1) |>
  mutate(k_lab = factor(k))

# Replicate hyperbolas across z panels
hyperbolas <- bind_rows(lapply(z_vals, function(zv) {
  hyperbolas |> mutate(z = zv,
                       z_lab = factor(zv, levels = z_vals, labels = z_labels))
}))

# --- Contour breaks: shared across panels for comparability ---
breaks_contour <- c(0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 25)

# --- Plot ---
p <- ggplot(grid, aes(x = alpha, y = psi)) +

  # Filled contours
  geom_contour_filled(
    aes(z = VarW),
    breaks = c(0, breaks_contour, Inf)
  ) +

  # Labelled contour lines
  geom_contour(
    aes(z = VarW),
    breaks  = breaks_contour,
    colour  = "white",
    linewidth = 0.35,
    alpha   = 0.7
  ) +

  # Hyperbola overlays
  geom_line(
    data      = hyperbolas,
    aes(x = alpha, y = psi, group = k_lab),
    colour    = "black",
    linetype  = "dashed",
    linewidth = 0.5,
    alpha     = 0.6
  ) +

  # Hyperbola labels at right edge
  geom_text(
    data = hyperbolas |>
      group_by(k, z_lab) |>
      filter(alpha == max(alpha)) |>
      ungroup(),
    aes(x = alpha + 0.15, y = psi,
        label = paste0("(1-\u03c8)\u03b1 = ", k)),
    hjust   = 0,
    size    = 2.5,
    colour  = "black",
    alpha   = 0.8
  ) +

  facet_wrap(~ z_lab, nrow = 1) +

  scale_x_continuous(
    name   = expression(alpha ~ "(generation interval shape)"),
    expand = expansion(mult = c(0, 0.12))
  ) +
  scale_y_continuous(
    name   = expression(psi ~ "(punctuation parameter)"),
    limits = c(0, 1),
    expand = expansion(mult = 0)
  ) +
  scale_fill_viridis_d(
    name   = "Var(W)",
    option = "plasma",
    direction = -1,
    # labels = function(x) breaks_contour   # label by lower break
  ) +

  labs(
    title    = expression("Variance of stochastic head-start " * italic(W)),
    subtitle = expression("Dashed curves: lines of constant (1-" * psi * ")" * alpha)
  ) +

  theme_minimal(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    strip.text      = element_text(face = "bold", size = 11),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm"),
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(colour = "grey40", size = 9)
  )

print(p)

# Save
ggsave("fig_varW_contour.pdf", p, width = 12, height = 4.5)
ggsave("fig_varW_contour.png", p, width = 12, height = 4.5, dpi = 300)