library(tidyverse)
library(patchwork)

# ==============================================================================
# Overlap rule visualisation
#
# Three-panel figure illustrating how punctuated infectiousness interacts
# with isolation-based interventions via the mode-centred overlap rule:
#   TE = η · P(det) · P(ε̃ > D₀)
#
# Panel A: Mode-centred individual profiles f_ψ(x + m_ψ) for varying ψ
# Panel B: Fraction averted S_ψ(m_ψ + D₀) vs isolation offset D₀
# Panel C: TE vs ψ for fixed isolation offsets D₀
# ==============================================================================

# ==============================================================================
# 1. Parameters
# ==============================================================================

T_gen    <- 5
popshape <- 10          # α
poprate  <- popshape / T_gen   # β = 2

psi_vals   <- c(0.1, 0.3, 0.5, 1.0)
psi_labels <- setNames(
  paste0("ψ = ", psi_vals),
  as.character(psi_vals)
)

# Color palette: red (spike) → blue (smooth)
psi_colors <- setNames(
  colorRampPalette(c("#C62828", "#E65100", "#1565C0", "#0D47A1"))(length(psi_vals)),
  as.character(psi_vals)
)

# Mode of Gamma(shape, rate) = max(0, (shape - 1)/rate)
gamma_mode <- function(shape, rate) {
  ifelse(shape >= 1, (shape - 1) / rate, 0)
}

# ==============================================================================
# 2. Panel A — Mode-centred density profiles
# ==============================================================================

x_grid <- seq(-3, 5, length.out = 1000)

df_A <- map_dfr(psi_vals, function(psi) {
  shape <- psi * popshape
  rate  <- poprate
  m_psi <- gamma_mode(shape, rate)
  tibble(
    x       = x_grid,
    density = dgamma(x + m_psi, shape = shape, rate = rate),
    psi     = factor(psi)
  )
})

panel_A <- ggplot(df_A, aes(x = x, y = density, colour = psi)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40",
             linewidth = 0.5) +
  annotate("text", x = -1, y = max(df_A$density, na.rm = TRUE) * 0.97,
           label = expression(D[0] == -1), hjust = 1.1, size = 3.2,
           colour = "grey30") +
  annotate("text", x = 1, y = max(df_A$density, na.rm = TRUE) * 0.97,
           label = expression(D[0] == +1), hjust = -0.1, size = 3.2,
           colour = "grey30") +
  scale_colour_manual(
    values = psi_colors,
    labels = function(b) parse(text = paste0("psi == ", b)),
    name   = NULL
  ) +
  labs(
    x = expression("Time relative to peak," ~~ tilde(epsilon) == epsilon - m[psi]),
    y = "Density"
  ) +
  coord_cartesian(xlim = c(-3, 5), ylim = c(0, NA)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        legend.key.width = unit(1.2, "cm"))

# ==============================================================================
# 3. Panel B — Fraction averted vs isolation offset D₀
# ==============================================================================

D0_grid <- seq(-3, 5, length.out = 1000)

df_B <- map_dfr(psi_vals, function(psi) {
  shape <- psi * popshape
  rate  <- poprate
  m_psi <- gamma_mode(shape, rate)
  # S_ψ(m_ψ + D₀) = P(Gamma(ψα, β) > m_ψ + D₀) = pgamma(m_ψ + D₀, shape, rate, lower=FALSE)
  tibble(
    D0             = D0_grid,
    frac_averted   = pgamma(m_psi + D0_grid, shape = shape, rate = rate,
                            lower.tail = FALSE),
    psi            = factor(psi)
  )
})

panel_B <- ggplot(df_B, aes(x = D0, y = frac_averted, colour = psi)) +
  annotate("rect", xmin = -3, xmax = 0, ymin = 0, ymax = 1,
           fill = "#E8F5E9", alpha = 0.5) +
  annotate("text", x = -1.5, y = 0.05, label = "Spike-favoured",
           size = 3, colour = "grey40", fontface = "italic") +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  scale_colour_manual(
    values = psi_colors,
    labels = function(b) parse(text = paste0("psi == ", b)),
    name   = NULL
  ) +
  labs(
    x = expression("Isolation offset from mode," ~~ D[0]),
    y = expression(theta / eta ~~ "(fraction averted)")
  ) +
  coord_cartesian(xlim = c(-3, 5), ylim = c(0, 1)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        legend.key.width = unit(1.2, "cm"))

# ==============================================================================
# 4. Panel C — TE vs ψ for fixed isolation offsets
# ==============================================================================

D0_fixed <- c(-1, 0, 0.5, 1.5)
D0_colors <- c("-1" = "#2E7D32", "0" = "grey30", "0.5" = "#E65100", "1.5" = "#C62828")

psi_fine <- seq(0.05, 1, length.out = 300)

df_C <- map_dfr(D0_fixed, function(d0) {
  map_dfr(psi_fine, function(psi) {
    shape <- psi * popshape
    rate  <- poprate
    m_psi <- gamma_mode(shape, rate)
    tibble(
      psi          = psi,
      frac_averted = pgamma(m_psi + d0, shape = shape, rate = rate,
                            lower.tail = FALSE),
      D0           = factor(d0)
    )
  })
})

# Label positions: at right edge of each curve
label_df <- df_C %>%
  filter(psi == max(psi_fine)) %>%
  mutate(label = paste0("D[0] == ", as.character(D0)))

panel_C <- ggplot(df_C, aes(x = psi, y = frac_averted, colour = D0)) +
  geom_line(linewidth = 0.8) +
  geom_text(data = label_df,
            aes(x = psi + 0.03, y = frac_averted, label = label),
            parse = TRUE, hjust = 0, size = 3.2, show.legend = FALSE) +
  scale_colour_manual(values = D0_colors, guide = "none") +
  labs(
    x = expression(psi),
    y = expression(theta / eta ~~ "(fraction averted)")
  ) +
  coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 1)) +
  theme_classic(base_size = 12)

# ==============================================================================
# 5. Combine and save
# ==============================================================================

fig <- panel_A + panel_B + panel_C +
  plot_layout(nrow = 1) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("figures/fig_overlap_rule.pdf", fig, width = 14, height = 4.5)
cat("Saved figures/fig_overlap_rule.pdf\n")
