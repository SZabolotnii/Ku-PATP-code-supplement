#!/usr/bin/env Rscript
# 02_visualizations.R -- theoretical PATP figures
# Output: PDF figures in ../figures/
# Date: 2026-05-11

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

theme_patp <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 11, face = "plain"),
        axis.title = element_text(size = 10),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(size = 9))

fig_dir <- "../figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ===================================================================
# Figure 1: p_i(alpha) curves for i in {2, 3, 4, 5}
# ===================================================================

p_i <- function(i, alpha) {
  1 / i + (4 - i - 3 / i) * alpha + (2 * i - 4 + 2 / i) * alpha^2
}

alpha_grid <- seq(0, 1, length.out = 201)
pi_df <- expand.grid(i = 2:5, alpha = alpha_grid) %>%
  mutate(p = p_i(i, alpha),
         i_lbl = paste0("i = ", i))

# Anchor points: (0, 1/i), (0.5, 1), (1, i)
anchors <- data.frame(i = rep(2:5, each = 3),
                      alpha = rep(c(0, 0.5, 1), 4),
                      p = c(1/2, 1, 2,  1/3, 1, 3,  1/4, 1, 4,  1/5, 1, 5)) %>%
  mutate(i_lbl = paste0("i = ", i))

fig1 <- ggplot(pi_df, aes(x = alpha, y = p, colour = i_lbl, group = i_lbl)) +
  geom_line(linewidth = 0.7) +
  geom_point(data = anchors, size = 1.8) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  scale_colour_brewer(palette = "Set1", name = "Index") +
  labs(title = expression("PATP exponent parameterization " * p[i](alpha)),
       x = expression(alpha), y = expression(p[i](alpha))) +
  theme_patp +
  coord_cartesian(ylim = c(0, 5.2))

ggsave(file.path(fig_dir, "fig1_p_i_alpha.pdf"), fig1,
       width = 5.5, height = 3.8, device = cairo_pdf)
cat("Saved fig1_p_i_alpha.pdf\n")

# ===================================================================
# Figure 2: g_2(alpha) for five canonical distributions
# ===================================================================

g2_df <- read.csv("results/theoretical_g2.csv")

# Keep only the five most illustrative distributions.
g2_subset <- g2_df %>%
  filter(distribution %in% c("Laplace", "Gaussian", "Uniform",
                              "GG(0.5)", "GG(4)")) %>%
  mutate(distribution = factor(distribution,
                               levels = c("GG(0.5)", "Laplace", "Gaussian",
                                          "Uniform", "GG(4)")))

fig2 <- ggplot(g2_subset, aes(x = alpha, y = g2,
                              colour = distribution, group = distribution)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0.5, linetype = "dotted", colour = "red", linewidth = 0.3) +
  annotate("text", x = 0.5, y = 0.05, label = "alpha = 1/2 (degenerate)",
           hjust = 0.5, size = 2.8, colour = "red") +
  scale_colour_brewer(palette = "Dark2", name = "Distribution") +
  labs(title = expression("Variance-reduction coefficient " * g[2](alpha) *
                          " for canonical distributions"),
       x = expression(alpha),
       y = expression(g[2](alpha) * " = " * Var[PATP] / Var[OLS])) +
  theme_patp +
  coord_cartesian(ylim = c(0, 1.1))

ggsave(file.path(fig_dir, "fig2_g2_alpha.pdf"), fig2,
       width = 6, height = 3.8, device = cairo_pdf)
cat("Saved fig2_g2_alpha.pdf\n")

# ===================================================================
# Figure 3: topographic plane (kappa, k)
# ===================================================================

# Canonical points from Section 2.4 of the paper.
topo_canonical <- data.frame(
  name = c("Gaussian", "Laplace", "Simpson", "Uniform",
           "Arcsine"),
  gamma_4 = c(0, 3, -0.6, -1.2, -1.5),
  k = c(2.0663, 1.9300, 2.0240, 1.7321, 1.1107)
) %>%
  mutate(kappa = 1 / sqrt(gamma_4 + 3))

# GG(beta) family as a curve.
gg_curve <- data.frame(beta = seq(0.5, 8, length.out = 200)) %>%
  mutate(
    # γ_4 = Γ(5/β)Γ(1/β)/Γ(3/β)² − 3
    gamma_4 = gamma(5 / beta) * gamma(1 / beta) / gamma(3 / beta)^2 - 3,
    # H = log(σ · Γ(1/β) · 2^{1/β}) − log(β) + 1/β + log(2) for GG with var=1
    # but k is invariant; compute numerically via density-based formula
    # H for GG(β) with c chosen so var=1:
    #   c = sqrt(Γ(1/β)/Γ(3/β))
    #   H = 1/β − log(β/(2 c Γ(1/β)))
    c_scale = sqrt(gamma(1 / beta) / gamma(3 / beta)),
    H = 1 / beta - log(beta / (2 * c_scale * gamma(1 / beta))),
    sigma = 1,
    Delta_e = exp(H) / 2,
    k = Delta_e / sigma,
    kappa = 1 / sqrt(gamma_4 + 3)
  )

fig3 <- ggplot() +
  # GG curve
  geom_path(data = gg_curve, aes(x = kappa, y = k),
            colour = "steelblue", linewidth = 0.6, alpha = 0.7) +
  # Highlight β = 0.5, 1, 2, 4 on GG curve
  geom_point(data = gg_curve %>%
               filter(beta %in% c(0.5, 1, 2, 4) | abs(beta - 0.5) < 0.02 |
                      abs(beta - 1) < 0.02 | abs(beta - 2) < 0.02 |
                      abs(beta - 4) < 0.02) %>%
               group_by(beta_round = round(beta * 2) / 2) %>%
               slice(1),
             aes(x = kappa, y = k), size = 1.2, colour = "steelblue") +
  # Canonical points
  geom_point(data = topo_canonical, aes(x = kappa, y = k),
             size = 2.8, colour = "darkred", shape = 16) +
  geom_text(data = topo_canonical, aes(x = kappa, y = k, label = name),
            size = 2.8, hjust = -0.1, vjust = -0.5, colour = "darkred") +
  # k_max horizontal line
  geom_hline(yintercept = 2.0663, linetype = "dashed",
             colour = "grey50", linewidth = 0.3) +
  annotate("text", x = 0.05, y = 2.0663, label = "k_max ≈ 2.066",
           hjust = 0, vjust = -0.4, size = 2.6, colour = "grey40") +
  # Annotations for GG curve
  annotate("text", x = 0.41, y = 1.9, label = "GG(β=1)", size = 2.5,
           colour = "steelblue", hjust = 1.1) +
  annotate("text", x = 0.745, y = 1.7, label = "GG(β→∞)", size = 2.5,
           colour = "steelblue", hjust = -0.05) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 2.3), breaks = seq(0, 2, 0.5)) +
  labs(title = "Topographic plane of symmetric distributions",
       subtitle = "Contrexcess kappa = 1/sqrt(gamma4 + 3); entropy coefficient k = exp(H)/(2 sigma)",
       x = expression(kappa),
       y = expression(k)) +
  theme_patp

ggsave(file.path(fig_dir, "fig3_topographic.pdf"), fig3,
       width = 6.5, height = 4.5, device = cairo_pdf)
cat("Saved fig3_topographic.pdf\n")

# ===================================================================
# Figure 4: PATP basis functions phi_2(xi; alpha) for alpha in {0, 0.25, 0.5, 0.75, 1}
# ===================================================================

xi_grid <- seq(-3, 3, length.out = 401)
alpha_vals <- c(0, 0.25, 0.5, 0.75, 1)

phi2_df <- expand.grid(xi = xi_grid, alpha = alpha_vals) %>%
  mutate(p = p_i(2, alpha),
         phi = sign(xi) * abs(xi)^p,
         alpha_lbl = factor(alpha,
                            levels = alpha_vals,
                            labels = c("alpha=0 (fractal)", "alpha=0.25",
                                       "alpha=0.5 (degenerate)", "alpha=0.75",
                                       "alpha=1 (power)")))

fig4 <- ggplot(phi2_df, aes(x = xi, y = phi,
                            colour = alpha_lbl, group = alpha_lbl)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  scale_colour_viridis_d(name = expression(alpha), option = "plasma") +
  labs(title = expression("PATP basis function " *
                          varphi[2](xi*";"*alpha) ==
                          sign(xi) %.% abs(xi)^{p[2](alpha)}),
       x = expression(xi), y = expression(varphi[2](xi*";"*alpha))) +
  theme_patp +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-9, 9))

ggsave(file.path(fig_dir, "fig4_phi2_basis.pdf"), fig4,
       width = 6, height = 3.8, device = cairo_pdf)
cat("Saved fig4_phi2_basis.pdf\n")

cat("\nAll 4 theoretical figures saved to", fig_dir, "\n")
