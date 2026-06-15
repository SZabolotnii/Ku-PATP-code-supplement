#!/usr/bin/env Rscript
# 04_results_viz.R -- experimental MC-result figures
# PATP paper, §6.2
# Date: 2026-05-11; rewritten 2026-06-15 (AJS revision):
#   - headline estimator switched from scalar proxy to full F_2^{-1}b
#   - distributions restricted to symmetric laws (Beta(2,5) -> Student-t(6))
#   - Fig 6 is now a CONVERGENCE validation (empirical g_2 -> closed-form theory)

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

# Symmetric illustration set (D6v2).
sym_levels <- c("Laplace", "GG(1.5)", "GG(4)", "Student-t(6)")

# Load data.
mc_df <- read.csv("results/monte_carlo.csv")
theo_df <- read.csv("results/theoretical_g2.csv")
theo_interp_df <- theo_df %>%
  filter(!is.na(g2)) %>%
  group_by(distribution, alpha) %>%
  summarise(g2 = mean(g2), .groups = "drop") %>%
  arrange(distribution, alpha)

# ===================================================================
# Figure 5: ARE vs N across distributions and alpha values
# (full F_2^{-1}b estimator)
# ===================================================================

fig5_df <- mc_df %>%
  mutate(distribution = factor(distribution, levels = sym_levels),
         alpha_lbl = factor(sprintf("α = %.2f", alpha)))

fig5 <- ggplot(fig5_df, aes(x = N, y = are_full,
                            colour = alpha_lbl, group = alpha_lbl)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_wrap(~ distribution, scales = "free_y") +
  scale_x_log10(breaks = c(50, 100, 200, 500),
                labels = c("50", "100", "200", "500")) +
  scale_colour_viridis_d(option = "plasma", name = expression(alpha),
                        end = 0.85) +
  labs(title = expression("Empirical ARE = " * Var[OLS] / Var[hat(mu)[PATP]] *
                          " (full " * F[2]^{-1} * bold(b) * " estimator, 1000 MC replications)"),
       x = "Sample size N",
       y = "ARE (PATP vs OLS)") +
  theme_patp

ggsave(file.path(fig_dir, "fig5_are_vs_N.pdf"), fig5,
       width = 6.5, height = 4.5, device = cairo_pdf)
cat("Saved fig5_are_vs_N.pdf\n")

# ===================================================================
# Figure 6: convergence of the full F_2^{-1}b estimator's empirical g_2 to the
# closed-form theory g_2(alpha). This is the DIRECT validation of the paper's
# central result (replaces the former proxy-diagnostic scatter).
# ===================================================================

conv_df <- read.csv("results/convergence_g2.csv")

# Theoretical g_2 at the two convergence endpoints, per distribution.
conv_theo <- conv_df %>%
  distinct(distribution, alpha) %>%
  rowwise() %>%
  mutate(g2_theo = {
    sub <- theo_interp_df %>% filter(distribution == .env$distribution)
    if (nrow(sub) == 0) NA_real_
    else approx(sub$alpha, sub$g2, xout = .data$alpha, rule = 2)$y
  }) %>%
  ungroup()

conv_plot <- conv_df %>%
  mutate(distribution = factor(distribution, levels = sym_levels),
         alpha_lbl = factor(sprintf("α = %.2f", alpha)))
conv_theo_plot <- conv_theo %>%
  mutate(distribution = factor(distribution, levels = sym_levels),
         alpha_lbl = factor(sprintf("α = %.2f", alpha)))

fig6 <- ggplot(conv_plot, aes(x = N, y = g2_empirical_full, colour = alpha_lbl)) +
  geom_hline(data = conv_theo_plot,
             aes(yintercept = g2_theo, colour = alpha_lbl),
             linetype = "dashed", linewidth = 0.45) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.9) +
  facet_wrap(~ distribution, scales = "free_y") +
  scale_x_log10(breaks = c(100, 250, 500, 1000, 2000, 4000)) +
  scale_colour_viridis_d(option = "plasma", end = 0.8, name = expression(alpha)) +
  labs(title = expression("Convergence of the full " * F[2]^{-1} * bold(b) *
                          " estimator's " * hat(g)[2](alpha) * " to closed-form theory"),
       subtitle = "Dashed: theoretical g_2(alpha) from the closed-form moment formula. Symmetric laws; M = 2000 replications.",
       x = "Sample size N (log scale)",
       y = expression(hat(g)[2](alpha) * " = " * Var[hat(mu)[PATP]] / Var[hat(mu)[OLS]])) +
  theme_patp

ggsave(file.path(fig_dir, "fig6_validation.pdf"), fig6,
       width = 6.5, height = 5, device = cairo_pdf)
cat("Saved fig6_validation.pdf\n")

# ===================================================================
# Figure 7: bias-variance decomposition (full F_2^{-1}b estimator)
# ===================================================================

bv_df <- mc_df %>%
  pivot_longer(cols = c(var_full, bias_full),
               names_to = "metric",
               values_to = "value") %>%
  mutate(
    value = ifelse(metric == "bias_full", abs(value), value),
    metric = recode(metric,
                    var_full  = "Var[mu_PATP]",
                    bias_full = "|Bias|"),
    distribution = factor(distribution, levels = sym_levels)
  )

fig7 <- ggplot(bv_df %>% filter(N == 200), aes(x = alpha, y = value,
                                                colour = metric)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_wrap(~ distribution, scales = "free_y") +
  scale_y_log10(labels = scales::label_scientific(digits = 1)) +
  scale_colour_manual(values = c("Var[mu_PATP]" = "#1f78b4",
                                 "|Bias|" = "#e31a1c"),
                      name = "") +
  labs(title = "Bias-variance decomposition for the PATP mean estimator (N = 200)",
       x = expression(alpha),
       y = "Value (log scale)") +
  theme_patp

ggsave(file.path(fig_dir, "fig7_bias_variance.pdf"), fig7,
       width = 6.5, height = 4.5, device = cairo_pdf)
cat("Saved fig7_bias_variance.pdf\n")

# ===================================================================
# Figure 8: scalar proxy vs full F_2^{-1}b vs closed-form theory (N = 500).
# Motivates the switch: the naive scalar M-estimator proxy departs from theory
# (badly at the signed-parity end), the full normal-equation estimator tracks it.
# ===================================================================

fig8_df <- mc_df %>%
  filter(N == 500) %>%
  select(distribution, alpha,
         g2_proxy = g2_empirical,
         g2_full  = g2_empirical_full) %>%
  pivot_longer(cols = c(g2_proxy, g2_full),
               names_to = "estimator", values_to = "g2_emp") %>%
  mutate(estimator = recode(estimator,
                            g2_proxy = "Scalar M-estimator (proxy)",
                            g2_full  = "Full F^{-1}b normal equations"),
         distribution = factor(distribution, levels = sym_levels))

fig8_theo <- theo_interp_df %>%
  filter(distribution %in% sym_levels) %>%
  mutate(distribution = factor(distribution, levels = sym_levels))

fig8 <- ggplot() +
  geom_line(data = fig8_theo,
            aes(x = alpha, y = g2),
            colour = "grey25", linewidth = 0.4) +
  geom_point(data = fig8_df,
             aes(x = alpha, y = g2_emp, colour = estimator, shape = estimator),
             size = 2.5) +
  facet_wrap(~ distribution, scales = "free_y") +
  scale_colour_manual(values = c("Scalar M-estimator (proxy)" = "#e31a1c",
                                  "Full F^{-1}b normal equations" = "#1f78b4"),
                      name = "Estimator") +
  scale_shape_manual(values = c("Scalar M-estimator (proxy)" = 17,
                                 "Full F^{-1}b normal equations" = 16),
                     name = "Estimator") +
  labs(title = expression("Empirical " * g[2](alpha) * " under proxy vs full PMM normal equations (N = 500)"),
       subtitle = "Solid grey line: closed-form theoretical g_2(alpha). The full estimator tracks theory; the proxy does not.",
       x = expression(alpha),
       y = expression(hat(g)[2](alpha) * " = " * Var[hat(mu)[PATP]] * " / " * Var[hat(mu)[OLS]])) +
  theme_patp +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "fig8_proxy_vs_full.pdf"), fig8,
       width = 6.8, height = 5.2, device = cairo_pdf)
cat("Saved fig8_proxy_vs_full.pdf\n")

cat("\nAll experimental figures saved to", fig_dir, "\n")
