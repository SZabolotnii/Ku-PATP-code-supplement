#!/usr/bin/env Rscript
# 04_results_viz.R -- experimental MC-result figures
# PATP paper, §6.2
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

fig_dir <- "../paper/figures"

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
# ===================================================================

fig5_df <- mc_df %>%
  mutate(distribution = factor(distribution,
                               levels = c("Laplace", "GG(1.5)", "GG(4)", "Beta(2,5)")),
         alpha_lbl = factor(sprintf("α = %.2f", alpha)))

fig5 <- ggplot(fig5_df, aes(x = N, y = are,
                            colour = alpha_lbl, group = alpha_lbl)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_wrap(~ distribution, scales = "free_y") +
  scale_x_log10(breaks = c(50, 100, 200, 500),
                labels = c("50", "100", "200", "500")) +
  scale_colour_viridis_d(option = "plasma", name = expression(alpha),
                        end = 0.85) +
  labs(title = expression("Empirical ARE = " * Var[OLS]/Var[PATP] * " over 1000 MC replications"),
       x = "Sample size N",
       y = "ARE (PATP vs OLS)") +
  theme_patp

ggsave(file.path(fig_dir, "fig5_are_vs_N.pdf"), fig5,
       width = 6.5, height = 4.5, device = cairo_pdf)
cat("Saved fig5_are_vs_N.pdf\n")

# ===================================================================
# Figure 6: empirical g_2 vs theoretical g_2 (proxy diagnostic scatter)
# ===================================================================

# Join MC and theoretical values by interpolation.
fig6_validation <- mc_df %>%
  filter(N >= 200) %>%
  rowwise() %>%
  mutate(
    g2_theo = {
      dist_name <- distribution
      sub <- theo_interp_df %>% filter(distribution == .env$dist_name)
      if (nrow(sub) == 0) {
        NA_real_
      } else {
        approx(sub$alpha, sub$g2, xout = .data$alpha, rule = 2)$y
      }
    }
  ) %>%
  ungroup() %>%
  filter(!is.na(g2_theo))

fig6 <- ggplot(fig6_validation, aes(x = g2_theo, y = g2_empirical,
                                     colour = distribution,
                                     shape = factor(N))) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_colour_brewer(palette = "Dark2", name = "Distribution") +
  scale_shape_manual(values = c(`200` = 16, `500` = 17), name = "N") +
  labs(title = expression("Proxy diagnostic for " * g[2](alpha) * ": closed form vs M-estimator"),
       subtitle = "Points near the diagonal would indicate agreement; departures show proxy mismatch; N = 200, 500",
       x = expression(g[2]^{"(theo)"} * "(alpha) from the closed-form moment formula"),
       y = expression(g[2]^{"(emp)"} * "(alpha) from 1000 MC replications")) +
  coord_equal(xlim = c(0, 2.65), ylim = c(0, 2.65)) +
  theme_patp

ggsave(file.path(fig_dir, "fig6_validation.pdf"), fig6,
       width = 6.5, height = 5, device = cairo_pdf)
cat("Saved fig6_validation.pdf\n")

# ===================================================================
# Figure 7 (bonus): bias-variance decomposition
# ===================================================================

bv_df <- mc_df %>%
  pivot_longer(cols = c(var_patp, bias),
               names_to = "metric",
               values_to = "value") %>%
  mutate(
    value = ifelse(metric == "bias", abs(value), value),
    metric = recode(metric,
                    var_patp = "Var[μ̂_PATP]",
                    bias     = "|Bias|"),
    distribution = factor(distribution,
                          levels = c("Laplace", "GG(1.5)", "GG(4)", "Beta(2,5)"))
  )

fig7 <- ggplot(bv_df %>% filter(N == 200), aes(x = alpha, y = value,
                                                colour = metric)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_wrap(~ distribution, scales = "free_y") +
  scale_y_log10(labels = scales::label_scientific(digits = 1)) +
  scale_colour_manual(values = c("Var[μ̂_PATP]" = "#1f78b4",
                                 "|Bias|" = "#e31a1c"),
                      name = "") +
  labs(title = "Bias-variance decomposition for the PATP mean estimator (N = 200)",
       x = expression(alpha),
       y = "Value (log scale)") +
  theme_patp

ggsave(file.path(fig_dir, "fig7_bias_variance.pdf"), fig7,
       width = 6.5, height = 4.5, device = cairo_pdf)
cat("Saved fig7_bias_variance.pdf\n")

cat("\nAll experimental figures saved to", fig_dir, "\n")
