#!/usr/bin/env Rscript
# 03_monte_carlo.R — Monte Carlo симуляції PATP проти OLS
# PATP paper, §6.2 експериментальна валідація
# Дата: 2026-05-11

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

set.seed(2026)

# ===================================================================
# 1. Генератори розподілів (centred, симетричні)
# ===================================================================

# GG(β) з var = 1 та mean = 0
rgg <- function(n, beta) {
  c_scale <- sqrt(gamma(1 / beta) / gamma(3 / beta))
  u <- runif(n)
  # Інверсія через GammaIncomplete: трансформація з Gamma
  # X = sign(U-0.5) * c * (Γ⁻¹(2|U-0.5|, 1/β))^{1/β}
  s <- sign(u - 0.5)
  z <- qgamma(2 * abs(u - 0.5), shape = 1 / beta, rate = 1)
  s * c_scale * z^(1 / beta)
}

# Beta(2, 5), centred and scaled to var = 1
rbeta25 <- function(n, a = 2, b = 5) {
  x <- rbeta(n, a, b)
  mu <- a / (a + b)
  var_b <- a * b / ((a + b)^2 * (a + b + 1))
  (x - mu) / sqrt(var_b)
}

# Test generators
stopifnot(abs(var(rgg(100000, 1.5)) - 1) < 0.05)
stopifnot(abs(var(rgg(100000, 4)) - 1) < 0.05)
stopifnot(abs(var(rbeta25(100000)) - 1) < 0.05)

# ===================================================================
# 2. PATP оцінка через uniroot
# ===================================================================

p_i <- function(i, alpha) {
  1 / i + (4 - i - 3 / i) * alpha + (2 * i - 4 + 2 / i) * alpha^2
}

# Score function: Z(mu; alpha) = (1/N) Σ [h1*(-1) + h2*(-p|xi|^{p-1})] for xi = x - mu
# Equivalent stationarity: dL/dmu = 0
# For S=2, optimal h*(mu) depends on mu via empirical moments;
# practical approach: weight h1 and h2 with theoretical coefficients
# from F_2^{-1} * b at OLS mu_0, then solve fixed-h Newton-Raphson.
# For simplicity, use empirical moment-matching: PATP estimator is fixed-point
# of one-step NR from OLS, weighted by theoretical h*.

# Simpler practical implementation: PATP estimator as M-estimator
# with psi-function psi(xi; alpha) = -1 + h2/h1 * (-p|xi|^{p-1}).
# Since psi(xi) = -1 - k * p * |xi|^{p-1} for some constant k tuned at OLS-residuals,
# we use a one-step correction approach:
#   1. mu_0 = sample mean
#   2. compute residual moments at mu_0
#   3. PATP estimate = mu_0 + adjustment
# For simplicity here, we directly solve via uniroot using psi defined below.

# psi function: derivative of L_S w.r.t. mu, set to zero.
# psi(mu) = -1 - lambda(alpha) * p * mean(|x - mu|^{p-1})
# where lambda depends on h_2*/h_1*; for fixed alpha and centred mu,
# practical PATP estimator solves: mean(psi_PATP(x_n - mu; alpha)) = 0,
# with psi_PATP(xi; alpha) = sign(xi) * |xi|^{p_2(alpha)} (signed-parity basis function).

# So: PATP-2 estimator solves mean(sign(x - mu) * |x - mu|^{p_2(alpha)}) = 0
# This is the M-estimator interpretation.

patp_estimator <- function(x, alpha, mu_init = NULL) {
  if (abs(alpha - 0.5) < 0.01) return(mean(x))  # degenerate -> OLS
  if (is.null(mu_init)) mu_init <- mean(x)
  p <- p_i(2, alpha)
  if (p <= 0) return(NA_real_)

  # Score: ψ(μ) = mean(sign(x - μ) * |x - μ|^p)
  score <- function(mu) {
    xi <- x - mu
    mean(sign(xi) * abs(xi)^p)
  }
  # Search bracket
  s_lo <- mu_init - 5 * sd(x)
  s_hi <- mu_init + 5 * sd(x)
  v_lo <- score(s_lo)
  v_hi <- score(s_hi)
  if (is.na(v_lo) || is.na(v_hi)) return(NA_real_)
  if (sign(v_lo) == sign(v_hi)) return(mu_init)
  tryCatch(
    uniroot(score, lower = s_lo, upper = s_hi, tol = 1e-8)$root,
    error = function(e) mu_init
  )
}

# ===================================================================
# 2A. Robust baseline estimators
# ===================================================================

winsor_mean <- function(x, trim = 0.1) {
  qs <- quantile(x, probs = c(trim, 1 - trim), names = FALSE, type = 8)
  mean(pmin(pmax(x, qs[1]), qs[2]))
}

huber_location <- function(x, k = 1.345, tol = 1e-8, max_iter = 50) {
  mu <- median(x)
  s <- mad(x, constant = 1.4826)
  if (!is.finite(s) || s <= 1e-12) s <- sd(x)
  if (!is.finite(s) || s <= 1e-12) return(mu)

  for (iter in seq_len(max_iter)) {
    r <- (x - mu) / s
    w <- pmin(1, k / pmax(abs(r), 1e-12))
    mu_new <- sum(w * x) / sum(w)
    if (!is.finite(mu_new)) return(mu)
    if (abs(mu_new - mu) < tol * max(1, abs(mu))) return(mu_new)
    mu <- mu_new
  }
  mu
}

median_of_means <- function(x, blocks = NULL) {
  n <- length(x)
  if (is.null(blocks)) blocks <- max(3, floor(sqrt(n)))
  blocks <- min(blocks, n)
  idx <- split(seq_len(n), cut(seq_len(n), breaks = blocks, labels = FALSE))
  median(vapply(idx, function(ii) mean(x[ii]), numeric(1)))
}

estimate_baseline <- function(x, estimator) {
  switch(estimator,
         mean = mean(x),
         median = median(x),
         trimmed_mean = mean(x, trim = 0.1),
         winsor_mean = winsor_mean(x, trim = 0.1),
         huber = huber_location(x),
         median_of_means = median_of_means(x),
         stop("Unknown estimator: ", estimator))
}

# ===================================================================
# 3. Monte Carlo sweep
# ===================================================================

distributions <- list(
  Laplace      = function(n) sqrt(0.5) * rexp(n) * sample(c(-1, 1), n, TRUE),
  `GG(1.5)`    = function(n) rgg(n, 1.5),
  `GG(4)`      = function(n) rgg(n, 4),
  `Beta(2,5)`  = function(n) rbeta25(n)
)

alpha_grid <- c(0.05, 0.30, 0.70, 0.95)
N_grid <- c(50, 100, 200, 500)
M <- 1000  # MC реплік

mc_results <- list()
iter <- 0
total_iter <- length(distributions) * length(N_grid) * (length(alpha_grid) + 1)

for (dist_name in names(distributions)) {
  rgen <- distributions[[dist_name]]
  for (N in N_grid) {
    # OLS reference
    iter <- iter + 1
    ols_means <- replicate(M, mean(rgen(N)))
    var_ols <- var(ols_means)

    for (alpha in alpha_grid) {
      iter <- iter + 1
      patp_means <- replicate(M, {
        x <- rgen(N)
        patp_estimator(x, alpha, mu_init = mean(x))
      })
      patp_means <- patp_means[!is.na(patp_means)]
      var_patp <- var(patp_means)
      bias_patp <- mean(patp_means)  # true mu = 0
      mse_patp <- var_patp + bias_patp^2

      mc_results[[length(mc_results) + 1]] <- data.frame(
        distribution = dist_name,
        N = N,
        alpha = alpha,
        var_ols = var_ols,
        var_patp = var_patp,
        bias = bias_patp,
        mse = mse_patp,
        are = var_ols / var_patp,
        g2_empirical = var_patp / var_ols
      )
    }
    cat(sprintf("[%d/%d] %s, N=%d done\n",
                iter, total_iter, dist_name, N))
  }
}

mc_df <- bind_rows(mc_results)

# ===================================================================
# 3A. Robust baselines on the same distributions
# ===================================================================

baseline_estimators <- c("mean", "median", "trimmed_mean", "winsor_mean",
                         "huber", "median_of_means")

baseline_results <- list()
for (dist_name in names(distributions)) {
  rgen <- distributions[[dist_name]]
  for (N in N_grid) {
    for (estimator in baseline_estimators) {
      estimates <- replicate(M, estimate_baseline(rgen(N), estimator))
      baseline_results[[length(baseline_results) + 1]] <- data.frame(
        distribution = dist_name,
        N = N,
        estimator = estimator,
        variance = var(estimates),
        bias = mean(estimates),
        mse = var(estimates) + mean(estimates)^2
      )
    }
  }
}

baseline_df <- bind_rows(baseline_results) %>%
  group_by(distribution, N) %>%
  mutate(relative_mse_vs_mean = mse / mse[estimator == "mean"]) %>%
  ungroup()

baseline_summary <- baseline_df %>%
  filter(N %in% c(100, 500)) %>%
  mutate(across(c(variance, bias, mse, relative_mse_vs_mean), \(x) round(x, 5))) %>%
  arrange(distribution, N, relative_mse_vs_mean)

# ===================================================================
# 3B. Alpha ablation and lightweight runtime diagnostics
# ===================================================================

alpha_ablation <- mc_df %>%
  group_by(distribution, N) %>%
  summarise(
    alpha_best_mse = alpha[which.min(mse)],
    mse_best = min(mse, na.rm = TRUE),
    alpha_best_are = alpha[which.max(are)],
    are_best = max(are, na.rm = TRUE),
    alpha_worst_mse = alpha[which.max(mse)],
    mse_worst = max(mse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(mse_best, are_best, mse_worst), \(x) round(x, 5)))

time_call <- function(fun, reps = 50) {
  invisible(fun())  # warm-up
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(reps)) invisible(fun())
  (proc.time()[["elapsed"]] - t0) / reps
}

runtime_N <- c(100, 1000, 10000, 100000)
runtime_results <- list()
runtime_estimators <- c("mean", "median", "huber", "median_of_means",
                        "patp_alpha_0.05", "patp_alpha_0.95")
for (N in runtime_N) {
  x <- distributions$Laplace(N)
  reps <- if (N <= 100) 1000 else if (N <= 1000) 500 else if (N <= 10000) 100 else 20
  for (estimator in runtime_estimators) {
    elapsed <- switch(estimator,
      mean = time_call(function() mean(x), reps = reps),
      median = time_call(function() median(x), reps = reps),
      huber = time_call(function() huber_location(x), reps = reps),
      median_of_means = time_call(function() median_of_means(x), reps = reps),
      patp_alpha_0.05 = time_call(function() patp_estimator(x, 0.05, mu_init = mean(x)), reps = reps),
      patp_alpha_0.95 = time_call(function() patp_estimator(x, 0.95, mu_init = mean(x)), reps = reps)
    )
    runtime_results[[length(runtime_results) + 1]] <- data.frame(
      distribution = "Laplace",
      N = N,
      estimator = estimator,
      reps = reps,
      elapsed_seconds = elapsed
    )
  }
}
runtime_df <- bind_rows(runtime_results)

# ===================================================================
# 4. Збереження результатів
# ===================================================================

dir.create("results", showWarnings = FALSE)
write.csv(mc_df, "results/monte_carlo.csv", row.names = FALSE)
write.csv(baseline_df, "results/robust_baselines.csv", row.names = FALSE)
write.csv(baseline_summary, "results/robust_baselines_summary.csv", row.names = FALSE)
write.csv(alpha_ablation, "results/alpha_ablation_summary.csv", row.names = FALSE)
write.csv(runtime_df, "results/runtime_summary.csv", row.names = FALSE)

# Зведена таблиця
summary_mc <- mc_df %>%
  filter(N %in% c(100, 500)) %>%
  mutate(across(c(var_patp, bias, mse, are, g2_empirical), \(x) round(x, 4))) %>%
  arrange(distribution, N, alpha)

cat("\n=== Підсумок MC (N = 100, 500) ===\n")
print(summary_mc)

cat("\n=== Robust baseline summary (N = 100, 500) ===\n")
print(baseline_summary)

cat("\n=== Alpha ablation summary ===\n")
print(alpha_ablation)

cat("\n=== Runtime summary (Laplace) ===\n")
print(runtime_df)

write.csv(summary_mc, "results/monte_carlo_summary.csv", row.names = FALSE)
cat("\nDone. CSV збережено у R/results/\n")
