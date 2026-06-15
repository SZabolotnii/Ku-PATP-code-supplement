#!/usr/bin/env Rscript
# 01_theoretical_g2.R — Теоретичні значення g_2(α) для канонічних розподілів
# PATP paper, §6.1 ілюстрації
# Дата: 2026-05-11

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ====================================================================
# 1. PATP exponent p_2(α) = 1/2 + (1/2)α + α²
# ====================================================================

p_i <- function(i, alpha) {
  1 / i + (4 - i - 3 / i) * alpha + (2 * i - 4 + 2 / i) * alpha^2
}

# ====================================================================
# 2. Дробні моменти ν_q = E[|X|^q] для канонічних розподілів
# (центровані, при необхідності — стандартизовані)
# ====================================================================

# Laplace(0, 1): f(x) = (1/2) exp(-|x|); var = 2
nu_q_laplace <- function(q) {
  if (q <= -1) return(NA_real_)
  gamma(q + 1)  # = ∫_0^∞ x^q e^{-x} dx
}
sigma_q_laplace <- function(q) 0  # симетричний

# Gaussian(0, 1): f(x) = φ(x); var = 1
nu_q_gauss <- function(q) {
  if (q <= -1) return(NA_real_)
  # E[|X|^q] for X ~ N(0,1):  (2/π)^{1/2} * 2^{q/2} * Γ((q+1)/2)
  (2 / pi)^(1/2) * 2^(q/2) * gamma((q + 1) / 2)
}
sigma_q_gauss <- function(q) 0  # симетричний

# Uniform[-1, 1]: density 1/2; var = 1/3
nu_q_uniform <- function(q) {
  if (q <= -1) return(NA_real_)
  1 / (q + 1)  # = ∫_0^1 x^q dx
}
sigma_q_uniform <- function(q) 0  # симетричний

# GG(β) standardised to var = 1; density f(x) ∝ exp(-|x/c|^β), c = sqrt(Γ(1/β)/Γ(3/β))
nu_q_gg <- function(q, beta) {
  if (q <= -1) return(NA_real_)
  c_scale <- sqrt(gamma(1/beta) / gamma(3/beta))
  c_scale^q * gamma((q + 1) / beta) / gamma(1 / beta)
}
sigma_q_gg <- function(q, beta) 0

# Student-t(ν) standardised to var = 1; symmetric, finite moments for q < ν
# (ν > 4 keeps ν_4 finite at α=1). With Y = X·sqrt((ν-2)/ν), X ~ t(ν):
#   E|Y|^q = ((ν-2)/ν)^{q/2}·ν^{q/2}·Γ((q+1)/2)Γ((ν-q)/2)/(Γ(1/2)Γ(ν/2)),  −1<q<ν.
nu_q_t <- function(q, nu) {
  if (q <= -1 || q >= nu) return(NA_real_)
  s <- sqrt((nu - 2) / nu)
  EXq <- nu^(q / 2) * gamma((q + 1) / 2) * gamma((nu - q) / 2) /
         (gamma(0.5) * gamma(nu / 2))
  s^q * EXq
}
sigma_q_t <- function(q, nu) 0  # симетричний

# Beta(2, 5) shifted to mean 0, scaled to var = 1 (numerical)
beta_density_centred <- function(x, a = 2, b = 5) {
  # X ~ Beta(a, b); E[X] = a/(a+b); Var = ab/((a+b)^2(a+b+1))
  mu <- a / (a + b)
  var_b <- a * b / ((a + b)^2 * (a + b + 1))
  sd_b <- sqrt(var_b)
  # Y = (X - mu)/sd has var 1
  y_to_x <- function(y) y * sd_b + mu
  # Density f_Y(y) = f_X(y_to_x(y)) * sd_b
  z <- y_to_x(x)
  ifelse(z > 0 & z < 1, dbeta(z, a, b) * sd_b, 0)
}

nu_q_beta25 <- function(q) {
  # Numerical integration; support is approximately (-0.84, 2.36) standardised
  integrand <- function(x) abs(x)^q * beta_density_centred(x)
  tryCatch(
    integrate(integrand, lower = -1, upper = 3,
              rel.tol = 1e-10, subdivisions = 1000)$value,
    error = function(e) NA_real_
  )
}

sigma_q_beta25 <- function(q) {
  integrand <- function(x) sign(x) * abs(x)^q * beta_density_centred(x)
  tryCatch(
    integrate(integrand, lower = -1, upper = 3,
              rel.tol = 1e-10, subdivisions = 1000)$value,
    error = function(e) NA_real_
  )
}

# ====================================================================
# 3. g_2(α) формула з §4.2 paper
# ====================================================================

g2_alpha <- function(alpha, c_2, nu_q_fn, sigma_q_fn) {
  if (abs(alpha - 0.5) < 0.005) return(1.0)  # вироджена точка
  p <- p_i(2, alpha)
  if (p - 1 <= -1) return(NA_real_)  # неіснування ν_{p-1}

  nu_pm1 <- nu_q_fn(p - 1)
  nu_pp1 <- nu_q_fn(p + 1)
  nu_2p  <- nu_q_fn(2 * p)
  sig_p  <- sigma_q_fn(p)

  if (any(is.na(c(nu_pm1, nu_pp1, nu_2p, sig_p)))) return(NA_real_)

  F22 <- nu_2p - sig_p^2
  num <- c_2 * F22 - nu_pp1^2
  den <- c_2 * (F22 - 2 * p * nu_pp1 * nu_pm1 + p^2 * c_2 * nu_pm1^2)

  if (abs(den) < 1e-12) return(NA_real_)
  num / den
}

# ====================================================================
# 4. Обчислення на сітці α
# ====================================================================

# Include exact endpoints alpha in {0, 1}: there is NO singularity there (only
# alpha = 1/2 is degenerate), so the summary can report the exact g2(0), g2(1)
# that match Sec. 4.2 analytically — e.g. Laplace g2(1) = 3/4 exactly, not the
# old 0.98-grid proxy 0.7438.
alpha_grid <- sort(unique(c(0, 1, seq(0.02, 0.98, by = 0.02))))
alpha_grid <- alpha_grid[abs(alpha_grid - 0.5) > 0.03]  # exclude near-degeneracy

distributions <- list(
  list(name = "Laplace",      c_2 = 2,
       nu_q_fn = nu_q_laplace,            sigma_q_fn = sigma_q_laplace),
  list(name = "Gaussian",     c_2 = 1,
       nu_q_fn = nu_q_gauss,              sigma_q_fn = sigma_q_gauss),
  list(name = "Uniform",      c_2 = 1/3,
       nu_q_fn = nu_q_uniform,            sigma_q_fn = sigma_q_uniform),
  list(name = "GG(0.5)",      c_2 = 1,
       nu_q_fn = function(q) nu_q_gg(q, 0.5), sigma_q_fn = function(q) 0),
  list(name = "GG(1.5)",      c_2 = 1,
       nu_q_fn = function(q) nu_q_gg(q, 1.5), sigma_q_fn = function(q) 0),
  list(name = "GG(4)",        c_2 = 1,
       nu_q_fn = function(q) nu_q_gg(q, 4), sigma_q_fn = function(q) 0),
  list(name = "Student-t(6)", c_2 = 1,
       nu_q_fn = function(q) nu_q_t(q, 6), sigma_q_fn = function(q) 0)
)

results <- list()
for (dst in distributions) {
  g2_vals <- sapply(alpha_grid, function(a) {
    g2_alpha(a, dst$c_2, dst$nu_q_fn, dst$sigma_q_fn)
  })
  results[[dst$name]] <- data.frame(
    distribution = dst$name,
    alpha = alpha_grid,
    g2 = g2_vals
  )
}

theoretical_g2 <- bind_rows(results)

# ====================================================================
# 5. Збереження
# ====================================================================

dir.create("results", showWarnings = FALSE, recursive = TRUE)
write.csv(theoretical_g2, "results/theoretical_g2.csv", row.names = FALSE)

# Підсумкова таблиця у трьох точках: α = 0, α = 1, оптимум
summary_tbl <- theoretical_g2 %>%
  group_by(distribution) %>%
  summarise(
    g2_at_0       = round(g2[which.min(abs(alpha - 0))], 4),   # exact fractal endpoint
    g2_at_1       = round(g2[which.min(abs(alpha - 1))], 4),   # exact signed-parity endpoint
    g2_min        = round(min(g2, na.rm = TRUE), 4),
    alpha_optimal = round(alpha[which.min(g2)], 3),
    .groups = "drop"
  )

cat("\n=== Підсумок g_2(α) для канонічних розподілів ===\n")
print(summary_tbl)
write.csv(summary_tbl, "results/theoretical_g2_summary.csv", row.names = FALSE)

cat("\nDone. CSV збережено в R/results/\n")
