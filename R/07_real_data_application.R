#!/usr/bin/env Rscript
# 07_real_data_application.R
# Real-data demonstration of the PATP *location* estimator (AJS revision).
#
# We use the fully validated full F_2^{-1}b estimator of 05_full_patp_estimator.R
# (patp_full) -- the estimator whose finite-sample g_2(alpha) was shown in
# 03_monte_carlo.R to converge to the closed-form theory. The task is robust
# location estimation of a heavy-tailed, near-symmetric real sample.
#
# Data: EuStockMarkets (base R) daily log-returns. Equity-index returns are
# leptokurtic and approximately symmetric -> the regime where the symmetric
# g_2(alpha) applies. We pick the series that is most symmetric and leptokurtic.
#
# Two numerical points, both principled:
#  (1) The sample is STANDARDISED to unit robust scale (z = (x-med)/MAD) before
#      estimation. g_2 is scale-invariant and location estimation is scale-
#      equivariant, so this does not change the target; it places the data in
#      the unit-scale regime in which patp_full was validated, avoiding the
#      scale-induced fractal-weight instability seen on raw returns.
#  (2) alpha* is chosen by the robust grid + bootstrap criterion of Sec.5
#      (minimum realised estimator variance), NOT the plug-in g_2 formula,
#      which is unstable on empirical fractional moments.
#
# Output: realised efficiency g_2_hat(alpha) = Var(PATP)/Var(mean) over the
# grid, the selected alpha*, the realised variance reduction, a comparison with
# the theoretical g_2 of a kurtosis-matched Student-t, and robust baselines.

suppressPackageStartupMessages({ library(dplyr) })

.this_dir <- if (file.exists("05_full_patp_estimator.R")) "." else "R"
source(file.path(.this_dir, "05_full_patp_estimator.R"))  # patp_full, p_i, empirical_moments
source(file.path(.this_dir, "01_theoretical_g2.R"))       # g2_alpha, nu_q_t

set.seed(2026)

.skew   <- function(x) { m <- mean(x); s <- sqrt(mean((x - m)^2)); mean((x - m)^3) / s^3 }
.exkurt <- function(x) { m <- mean(x); v <- mean((x - m)^2); mean((x - m)^4) / v^2 - 3 }

# Simple Huber location (robust baseline; same convention as 03_monte_carlo.R)
huber_loc <- function(x, k = 1.345, tol = 1e-8, max_iter = 50) {
  mu <- median(x); s <- mad(x); if (!is.finite(s) || s <= 1e-12) s <- sd(x)
  for (i in seq_len(max_iter)) {
    r <- (x - mu) / s; w <- pmin(1, k / pmax(abs(r), 1e-12))
    mn <- sum(w * x) / sum(w)
    if (abs(mn - mu) < tol * max(1, abs(mu))) return(mn)
    mu <- mn
  }
  mu
}

# ====================================================================
# 1. Data: choose the most symmetric, leptokurtic return series
# ====================================================================
data(EuStockMarkets, package = "datasets")
logret <- apply(as.matrix(EuStockMarkets), 2, function(p) diff(log(p)))
nms <- colnames(logret)

shape <- data.frame(series = nms,
                    n = nrow(logret),
                    skew = apply(logret, 2, .skew),
                    exkurt = apply(logret, 2, .exkurt))
shape$abs_skew <- abs(shape$skew)
shape <- shape %>% arrange(abs_skew)
cat("=== EuStockMarkets daily log-return shape ===\n")
print(shape[, c("series", "n", "skew", "exkurt")], row.names = FALSE, digits = 3)

sel <- shape$series[which.min(shape$abs_skew)]
x   <- logret[, sel]
n   <- length(x)
cat(sprintf("\nSelected series: %s  (skew=%.3f, excess kurtosis=%.2f, n=%d)\n",
            sel, .skew(x), .exkurt(x), n))

# standardise to unit robust scale (validated regime; g_2 is scale-invariant)
z  <- (x - median(x)) / mad(x)
ek <- .exkurt(z)

# ====================================================================
# 2. Realised efficiency g_2_hat(alpha) by bootstrap; alpha* = argmin variance
# ====================================================================
ag <- setdiff(seq(0.05, 0.95, 0.05), 0.5)
B  <- 2000

est_mean  <- numeric(B)
est_med   <- numeric(B)
est_huber <- numeric(B)
est_patp  <- matrix(NA_real_, B, length(ag))
for (b in seq_len(B)) {
  zb <- z[sample.int(n, n, replace = TRUE)]
  m0 <- mean(zb)
  est_mean[b]  <- m0
  est_med[b]   <- median(zb)
  est_huber[b] <- huber_loc(zb)
  for (k in seq_along(ag)) est_patp[b, k] <- patp_full(zb, ag[k], mu_init = m0)
}

v_mean   <- var(est_mean)
g2_real  <- apply(est_patp, 2, function(v) var(v[is.finite(v)])) / v_mean
k_star   <- which.min(g2_real)
alpha_st <- ag[k_star]
g2_st    <- g2_real[k_star]
g2_med   <- var(est_med)   / v_mean
g2_huber <- var(est_huber) / v_mean

# theoretical g_2 of a kurtosis-matched Student-t (ex.kurtosis = 6/(nu-4))
nu_fit <- 4 + 6 / ek
g2_t_theory <- if (nu_fit > 4) g2_alpha(alpha_st, 1, function(q) nu_q_t(q, nu_fit),
                                        function(q) 0) else NA_real_

cat(sprintf("\nStandardised sample: excess kurtosis = %.2f -> kurtosis-matched t(nu=%.1f)\n",
            ek, nu_fit))
cat("\n=== Realised efficiency g_2_hat(alpha) = Var(PATP)/Var(mean), bootstrap B=2000 ===\n")
print(data.frame(alpha = ag, g2_realized = round(g2_real, 4)), row.names = FALSE)
cat(sprintf("\nalpha* (argmin realised variance) = %.2f\n", alpha_st))
cat(sprintf("realised g_2(alpha*)              = %.4f  (%.1f%% variance reduction vs sample mean)\n",
            g2_st, 100 * (1 - g2_st)))
cat(sprintf("theoretical g_2(alpha*) for matched t(%.1f) = %.4f\n", nu_fit, g2_t_theory))
cat(sprintf("baselines vs mean: median g_2 = %.4f,  Huber g_2 = %.4f\n", g2_med, g2_huber))

# ====================================================================
# 3. Save
# ====================================================================
dir.create("results", showWarnings = FALSE)
write.csv(data.frame(alpha = ag, g2_realized = g2_real),
          "results/real_data_g2_curve.csv", row.names = FALSE)
write.csv(shape[, c("series", "n", "skew", "exkurt")],
          "results/real_data_series.csv", row.names = FALSE)
out <- data.frame(
  dataset = "EuStockMarkets daily log-returns",
  series = sel, n = n,
  skew = round(.skew(x), 4), exkurt = round(.exkurt(x), 3),
  alpha_star = alpha_st,
  g2_realized = round(g2_st, 4),
  var_reduction_pct = round(100 * (1 - g2_st), 1),
  nu_matched_t = round(nu_fit, 2),
  g2_theory_matched_t = round(g2_t_theory, 4),
  g2_median = round(g2_med, 4),
  g2_huber = round(g2_huber, 4)
)
write.csv(out, "results/real_data_application.csv", row.names = FALSE)
cat("\nSaved results/real_data_application.csv, real_data_g2_curve.csv, real_data_series.csv\n")
