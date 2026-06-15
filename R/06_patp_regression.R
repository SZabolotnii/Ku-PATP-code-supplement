#!/usr/bin/env Rscript
# 06_patp_regression.R
# Form-B PATP estimator for linear regression y = X beta + xi (S = 2).
# Generalises the location estimator of 05_full_patp_estimator.R from the mean
# to regression coefficients. For a linear model with i.i.d. symmetric errors,
# the efficient PATP score
#       psi*(r) = h1* * r + h2* * sign(r)|r|^p,     h* = F_2^{-1} b
# (with F_2, b built from residual moments) yields a one-step Newey estimator
# whose coefficient covariance is g_2(alpha) * sigma^2 (X'X)^{-1}, i.e. the SAME
# variance-reduction factor g_2(alpha) derived for the location case multiplies
# every regression coefficient. This script (a) defines the estimator and
# (b) validates by Monte Carlo that the empirical slope-variance ratio matches
# the closed-form g_2(alpha) for symmetric error laws (D6v2 set).

suppressPackageStartupMessages({
  library(dplyr)
})

.this_dir <- if (file.exists("05_full_patp_estimator.R")) "." else "R"
source(file.path(.this_dir, "05_full_patp_estimator.R"))   # p_i, empirical_moments, build_F2_b_hstar
source(file.path(.this_dir, "01_theoretical_g2.R"))        # g2_alpha, nu_q_*, distributions (theory)

# ====================================================================
# 1. PATP regression estimator (one-step Newey, optionally re-iterated)
# ====================================================================

#' Fit y = X beta + xi by the Form-B PATP normal-equation estimator.
#' @param X numeric design matrix (n x k), include an intercept column if wanted.
#' @param y numeric response (length n).
#' @param alpha PATP parameter in [0,1] (alpha = 1/2 -> OLS).
#' @return list(beta, iter, converged, cond_F2, fallback)
patp_regression <- function(X, y, alpha, max_iter = 10, tol = 1e-7,
                            cond_max = 1e10, eps = 0) {
  X <- as.matrix(X)
  n <- nrow(X); k <- ncol(X)
  XtX <- crossprod(X)                       # X'X
  XtX_inv <- tryCatch(solve(XtX), error = function(e) NULL)
  if (is.null(XtX_inv)) stop("X'X singular")
  beta <- as.vector(XtX_inv %*% crossprod(X, y))   # OLS start
  s0   <- sqrt(mean((y - X %*% beta)^2))           # OLS residual scale
  cap  <- 5 * s0 * sqrt(diag(XtX_inv))             # per-coef step cap (~5 OLS SEs)

  if (abs(alpha - 0.5) < 0.01) {
    return(list(beta = beta, iter = 0, converged = TRUE,
                fallback = "OLS_degenerate", cond_F2 = NA_real_))
  }
  p <- p_i(2, alpha)
  if (p <= 0) return(list(beta = beta, iter = 0, converged = FALSE,
                          fallback = "invalid_p", cond_F2 = NA_real_))

  cond_last <- NA_real_
  for (it in seq_len(max_iter)) {
    r <- as.vector(y - X %*% beta)
    fbh <- build_F2_b_hstar(r, 0, p, cond_max = cond_max, eps = eps)  # moments of residuals
    cond_last <- fbh$cond_F2
    if (!fbh$ok) {
      return(list(beta = beta, iter = it, converged = TRUE,
                  fallback = "singular_F2_to_OLS", cond_F2 = cond_last))
    }
    h <- fbh$h_star
    w_r <- if (eps > 0) sqrt(r^2 + eps^2) else abs(r)        # Sec.5 smoothing safeguard
    psi_star <- h[1] * r + h[2] * sign(r) * w_r^p
    Eprime   <- sum(h * fbh$b)               # E[psi*'] = h1 + h2 * p * nu_{p-1} = h . b
    if (!is.finite(Eprime) || abs(Eprime) < 1e-14) break
    delta <- as.vector(XtX_inv %*% crossprod(X, psi_star)) / Eprime
    delta <- pmax(pmin(delta, cap), -cap)          # damp to prevent blow-ups (cf. patp_full max_step)
    beta_new <- beta + delta
    if (max(abs(delta)) < tol * max(1, max(abs(beta)))) {
      return(list(beta = beta_new, iter = it, converged = TRUE,
                  fallback = "none", cond_F2 = cond_last))
    }
    beta <- beta_new
  }
  list(beta = beta, iter = max_iter, converged = TRUE,
       fallback = "max_iter", cond_F2 = cond_last)
}

# ====================================================================
# 2. Monte Carlo validation: empirical slope-variance ratio vs g_2(alpha)
# ====================================================================

# Symmetric error generators (var = 1), matching the D6v2 illustration set.
err_gen <- list(
  Laplace        = function(n) sqrt(0.5) * rexp(n) * sample(c(-1, 1), n, TRUE),
  `GG(4)`        = function(n) {
    beta <- 4; c_scale <- sqrt(gamma(1/beta) / gamma(3/beta))
    u <- runif(n); s <- sign(u - 0.5)
    s * c_scale * qgamma(2 * abs(u - 0.5), shape = 1/beta, rate = 1)^(1/beta)
  },
  `Student-t(6)` = function(n) rt(n, df = 6) * sqrt((6 - 2) / 6)
)

# Theoretical g_2(alpha) for the symmetric error law (c_2 = 1; sigma_p = 0).
g2_theory_named <- function(alpha, name) {
  fn <- switch(name,
    Laplace        = list(c2 = 2, nu = function(q) nu_q_laplace(q)),
    `GG(4)`        = list(c2 = 1, nu = function(q) nu_q_gg(q, 4)),
    `Student-t(6)` = list(c2 = 1, nu = function(q) nu_q_t(q, 6)),
    stop("unknown law"))
  g2_alpha(alpha, fn$c2, fn$nu, function(q) 0)
}

run_regression_mc <- function(M = 2000, n = 400, alpha_grid = c(0.05, 0.50, 0.95),
                              beta_true = c(1.0, 2.0), seed = 2026) {
  set.seed(seed)
  out <- list()
  for (law in names(err_gen)) {
    egen <- err_gen[[law]]
    for (alpha in alpha_grid) {
      b_ols  <- matrix(NA_real_, M, length(beta_true))
      b_patp <- matrix(NA_real_, M, length(beta_true))
      for (m in seq_len(M)) {
        x1 <- rnorm(n)
        X  <- cbind(1, x1)
        y  <- as.vector(X %*% beta_true) + egen(n)
        b_ols[m, ]  <- solve(crossprod(X), crossprod(X, y))
        b_patp[m, ] <- patp_regression(X, y, alpha)$beta
      }
      # slope = coefficient 2
      v_ols  <- var(b_ols[, 2])
      v_patp <- var(b_patp[, 2])
      out[[length(out) + 1]] <- data.frame(
        law = law, alpha = alpha, n = n, M = M,
        slope_var_ols  = v_ols,
        slope_var_patp = v_patp,
        g2_emp  = v_patp / v_ols,
        g2_theo = g2_theory_named(alpha, law),
        bias_slope_patp = mean(b_patp[, 2]) - beta_true[2]
      )
    }
    cat(sprintf("regression MC: %s done\n", law))
  }
  bind_rows(out)
}

# ====================================================================
# 3. Run when invoked as a script
# ====================================================================

if (sys.nframe() == 0L) {
  if (!exists(".patp_regression_done")) {
    cat("=== PATP regression Monte Carlo validation (symmetric errors) ===\n")
    reg_mc <- run_regression_mc(M = 2000, n = 400)
    reg_mc_out <- reg_mc %>%
      mutate(across(c(slope_var_ols, slope_var_patp, g2_emp, g2_theo,
                      bias_slope_patp), \(x) round(x, 4)))
    print(reg_mc_out)
    dir.create("results", showWarnings = FALSE)
    write.csv(reg_mc_out, "results/regression_mc.csv", row.names = FALSE)
    cat("\nSaved results/regression_mc.csv\n")
    .patp_regression_done <- TRUE
  }
}
