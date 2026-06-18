#!/usr/bin/env Rscript
# run_all.R - PATP code-supplement pipeline driver.
#   01 -> theoretical g_2(alpha) for the symmetric canonical laws (CSV)
#   02 -> theoretical figures fig1-fig4
#   05 -> full F^{-1}b estimator definitions (Algorithm 1)
#   03 -> Monte Carlo: full estimator (primary) + proxy comparison + convergence
#   06 -> regression-coefficient validation (Monte Carlo)
#   07 -> real-data application (EuStockMarkets daily log-returns) + bootstrap
#   04 -> experimental figures fig5-fig9
#
# Scripts 01-03 execute their analysis when sourced; 05 only defines functions.
# Scripts 06 and 07 are invoked as separate Rscript processes. 04 is sourced
# LAST because Fig. 9 consumes the bootstrap CSVs written by 07.

scripts_sourced <- c("01_theoretical_g2.R",
                     "02_visualizations.R",
                     "05_full_patp_estimator.R",
                     "03_monte_carlo.R")
scripts_spawned <- c("06_patp_regression.R",
                     "07_real_data_application.R")
scripts_figures <- c("04_results_viz.R")

# Resolve script directory so `Rscript R/run_all.R` from repo root works
# the same as `cd R && Rscript run_all.R`.
script_dir <- {
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep("^--file=", args, value = TRUE)
  if (length(match) > 0) {
    dirname(sub("^--file=", "", match[1]))
  } else if (file.exists("01_theoretical_g2.R")) {
    "."
  } else {
    "R"
  }
}
old_wd <- getwd()
setwd(script_dir)
on.exit(setwd(old_wd), add = TRUE)

t_start <- Sys.time()
for (s in scripts_sourced) {
  cat("\n=== Running", s, "===\n")
  t0 <- Sys.time()
  source(s)
  cat(sprintf("[%.1fs] %s done\n",
              as.numeric(Sys.time() - t0, units = "secs"), s))
}
for (s in scripts_spawned) {
  cat("\n=== Running", s, "(separate process) ===\n")
  t0 <- Sys.time()
  st <- system2("Rscript", s)
  if (st != 0L) stop("failed: ", s)
  cat(sprintf("[%.1fs] %s done\n",
              as.numeric(Sys.time() - t0, units = "secs"), s))
}
for (s in scripts_figures) {
  cat("\n=== Running", s, "===\n")
  t0 <- Sys.time()
  source(s)
  cat(sprintf("[%.1fs] %s done\n",
              as.numeric(Sys.time() - t0, units = "secs"), s))
}
cat(sprintf("\nTotal: %.1fs\n",
            as.numeric(Sys.time() - t_start, units = "secs")))
