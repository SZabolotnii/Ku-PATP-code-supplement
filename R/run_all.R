#!/usr/bin/env Rscript
# run_all.R — Pipeline driver для PATP
# Запускає всі скрипти послідовно:
#   01 → теоретичні g_2(α)
#   02 → теоретичні графіки fig1-fig4
#   05 → визначення повного F^{-1}b PATP естиматора (Phase B)
#   03 → Monte Carlo симуляції (proxy + full estimator)
#   04 → експериментальні графіки fig5-fig8

# Note: 05 is sourced inside 03 as well; listing it separately here
# ensures self-test output appears in the run log.
scripts <- c("01_theoretical_g2.R",
             "02_visualizations.R",
             "05_full_patp_estimator.R",
             "03_monte_carlo.R",
             "04_results_viz.R")

t_start <- Sys.time()
for (s in scripts) {
  cat("\n=== Running", s, "===\n")
  t0 <- Sys.time()
  source(s)
  cat(sprintf("[%.1fs] %s done\n",
              as.numeric(Sys.time() - t0, units = "secs"), s))
}
cat(sprintf("\nTotal: %.1fs\n",
            as.numeric(Sys.time() - t_start, units = "secs")))
