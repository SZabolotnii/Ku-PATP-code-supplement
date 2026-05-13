#!/usr/bin/env Rscript
# run_all.R -- pipeline driver for the PATP numerical supplement.
# Runs all scripts in sequence:
#   01 -> theoretical g_2(alpha)
#   02 -> theoretical figures fig1-fig4
#   03 -> Monte Carlo simulations
#   04 -> experimental figures fig5-fig7

args_file <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_file[grepl(file_arg, args_file)][1])
if (!is.na(script_path)) {
  setwd(dirname(normalizePath(script_path)))
}

scripts <- c("01_theoretical_g2.R",
             "02_visualizations.R",
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
