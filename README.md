# PATP Code Supplement

This repository contains the reproducibility code for the paper:

**Parametrically Adaptive Transition Polynomial: a Signed-Parity Continuous-alpha Extension of Kunchenko Stochastic Polynomials**

The repository is intentionally separate from the manuscript repository. It contains only code, computed CSV summaries, rendered figures, and Lean formalization files that are suitable for public release.

## Contents

- `R/` - numerical pipeline for the theoretical curves, the Monte Carlo experiments, the large-N convergence validation, robust baselines, alpha ablations, runtime diagnostics, the regression-coefficient check, and the real-data application. The illustrations are restricted to symmetric error laws (Laplace, GG(1.5), GG(4), Student-t(6)), the regime in which the closed-form g_2(alpha) is exact.
  - `R/05_full_patp_estimator.R` - reference implementation of Algorithm 1 (full F_2^{-1} b normal-equation solver via one-step Newey linearisation); this is the headline estimator used throughout.
  - `R/03_monte_carlo.R` - Monte Carlo driver: the full estimator is primary; the scalar M-estimator proxy is reported alongside for comparison (fig8); a large-N convergence sweep (M=5000, N up to 4000) writes `convergence_g2.csv`.
  - `R/06_patp_regression.R` - confirms by Monte Carlo that the same g_2(alpha) factor governs the variance of regression coefficients (`regression_mc.csv`).
  - `R/07_real_data_application.R` - real-data location example on `EuStockMarkets` daily log-returns; selects alpha* by a bootstrap criterion and reports the realised variance reduction (`real_data_application.csv`, `real_data_g2_curve.csv`).
- `R/results/` - CSV outputs used by the manuscript figures and tables; includes `monte_carlo*.csv` (with `*_full` columns), `convergence_g2.csv`, `regression_mc.csv`, and `real_data_*.csv`.
- `figures/` - rendered PDF figures `fig1`-`fig8`; `fig6_validation.pdf` shows the full estimator's empirical g_2 converging to the closed form as N grows, and `fig8_proxy_vs_full.pdf` contrasts the full estimator with the scalar proxy against the theoretical curve.
- `Lean/` - Lean 4 formalization of the PATP exponent map, signed-parity basis facts, degeneracy at `alpha = 1/2`, entropy-coefficient algebra, and the compact algebraic step behind the `g_2(alpha)` formula.
- `lakefile.lean`, `lean-toolchain` - Lean project configuration.

## Reproduce the R outputs

Install the R packages used by the scripts:

```r
install.packages(c("ggplot2", "dplyr", "tidyr"))
```

Run the full numerical pipeline from the repository root:

```bash
Rscript R/run_all.R
```

The scripts regenerate:

- CSV summaries under `R/results/` (theoretical, Monte Carlo, convergence,
  robust baselines, alpha ablation, runtime, regression, and real-data)
- PDF figures `fig1`-`fig8` under `figures/`

Wall-clock on a single thread is a few minutes end-to-end; the dominant costs
are the convergence sweep in `03_monte_carlo.R` (5000 replicates up to N=4000)
and the bootstrap in `07_real_data_application.R`. Scripts `06` and `07` carry a
top-level run guard, so `run_all.R` launches them as separate `Rscript`
processes.

## Check the Lean formalization

Install Lean via `elan`, then run:

```bash
lake update
lake build
```

The project is pinned to `leanprover/lean4:v4.26.0` and mathlib `v4.26.0`.

## Citation

If this code is useful, cite the associated PATP manuscript and use the metadata in `CITATION.cff`.

## License

Code and reproducibility scripts are released under the MIT License. Rendered figures and CSV summaries are included to support reproducibility of the manuscript.
