# PATP Code Supplement

This repository contains the reproducibility code for the paper:

**Parametrically Adaptive Transition Polynomial: a Signed-Parity Continuous-alpha Extension of Kunchenko Stochastic Polynomials**

The repository is intentionally separate from the manuscript repository. It contains only code, computed CSV summaries, rendered figures, and Lean formalization files that are suitable for public release.

## Contents

- `R/` - numerical pipeline for the theoretical curves, Monte Carlo experiments, robust baselines, alpha ablations, and runtime diagnostics.
  - `R/05_full_patp_estimator.R` - reference implementation of Algorithm 1 (full F_2^{-1} b normal-equation solver via one-step Newey linearisation, with automatic scalar-proxy fallback when F_2 is ill-conditioned).
  - `R/03_monte_carlo.R` - Monte Carlo driver that runs both the proxy and the full estimator in parallel and emits side-by-side columns.
- `R/results/` - CSV outputs used by the manuscript figures and tables; includes `monte_carlo*.csv` with `*_full` columns for the full estimator.
- `figures/` - rendered PDF figures `fig1`-`fig8`, including the new `fig8_proxy_vs_full.pdf` comparing the empirical g_2 of the proxy and the full estimator against the closed-form theoretical curve.
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

- CSV summaries under `R/results/`
- PDF figures `fig1`-`fig8` under `figures/`

Wall-clock on a single thread is approximately 40 seconds end-to-end; the
dominant cost is the 1000-replicate Monte Carlo sweep in
`03_monte_carlo.R`, which now executes both the proxy and the full
F_2^{-1} b estimator for each grid point.

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
