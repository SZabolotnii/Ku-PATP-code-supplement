# Reproducibility Notes

## Numerical experiments

The R pipeline is deterministic where seeds are used inside the Monte Carlo scripts. Running

```bash
Rscript R/run_all.R
```

from the repository root regenerates the CSV summaries (including the
convergence, regression-coefficient, and real-data outputs) and the eight
manuscript figures. Scripts `06` and `07` carry a top-level run guard and are
launched as separate `Rscript` processes by `run_all.R`.

## Formal checks

The Lean files are a compact formal supplement, not a full mechanization of the statistical theory. They cover structural algebraic facts used in the manuscript:

- the quadratic exponent map `p_i(alpha)`;
- the three special regimes `alpha = 0`, `alpha = 1/2`, and `alpha = 1`;
- signed-parity basis identities;
- degeneracy at `alpha = 1/2`;
- entropy-coefficient algebra;
- a compact algebraic step for the `g_2(alpha)` formula.

Run:

```bash
lake update
lake build
```

to verify the Lean project.
