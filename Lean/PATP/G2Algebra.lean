import Mathlib.Tactic

/-!
# PATP.G2Algebra -- algebraic core of the S = 2 variance-reduction formula

This module verifies only the symbolic algebra behind the displayed
`g_2(alpha)` formula.  The moment symbols are abstract real variables;
measure-theoretic existence of the corresponding expectations is part of the
manuscript assumptions, not of this Lean file.
-/

namespace PATP

/-- Determinant of the PATP S = 2 centered correlant matrix. -/
noncomputable def delta2
    (c2 nu2p sigp nupp1 : ℝ) : ℝ :=
  c2 * (nu2p - sigp ^ 2) - nupp1 ^ 2

/-- Numerator of `b^T F^{-1} b` for the S = 2 PATP matrix. -/
noncomputable def bFbNumerator
    (p c2 nupm1 nupp1 nu2p sigp : ℝ) : ℝ :=
  (nu2p - sigp ^ 2)
    - 2 * p * nupp1 * nupm1
    + p ^ 2 * c2 * nupm1 ^ 2

/-- Algebraic closed form of the S = 2 PATP variance-reduction coefficient. -/
noncomputable def g2Formula
    (p c2 nupm1 nupp1 nu2p sigp : ℝ) : ℝ :=
  delta2 c2 nu2p sigp nupp1
    / (c2 * bFbNumerator p c2 nupm1 nupp1 nu2p sigp)

/-- Pure algebraic identity: if `x = A / Δ`, then `1 / (c · x) = Δ / (c · A)`.

This is the internal helper for `g2Formula_eq_from_bFb`. Despite its position in this file,
it is not specific to `g_2(α)` — it is the dividend-inversion identity used to obtain the
boxed manuscript formula from the quadratic-form representation. -/
theorem inv_c_mul_div_eq
    (c A Δ : ℝ) (hc : c ≠ 0) (hA : A ≠ 0) (hΔ : Δ ≠ 0) :
    1 / (c * (A / Δ)) = Δ / (c * A) := by
  field_simp [hc, hA, hΔ]

/-- **Manuscript step §4.2 (closed form of `g_2(α)`).**

Given the quadratic-form representation `b^T F^{-1} b = bFbNumerator / delta2`, the
variance-reduction coefficient `g_2(α) := 1 / (c2 · b^T F^{-1} b)` equals
`g2Formula = delta2 / (c2 · bFbNumerator)`.

The three nonvanishing hypotheses correspond to the following manuscript conditions in §4.2:
* `hc2 : c2 ≠ 0` — second cumulant of the noise is nonzero (excludes degenerate noise);
* `hA : bFbNumerator ≠ 0` — the PATP quadratic form is nondegenerate at the chosen `α`,
  which fails exactly at the midpoint `α = 1/2` (see `PATP.Degeneracy`);
* `hΔ : delta2 ≠ 0` — the centered correlant matrix is invertible (positive-definite under
  the standard PMM moment assumptions of §2.4). -/
theorem g2Formula_eq_from_bFb
    (p c2 nupm1 nupp1 nu2p sigp : ℝ)
    (hc2 : c2 ≠ 0)
    (hA : bFbNumerator p c2 nupm1 nupp1 nu2p sigp ≠ 0)
    (hDelta : delta2 c2 nu2p sigp nupp1 ≠ 0) :
    1 / (c2 *
        (bFbNumerator p c2 nupm1 nupp1 nu2p sigp
          / delta2 c2 nu2p sigp nupp1))
      = g2Formula p c2 nupm1 nupp1 nu2p sigp := by
  unfold g2Formula
  exact inv_c_mul_div_eq
    c2
    (bFbNumerator p c2 nupm1 nupp1 nu2p sigp)
    (delta2 c2 nu2p sigp nupp1)
    hc2 hA hDelta

end PATP
