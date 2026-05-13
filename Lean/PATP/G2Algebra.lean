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

/-- If `b^T F^{-1} b = A / Delta`, then
`1 / (c2 * b^T F^{-1} b) = Delta / (c2 * A)`.

This is the final algebraic step used in the manuscript to obtain the boxed
formula for `g_2(alpha)` from the quadratic form. -/
theorem g2_from_quadratic_form
    (c2 A Delta : ℝ) (hc2 : c2 ≠ 0) (hA : A ≠ 0) (hDelta : Delta ≠ 0) :
    1 / (c2 * (A / Delta)) = Delta / (c2 * A) := by
  field_simp [hc2, hA, hDelta]

/-- Expanded version of the previous theorem with the PATP moment symbols. -/
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
  exact g2_from_quadratic_form
    c2
    (bFbNumerator p c2 nupm1 nupp1 nu2p sigp)
    (delta2 c2 nu2p sigp nupp1)
    hc2 hA hDelta

end PATP
