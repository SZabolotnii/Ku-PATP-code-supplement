import Mathlib.Tactic
import PATP.Param

/-!
# PATP.ExponentSeparation — algebraic separation of PATP exponents

This module records the algebraic layer behind the manuscript claim that the
quadratic PATP bridge has a common exponent-collision point at `α = 1/2`.
It deliberately does not formalize measure-theoretic linear independence or
positive definiteness of moment matrices.
-/

namespace PATP

/-- Auxiliary quadratic factor controlling `p_i(α) - p_j(α)`.

This is a pure polynomial expression in `α` with rational coefficients in `i, j`; it does not
involve `1/i`, so it stays computable even though `PATP.piAlpha` is `noncomputable`. -/
def exponentCollisionFactor (i j : ℕ) (α : ℝ) : ℝ :=
  1 + (((i : ℝ) * (j : ℝ)) - 3) * α
    - (2 * ((i : ℝ) * (j : ℝ)) - 2) * α ^ 2

/-- Factorization of the exponent difference.

For nonzero indices,
`p_i(α) - p_j(α)` is the product of `(j-i)/(ij)` and a single quadratic
factor. This is the algebraic statement used in the paper before reducing
the root analysis to the interval `[0,1]`.
-/
theorem piAlpha_sub_eq_factor (i j : ℕ) (hi : i ≠ 0) (hj : j ≠ 0) (α : ℝ) :
    piAlpha i α - piAlpha j α =
      (((j : ℝ) - (i : ℝ)) / ((i : ℝ) * (j : ℝ))) *
        exponentCollisionFactor i j α := by
  unfold piAlpha exponentCollisionFactor
  have hi' : (i : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hi
  have hj' : (j : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hj
  field_simp [hi', hj']
  ring

/-- The common midpoint collision: all nonzero exponents equal at `α = 1/2`. -/
theorem piAlpha_sub_at_half_eq_zero (i j : ℕ) (hi : i ≠ 0) (hj : j ≠ 0) :
    piAlpha i (1 / 2) - piAlpha j (1 / 2) = 0 := by
  rw [piAlpha_at_half i hi, piAlpha_at_half j hj]
  ring

/-- The collision factor vanishes at `α = 1/2` for any indices. -/
theorem exponentCollisionFactor_at_half (i j : ℕ) :
    exponentCollisionFactor i j (1 / 2) = 0 := by
  unfold exponentCollisionFactor
  ring

/-- Explicit second algebraic root of the collision factor when `ij ≠ 1`. -/
theorem exponentCollisionFactor_at_second_root
    (i j : ℕ) (hij : (i : ℝ) * (j : ℝ) ≠ 1) :
    exponentCollisionFactor i j (-(1 / (((i : ℝ) * (j : ℝ)) - 1))) = 0 := by
  unfold exponentCollisionFactor
  -- `field_simp` needs the denominator `ij − 1` in `≠ 0` form, not `ij ≠ 1`.
  have hij' : (i : ℝ) * (j : ℝ) - 1 ≠ 0 := sub_ne_zero.mpr hij
  field_simp [hij']
  ring

end PATP
