import Mathlib.Tactic
import PATP.Basis

/-!
# PATP.Degeneracy — algebraic collapse at α = 1/2

This module captures the basis-structural part of the degeneracy claim:
at the midpoint `α = 1/2`, every Form-B PATP basis function equals the
linear residual `ξ`. Matrix rank and probabilistic positive-definiteness
claims are left to the manuscript-level assumptions.
-/

namespace PATP

/-- Any two nonzero PATP basis functions coincide pointwise at `α = 1/2`. -/
theorem patpBasis_pair_eq_at_half
    (i j : ℕ) (hi : i ≠ 0) (hj : j ≠ 0) (ξ : ℝ) :
    patpBasis i (1 / 2) ξ = patpBasis j (1 / 2) ξ := by
  rw [patpBasis_at_half i hi ξ, patpBasis_at_half j hj ξ]

/-- Any nonzero PATP basis function is collinear with the linear basis at `α = 1/2`. -/
theorem patpBasis_collinear_with_linear_at_half
    (i : ℕ) (hi : i ≠ 0) (ξ c : ℝ) :
    c * patpBasis i (1 / 2) ξ = c * ξ := by
  rw [patpBasis_at_half i hi ξ]

/-- The difference of any two nonzero PATP basis functions vanishes at `α = 1/2`. -/
theorem patpBasis_pair_sub_eq_zero_at_half
    (i j : ℕ) (hi : i ≠ 0) (hj : j ≠ 0) (ξ : ℝ) :
    patpBasis i (1 / 2) ξ - patpBasis j (1 / 2) ξ = 0 := by
  rw [patpBasis_pair_eq_at_half i j hi hj ξ]
  ring

end PATP
