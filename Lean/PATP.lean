import PATP.Param
import PATP.Basis
import PATP.ExponentSeparation
import PATP.Degeneracy
import PATP.Entropy
import PATP.G2Algebra

/-!
# PATP — Parametrically-Adaptive Transition Polynomial

Entry point for the Lean formalization of the PATP apparatus.

Поточний скоуп:

* `PATP.Param` — функція показників `p_i(α)` та теореми про три
  граничні випадки $\alpha \in \{0, 1/2, 1\}$ (закриває §3.2-§3.3 paper).
* `PATP.Basis` — знакозбережна сім'я $\varphi_i(\xi; \alpha)$
* `PATP.ExponentSeparation` — алгебраїчна факторизація `p_i(α)-p_j(α)`
* `PATP.Degeneracy` — колапс Form-B базису до `ξ` при `α=1/2`
* `PATP.Entropy` — алгебра ентропійного коефіцієнта `k`
* `PATP.G2Algebra` — алгебраїчний крок формули $g_2(\alpha)$ (§4.2)
-/
