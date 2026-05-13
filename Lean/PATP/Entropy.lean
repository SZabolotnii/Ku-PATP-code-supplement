import Mathlib.Tactic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecialFunctions.Log.Basic

/-!
# PATP.Entropy — ентропійний коефіцієнт як другий параметр форми

Формалізує метрологічну конструкцію Новицького-Зографа~(1991) щодо
ентропійного коефіцієнта $k$ симетричних розподілів, що використовується
у PATP-калібруванні як третій параметр-disambiguator поряд із
$(\gamma_3, \gamma_4)$ (див.~§2.4 paper).

## Означення

* `entropicError` $H \mapsto \Delta_{\text{е}} = e^H / 2$ (Новицький, §3.4)
* `entropyCoefficient` $(H, \sigma) \mapsto k = \Delta_{\text{е}} / \sigma = e^H / (2\sigma)$

## Доведені властивості

* `entropyCoefficient_via_error`: $k = \Delta_{\text{е}} / \sigma$
* `entropyCoefficient_scale_invariant`: $k(c\xi) = k(\xi)$ для $c > 0$
* `entropyCoefficient_pos`: $k > 0$ для $\sigma > 0$

## Канонічні розподіли (значення $k$)

* `entropyCoefficient_gaussian`:  $k_{\mathrm{Gauss}} = \sqrt{2\pi e}/2 \approx 2.0663$
* `entropyCoefficient_uniform`:   $k_{\mathrm{Unif}}  = \sqrt{3} \approx 1.7321$
* `entropyCoefficient_laplace`:   $k_{\mathrm{Lap}}   = e/\sqrt{2} \approx 1.9230$

Конкретні значення ентропій $H$ для цих розподілів узяті з табличних
формул інформаційної теорії~(\citealp{coverthomas2006}, гл.~12) і
тут постульовані як аксіоматичні \texttt{def}-визначення; виведення
$H$ через інтеграл $-\int f \ln f$ потребує \texttt{MeasureTheory} і
винесено в майбутній модуль \texttt{PATP.EntropyIntegral}.
-/

namespace PATP

open Real

/-! ## 1. Базові означення -/

/-- Ентропійне значення похибки за Новицьким:
$\Delta_{\text{е}} = e^H / 2$. -/
noncomputable def entropicError (H : ℝ) : ℝ :=
  Real.exp H / 2

/-- Ентропійний коефіцієнт: $k = e^H / (2\sigma)$.

Безрозмірна величина, що характеризує форму розподілу (інваріант до
масштабу і зсуву) і виокремлюється поряд з ексцесом $\gamma_4$ як
другий незалежний параметр форми (Новицький-Зограф 1991). -/
noncomputable def entropyCoefficient (H σ : ℝ) : ℝ :=
  Real.exp H / (2 * σ)

/-! ## 2. Алгебраїчні властивості -/

/-- Альтернативна форма через ентропійну похибку: $k = \Delta_{\text{е}}/\sigma$. -/
theorem entropyCoefficient_via_error (H σ : ℝ) (hσ : σ ≠ 0) :
    entropyCoefficient H σ = entropicError H / σ := by
  unfold entropyCoefficient entropicError
  field_simp

/-- Додатність ентропійного коефіцієнта при $\sigma > 0$. -/
theorem entropyCoefficient_pos (H σ : ℝ) (hσ : 0 < σ) :
    0 < entropyCoefficient H σ := by
  unfold entropyCoefficient
  positivity

/-- **Інваріантність до масштабу.** Якщо $\xi' = c\xi$ з $c > 0$, то
$H(\xi') = H(\xi) + \ln c$ і $\sigma(\xi') = c\sigma(\xi)$. Тоді
$k(\xi') = k(\xi)$.

Це фундаментальна властивість, що дозволяє розглядати $k$ як параметр
форми, а не масштабу. -/
theorem entropyCoefficient_scale_invariant
    (H σ c : ℝ) (hσ : 0 < σ) (hc : 0 < c) :
    entropyCoefficient (H + Real.log c) (c * σ) = entropyCoefficient H σ := by
  unfold entropyCoefficient
  rw [Real.exp_add, Real.exp_log hc]
  field_simp

/-! ## 3. Ентропії канонічних розподілів

Ці визначення подають \emph{табличні} формули $H$ для канонічних
розподілів без проходження через інтеграл $-\int f \ln f$. Виведення
з інтегрального визначення винесено у майбутній модуль
\texttt{PATP.EntropyIntegral}.
-/

/-- Диференціальна ентропія гаусівського розподілу $\mathcal{N}(0, \sigma^2)$:
$$H_{\mathrm{Gauss}}(\sigma) = \tfrac{1}{2}\ln(2\pi e \sigma^2).$$
-/
noncomputable def gaussianEntropy (σ : ℝ) : ℝ :=
  (1/2) * Real.log (2 * Real.pi * Real.exp 1 * σ^2)

/-- Диференціальна ентропія рівномірного розподілу $U[-a/2, a/2]$
зі стандартним відхиленням $\sigma = a/\sqrt{12}$:
$$H_{\mathrm{Unif}}(a) = \ln a.$$
-/
noncomputable def uniformEntropy (a : ℝ) : ℝ :=
  Real.log a

/-- Диференціальна ентропія розподілу Лапласа з параметром $b$
($f(x) = (1/2b)\,e^{-|x|/b}$), для якого $\sigma = b\sqrt{2}$:
$$H_{\mathrm{Lap}}(b) = \ln(2 b e).$$
-/
noncomputable def laplaceEntropy (b : ℝ) : ℝ :=
  Real.log (2 * b * Real.exp 1)

/-! ## 4. Допоміжна лема: $e^{(1/2) \ln x} = \sqrt{x}$ для $x > 0$ -/

/-- $\exp(\tfrac{1}{2}\ln x) = \sqrt{x}$ для додатних $x$. -/
theorem exp_half_log (x : ℝ) (hx : 0 < x) :
    Real.exp ((1/2 : ℝ) * Real.log x) = Real.sqrt x := by
  have h1 : (1/2 : ℝ) * Real.log x = Real.log x / 2 := by ring
  rw [h1]
  rw [show Real.log x / 2 = Real.log (Real.sqrt x) by
        rw [Real.log_sqrt hx.le]]
  exact Real.exp_log (Real.sqrt_pos.mpr hx)

/-! ## 5. Канонічні значення $k$ -/

/-- **Ентропійний коефіцієнт гаусівського розподілу:**
$$k_{\mathrm{Gauss}} = \frac{\sqrt{2\pi e}}{2} \approx 2.0663.$$
Максимальне можливе значення $k$ серед усіх розподілів з фіксованою
дисперсією (наслідок принципу максимуму ентропії). -/
theorem entropyCoefficient_gaussian (σ : ℝ) (hσ : 0 < σ) :
    entropyCoefficient (gaussianEntropy σ) σ
      = Real.sqrt (2 * Real.pi * Real.exp 1) / 2 := by
  unfold entropyCoefficient gaussianEntropy
  have h_pos : (0 : ℝ) < 2 * Real.pi * Real.exp 1 * σ^2 := by positivity
  rw [exp_half_log _ h_pos]
  -- Goal: sqrt(2 π e σ²) / (2σ) = sqrt(2 π e) / 2
  rw [show (2 * Real.pi * Real.exp 1 * σ^2)
        = (2 * Real.pi * Real.exp 1) * σ^2 from by ring]
  rw [Real.sqrt_mul (by positivity)]
  rw [Real.sqrt_sq hσ.le]
  field_simp

/-- **Ентропійний коефіцієнт рівномірного розподілу $U[-a/2, a/2]$:**
$$k_{\mathrm{Unif}} = \sqrt{3} \approx 1.7321.$$
Обчислюється з $H = \ln a$ і $\sigma = a/\sqrt{12}$. -/
theorem entropyCoefficient_uniform (a : ℝ) (ha : 0 < a) :
    entropyCoefficient (uniformEntropy a) (a / Real.sqrt 12)
      = Real.sqrt 3 := by
  unfold entropyCoefficient uniformEntropy
  have h_sqrt12_pos : 0 < Real.sqrt 12 := Real.sqrt_pos.mpr (by norm_num)
  rw [Real.exp_log ha]
  -- Goal: a / (2 * (a / sqrt 12)) = sqrt 3
  rw [show 2 * (a / Real.sqrt 12) = 2 * a / Real.sqrt 12 by ring]
  rw [div_div_eq_mul_div]
  rw [show a * Real.sqrt 12 / (2 * a) = Real.sqrt 12 / 2 from by
    field_simp]
  -- sqrt 12 / 2 = sqrt 3 ?
  rw [show (12 : ℝ) = 4 * 3 from by norm_num]
  rw [Real.sqrt_mul (by norm_num)]
  rw [show Real.sqrt 4 = 2 from by
    rw [show (4 : ℝ) = 2^2 from by norm_num]
    exact Real.sqrt_sq (by norm_num)]
  field_simp

/-- **Ентропійний коефіцієнт розподілу Лапласа з параметром $b$:**
$$k_{\mathrm{Lap}} = \frac{e}{\sqrt{2}} \approx 1.9230.$$
Обчислюється з $H = \ln(2be)$ і $\sigma = b\sqrt{2}$. -/
theorem entropyCoefficient_laplace (b : ℝ) (hb : 0 < b) :
    entropyCoefficient (laplaceEntropy b) (b * Real.sqrt 2)
      = Real.exp 1 / Real.sqrt 2 := by
  unfold entropyCoefficient laplaceEntropy
  have h_2be_pos : (0 : ℝ) < 2 * b * Real.exp 1 := by positivity
  have h_sqrt2_pos : 0 < Real.sqrt 2 := Real.sqrt_pos.mpr (by norm_num)
  rw [Real.exp_log h_2be_pos]
  -- Goal: (2 * b * exp 1) / (2 * (b * sqrt 2)) = exp 1 / sqrt 2
  field_simp

/-! ## 6. Конкретні приклади -/

/-- **Приклад:** $k_{\mathrm{Gauss}}(\sigma = 1) = \sqrt{2\pi e}/2$. -/
example : entropyCoefficient (gaussianEntropy 1) 1
        = Real.sqrt (2 * Real.pi * Real.exp 1) / 2 :=
  entropyCoefficient_gaussian 1 one_pos

/-- **Приклад:** $k_{\mathrm{Unif}}(a = 1) = \sqrt{3}$. -/
example : entropyCoefficient (uniformEntropy 1) (1 / Real.sqrt 12) = Real.sqrt 3 :=
  entropyCoefficient_uniform 1 one_pos

/-- **Приклад:** $k_{\mathrm{Lap}}(b = 1) = e/\sqrt{2}$. -/
example : entropyCoefficient (laplaceEntropy 1) (1 * Real.sqrt 2)
        = Real.exp 1 / Real.sqrt 2 :=
  entropyCoefficient_laplace 1 one_pos

end PATP
