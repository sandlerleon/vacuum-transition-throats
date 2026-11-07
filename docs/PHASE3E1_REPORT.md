# Phase 3E.1: field-dependent Gauss–Bonnet couplings, nested throats and effective causality (branch `phase3-recycling`)

Status: **stopped after 3E.1 as instructed. Verdict: NO-GO for all three ideas as tested.** (1) A one-loop-type, field-dependent Gauss–Bonnet coupling does not change the flat-collapse no-go (proved for arbitrary coupling functions, potentials and any number of canonical scalars), and the only closed-collapse bounces need the opposite coupling sign, a relative Gauss–Bonnet strength ≥ 4, and are all gradient-unstable (c_T² < 0). (2) A second throat requires the scalar kinetic function to vanish twice, enlarges the local NEC violation, and worsens the test-field instability. (3) Gauss–Bonnet couplings do modify the tensor propagation speed, c_T² = F_T/G_T, from an action, but this does not supply a source. 37 automated checks pass.

Reproduce (`code/`): `phase3e1_corrections.py` → `phase3e1_corrections.json` (≈5 min), `tests/test_all.py`. κ = 1 (M_p² = 1/κ) in the numerical scans; λξ is the Gauss–Bonnet coefficient.

## Part A. One-loop-type correction

**Literature (search-verified secondary sources, not read in full).** In the loop-corrected heterotic action the Gauss–Bonnet term is multiplied by a function of a modulus (Antoniadis–Rizos–Tamvakis, Nucl. Phys. B415 (1994), hep-th/9305025); the tree-level dilaton–Gauss–Bonnet term alone does not violate the energy conditions and cannot give singularity-free solutions; singularity-free solutions are reported for flat backgrounds and a definite sign of the coupling, and Kanti–Rizos–Tamvakis (gr-qc/9806085) find them for ξ = φⁿ with even n. The exact coefficient and form of the modulus function were **not** verified, so I did not adopt a specific one; Part A is proved for **arbitrary** ξ(φ, χ) and V(φ, χ).

**Setup.** Einstein-frame S = ∫√−g[R/2κ − ½Σ(∂φ_i)² − V + λξ(φ_i)G], two canonical scalars (dilaton and modulus). Euler–Lagrange equations for N, a, φ, χ derived symbolically for FLRW(k).

**Exact results (symbolic).**
- Constraint at H = 0: ½(φ̇² + χ̇²) + V = 3k/(κa²). All Gauss–Bonnet terms carry a factor of H (or ȧ), so they drop out of the constraint at the turning point.
- **Flat collapse (k = 0): ä/a = κ(−φ̇² − χ̇² + 2V)/4, which equals −(κ/2)(φ̇² + χ̇²) ≤ 0 on the constraint surface (residual 0).** So **no bounce is possible for any ξ, any V, any number of canonical scalars**: this includes the modulus-dependent loop coupling whatever its form. The singularity-free solutions of the literature are therefore not collapse-to-expansion bounces; consistent with the 3E.0 result, and with nonsingular meaning "no initial singularity" (flat asymptotics), a different property.
- **Closed collapse (k = +1):** on the constraint surface
 ä/a = a²[2k − κa²|φ̇|² − 8κλ k ξ_{ab}φ̇^aφ̇^b] / (2(a⁴ + 96k²κλ²|∇ξ|²)).
 A bounce (ä > 0) needs 2 > κa²|φ̇|² + 8κλ ξ_{ab}φ̇^aφ̇^b. With V ≈ 3/(κa²), small φ̇, it is the classical closed de Sitter bounce of Phases 1–2 (no Gauss–Bonnet role, and the flatness-tuning problem found there). Driven by Gauss–Bonnet instead (V = 0, φ̇² = 6/(κa²)) it requires **λξ_{φφ} < 0 and 48|λξ_{φφ}|/a² > 4, i.e. a relative strength ε_GB ≳ 1/12 of order one**, again where the expansion is uncontrolled.

**Numerical test (single field, ξ = φ², V = 0; 720 contracting initial data, both signs of λ, flat and closed).**

| case | runs | bounces | tensor-sector strength ε_T at bounce | bounces with ε_T < 0.3 | with G_T ≤ 0 or F_T ≤ 0 |
|---|---|---|---|---|---|
| k = 0, λ > 0 | 180 | **0** | – | – | – |
| k = 0, λ < 0 | 180 | **0** | – | – | – |
| k = +1, λ > 0 (heterotic sign with ξ'' > 0) | 180 | **0** | – | – | – |
| k = +1, λ < 0 | 180 | 76 | 4.2 to 4×10⁴ | **0** | **76 of 76** |

Every Gauss–Bonnet-driven closed bounce occurs only for λξ'' < 0 (the sign the literature says is "definite"; whether the actual heterotic one-loop function has this sign is not verified here), at ε_T ≥ 4 (tensor corrections of order one or far larger), and with F_T < 0 (c_T² from −39 to +0.95, median −4.3): **gradient (Laplacian) instability and loss of hyperbolicity in the tensor sector** (G_T > 0, so no ghost). Gradient and scalar-sector stability beyond the tensor sector were not derived.

**Comparison with 3E.0.** The one-loop-type correction does not change the no-go for k = 0 (proved) and gives no controlled closed bounce. Omitted comparable terms: all higher-order α′, loop corrections to other couplings, H-flux, compactification moduli stability, string-frame vs Einstein-frame matter couplings (the Einstein frame was used as physical; a different matter coupling could change the interpretation of proper times).

**Independent check added in Phase 3F.0 (after review).** The flat-bounce statement was checked against the published Guo–Schwarz equations (Phys. Rev. D 81, 123520: 6H² = φ̇² + 2V + 24ξ̇H³; 2Ḣ = −φ̇² + 4ξ̈H² + 4ξ̇H(2Ḣ − H²), for the action R/2 − (∂φ)²/2 − V − ξR²_GB/2): the derived equations equal them exactly (E_N = (a³/2)·GS_F, E_a = 3a²·GS_R + (3a²/2)·GS_F, E_φ = −a³·GS_S), the Noether identity holds with residual 0, and the coefficient of Ḣ at H = 0 is 2 (no degeneracy), so Ḣ(H = 0) = −φ̇²/2 exactly. The theorem's scope remains: canonical scalars, arbitrary ξ and V, flat FLRW. The 37 tests check consistency with these equations; they do not by themselves establish the theorem, which rests on the symbolic derivation above and this literature comparison.

## Part B. Nested throats

- **Topology.** In spherical symmetry two local minima of the areal radius need an intervening maximum. A profile R(r) = √(r² + a²) + β e^{−(r/w)²} (a = 1, β = 1.5, w = 0.8) has minima at r = ±1.117 and a maximum at r = 0 (the belly); R″ changes sign at r = ±0.498.
- **Source requirement (exact identity for any f(r), R(r)):** G^r_r − G^t_t = −2fR″/R. For scalar sources κεfφ′² = −2fR″/R, so the kinetic sign ε = −sign(R″): phantom at each throat, canonical at the belly. **A single scalar needs its kinetic function h(φ) to vanish at each inflection point (two for this nested profile): the perturbation kinetic term vanishes there (infinite effective coupling).** This is not a "stronger inner wormhole"; it is a source with two points of strong coupling. Multi-field sources avoid the zero only by adding more unsupported fields.
- **NEC budget** κ∫(−T_kk)R²dr over the NEC-violating regions: 5.6 (single throat) vs **16.0** (nested), i.e. ×2.8 local violation; the belly region satisfies the NEC with a budget of 11.0, but local NEC violation at the throats cannot be cancelled by it, since the null focusing is local.
- **Coupled-stability test (test-field toy; the only computable one).** The ξR-coupled test scalar (ℓ = 0, m = 0) has a growing bound state above ξ_c. Single throat: ξ_c = 0.511, γ(ξ = 12) = 4.15; nested (β = 1.5): ξ_c = 0.322, γ(12) = 6.67; nested (β = 3, w = 0.6): ξ_c = 0.195, γ(12) = 13.4. **The inner throat lowers the threshold and increases the growth rate: coupling between the throats does not stabilize, it worsens the test-field instability** (the cavity between them binds modes). The coupled perturbation operator of the *source* was not derived because no admissible source exists; the mechanical two-oscillator analogy is not a gravitational result, and I did not substitute it for one.

## Part C. Effective propagation speeds (from an action, not an adjustable c)

For λξ(φ)G in Einstein–scalar–Gauss–Bonnet the tensor modes satisfy G_T ḣ² − F_T(∇h)²/a² with **G_T = 1/κ + 8λHξ̇ and F_T = 1/κ + 8λξ̈**, so **c_T² = F_T/G_T = (1 + 8κλξ̈)/(1 + 8κλHξ̇)**. G_T was derived independently from the exact Bianchi-I shear action (difference 0) and agrees with the Horndeski-mapped value; F_T is the Horndeski-mapped value (checked symbolically, residual 0), not derived from a spatially inhomogeneous action here. The tensor characteristic cone differs from the matter light cone wherever λ(ξ̈ − Hξ̇) ≠ 0: this is a physical change of a dimensionless ratio (c_T/c_matter), not a change of units. Effects found: c_T² > 1 (superluminal gravitational characteristics) and c_T² < 0 (loss of hyperbolicity) occur in the closed Gauss–Bonnet bounces above; the scalar sector's cone was not derived. Modified cones do not by themselves make horizons traversable or supply NEC violation.

## Part D–E. Goldilocks feasibility and decision gate

| requirement | result |
|---|---|
| consistent field equations and constraints | **yes** (symbolic, residual-checked) |
| controlled string coupling and curvature expansion | **no**: the Gauss–Bonnet-driven bounce needs ε_GB ≳ 1/12 and ε_T ≥ 4 |
| regular curvature, finite tidal forces | not reached (no admissible bounce) |
| no ghost or loss of hyperbolicity | **no**: all 76 bounces have F_T < 0 (c_T² < 0) |
| dynamically stable or finite-time-controlled | **no** (gradient instability; nested test-field growth enlarged) |
| continuation from collapse to expansion | only in the unstable/uncontrolled closed case |

1. Does the one-loop correction change the 3E.0 no-go? **No** for flat collapse (theorem for any coupling), **no controlled change** for closed collapse.
2. Are two-throat geometries supported by the theory? **No**: need two zeros of a kinetic function and a ×2.8 larger local NEC violation.
3. Does coupling between throats improve stability? **No** in the computable toy (ξ_c halves or lower, growth rate ×1.6–3.2).
4. Does any field have a physically modified cone? **Yes**, tensor modes (c_T² ≠ 1) from the Gauss–Bonnet action, with c_T² < 0 where the bounce occurs.
5. Is there a controlled candidate? **No.**

## What this does and does not exclude

It excludes controlled bounces in Einstein–scalar–Gauss–Bonnet gravity with canonical scalars for flat collapse (proved) and for closed collapse in the scanned region (opposite-sign, ε_T ≥ 4, c_T² < 0). It does not exclude: the full ten-dimensional heterotic theory, other string-derived terms (H-flux, higher orders, non-canonical moduli kinetic terms, potentials with negative regions combined with extra structure), inhomogeneous collapse, or the scalar-sector analysis; and the exact one-loop modulus function was not implemented.

## Recommended next step

Stop the Gauss–Bonnet family here. If the programme continues, the only structurally different direction left in the census is S6 from 3D.3 (a derivable limiting-curvature or Horndeski model): pick one, derive its homogeneous and spherical collapse, and test hard constraints H1–H7 including c_T² > 0 and the Einstein-frame bounce criterion before any scan.
