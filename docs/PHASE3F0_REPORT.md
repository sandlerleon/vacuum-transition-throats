# Phase 3F.0: causal relativistic dissipative fluids and matter-recycling feasibility (branch `phase3-recycling`)

Status: **stopped before expensive evolution, as instructed; a limited feasibility test, not a large program. Verdict: a causal dissipative fluid cannot produce the transition. Bulk viscosity makes a bounce harder; shear viscosity can damp the growth of spacetime anisotropy, but only at the edge of the hydrodynamic regime, and it still leaves the collapse singular; the dissipative stress keeps the null energy condition. Its useful role is to regulate how matter and anisotropy behave once a gravitational bounce mechanism exists.** 38 automated checks pass.

Reproduce (`code/`): `phase3f0_fluid.py` → `phase3f0_fluid.json` (≈30 s), `tests/test_all.py`. κ = 8πG = 1 in numerical parts.

## Part 0. Independent check of the Phase 3E.1 flat-bounce statement (requested)

Compared with the published Einstein–scalar–Gauss–Bonnet cosmology (Guo–Schwarz, Phys. Rev. D 81, 123520: 6H² = φ̇² + 2V + 24ξ̇H³; 2Ḣ = −φ̇² + 4ξ̈H² + 4ξ̇H(2Ḣ − H²); search-verified form, normalization of that action): my symbolic Euler–Lagrange equations equal them **exactly** (E_N = (a³/2)·GS_F; E_a = 3a²·GS_R + (3a²/2)·GS_F; E_φ = −a³·GS_S), the reparametrisation (Noether) identity ȧE_a + φ̇E_φ − N Ė_N = 0 holds with residual 0, and the coefficient of Ḣ in the Raychaudhuri equation at H = 0 is 2 (no degeneracy). Hence Ḣ(H = 0) = −φ̇²/2 for any ξ(φ), V(φ): the extra Gauss–Bonnet terms are all proportional to H. The Phase 3E.1 statement therefore stands, **scoped to canonical scalars, arbitrary ξ and V, flat FLRW**; the 37 earlier tests check consistency, not the theorem, which rests on the derivation plus this literature comparison. The 3E.1 report has been amended accordingly.

## Step 1. Fluid theory

GR coupled to a fluid with T^{μν} = (ρ + p + Π)u^μu^ν + (p + Π)g^{μν} + π^{μν} + 2q^{(μ}u^{ν)}, causal Israel–Stewart-type relaxation: bulk τ_Π Π̇ + Π = −ζΘ (truncated: the further IS coupling terms are omitted), shear τ_π π̇ + π = −2ησ. Causality requires the characteristic speeds c_bulk² = c_s² + ζ/(τ_Π(ρ + p)) ≤ 1 and c_shear² = η/(τ_π(ρ + p)) ≤ 1; hydrodynamic validity requires the Knudsen number K = τ|Θ| ≲ 1. In all runs η and ζ are tied to τ so these hold by construction (c², K are the controls). Entropy production is non-negative by the sign convention ζ, η ≥ 0 (not separately integrated).

## Step 2. Collapse dynamics

**Bulk viscosity, closed universe with Λ (the Lemaître bounce of Phases 1–2).** Dust + Λ = 0.02 (R_dS = 12.25) + k = +1 contracting from a₀ = 30. The inviscid critical dust mass below which a bounce occurs is M_c = 14.142 (analytic Lemaître threshold). With causal bulk viscosity:

| K | c_b² = 0 | c_b² = 0.25 | c_b² = 1 |
|---|---|---|---|
| 0.1 | 14.14 | 12.04 | 7.67 |
| 0.3 | 14.14 | 10.16 | 4.44 |
| 1.0 | 14.14 | 8.41 | **2.69** |

**M_c falls monotonically with the viscosity: bulk viscosity makes the bounce harder (the allowed matter content shrinks by a factor of up to 5).** Derivation: after contraction Π = −ζ∫Θ(t′)e^{−(t−t′)/τ}dt′/τ ≥ 0 (Θ < 0), also at the turning point (the lag keeps Π positive), so Ḣ = −(κ/2)(ρ + p + Π) + k/a² is smaller than in the ideal case: viscosity adds positive effective pressure that gravitates. It cannot create a bounce, and a viscous slowing of the collapse is not a bounce.

**Shear viscosity, axisymmetric Bianchi I collapse (spacetime shear versus fluid shear).** The exact equation in these variables is d(α̇ − β̇)/dt = −Θ(α̇ − β̇) + (3/2)S with S = −(4/3)η(α̇ − β̇) quasi-statically, so σ̇/σ = |Θ| − 2κη. With Θ̇/Θ = (1+w)|Θ|/2 in a power-law contraction, **d ln Σ²/dt = (1 − w)|Θ| − 4κη**, where Σ² = σ²/(Θ²/3). Causality (η = c²τ(ρ+p)), Knudsen (τ|Θ| = K) and Θ² = 3κρ give **net damping of Σ² iff c²K > 3(1 − w)/(4(1 + w))** (quasi-static, K ≲ 1). *A first estimate I made (damping needs κη ≥ |Θ|, threshold (1+w)c²K ≥ 3, "impossible for any causal EOS") was contradicted by the numerical test and corrected; the corrected threshold is milder.*

Numerical test (20 runs; start Σ² = 0.0025, collapse to V = 10⁻³; constraint residual < 10⁻⁵):

| w | inviscid Σ²: start → end | viscous, K ≤ 1, c² ≤ 1 |
|---|---|---|
| 1/3 (radiation) | 0.0025 → **0.200** (×80) | damped only at the edge K = 1, c² = 1 (→ 0.0017); K = 0.3, c² = 1 gives 0.021; K = 1, c² = 0.25 gives 0.075 |
| 1 (stiff) | 0.0025 → 0.0025 (marginal) | damped in all 6 admissible cases (→ 0.0017 … 0.0000) |

The refined inequality agrees with 12 of 12 admissible runs and 18 of 20 overall (K > 1 lies outside the quasi-static derivation). So **viscosity can suppress the anisotropic growth that decides the spin/torsion bounce of Phase 3C.2 (the Σ requirement B N₀² > Σ), but only when the relaxation time and shear-mode speed are both near their maximal values** (c²K above 0.375 for radiation); this is an edge of the hydrodynamic regime, not a safe interior. It does not remove the singularity: isotropic radiation collapse remains singular.

## Step 3–4. Stability and a Goldilocks region

The relevant timescale ordering is τ_relax ≲ τ_collapse (K ≲ 1) with damping rate above the anisotropy growth rate: the window is c²K ∈ (0.375, 1] for radiation, a narrow region at the hydrodynamic boundary; the threshold 3(1 − w)/(4(1 + w)) is below the causal-hydrodynamic maximum c²K = 1 only for w > −1/7, and it widens with stiffness (threshold 0 at w = 1). The coupled matter-plus-gravity perturbation spectrum (Step 3) was not derived; only the homogeneous shear mode and the bulk contribution to the bounce condition were computed. Damping of fluid perturbations is not stability of the gravitational background, and I make no such claim.

## Step 5. Regular transition

No run produced a regular minimum of the areal radius caused by dissipation. **The null energy condition holds in every run**: the minimum over the principal directions of (ρ + p + π_eigenvalue)/(ρ + p) is ≥ 0.979 across all 20 shear runs and Π ≥ 0 in all bulk runs. A causal, hydrodynamic dissipative stress has |Π|, |π|, |q| well below ρ + p (here |Π|/ρ ≤ 0.36 at K = 0.3, c_b² = 1), so T_kk ≥ 0 and the Penrose theorem applies to any spherical collapse of such a fluid that forms a trapped surface: no regular transition without an additional gravitational mechanism supplying NEC violation (which none of Phases 3C–3E.1 delivered in a controlled, stable form). Heat transport in an inhomogeneous core–envelope collapse (energy flux inside a horizon) was not simulated; it can redistribute mass-energy within the trapped region (the loophole noted after 3D.0) but does not change the sign of T_kk in the hydrodynamic regime.

## Step 6. Comparison with 3E.1 and applicability

- The 3E.1 tensor-gradient instability (F_T < 0, c_T² < 0) is a property of the gravitational action; fluid damping does not touch it, as the proposal itself states. No claim is made otherwise.
- Hydrodynamic validity: the useful range K ≲ 1 sits exactly where the strongest viscous effects (K → 1) occur; near the bounce densities (Planckian k_F in 3D.0) the continuum fluid description is not justified. This is a conditional statement, not extrapolated.

## Decision gate

| question | answer |
|---|---|
| does a causal dissipative fluid improve collapse behavior? | **partially**: it can suppress shear growth at the edge of the regime (Σ² 0.0025 → 0.0017 instead of 0.20) |
| does it produce or help a bounce? | **no**: bulk viscosity lowers M_c by up to ×5; NEC preserved (min ratio 0.979) |
| does it cure the 3E.1 instability? | **no** (separate sector) |
| is a regular transition possible with ordinary dissipative matter? | **no** (NEC ⇒ Penrose, trapped surface ⇒ singularity) |
| useful as a regulator for a future gravitational source? | **possibly**: it can relax the anisotropy obstacle that selects against spin/torsion bounces, provided a derived bounce source exists |

## Limits

Truncated Israel–Stewart (no nonlinear IS terms), homogeneous models only (FLRW and axisymmetric Bianchi I), one-parameter closures τ = K/|Θ|, η = c²τ(ρ+p); entropy production not integrated separately; no inhomogeneous core–envelope or turbulent flow; no perturbation spectrum; first-order viscous (Navier–Stokes) theory not evolved because of acausality; the inequality for shear damping is quasi-static and has been checked only for two equations of state.

## Recommended next step

If a derived gravitational bounce mechanism is ever established, couple it to this causal fluid and test whether the shear obstruction identified in 3C.2 is relieved within K ≲ 1; until then fluid dynamics supplies no new source and the programme's open item remains S6 of the 3D.3 census (one concrete limiting-curvature or Horndeski model with c_T² > 0 and an Einstein-frame bounce).
