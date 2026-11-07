# Phase 3E.0: string-inspired high-curvature collapse (branch `phase3-recycling`)

Status: **stopped before any ten-dimensional or axisymmetric evolution, as instructed. Verdict: NO-GO for the model actually tested (heterotic tree-level string frame, 4D reduction with fixed internal volume, dilaton + Gauss–Bonnet). No physical regular transition was found: flat and closed collapse cannot bounce at all (two exact results), and the open-slicing bounces that do occur are leading-order string-frame artifacts that survive when the α′ term is switched off, vanish in the Einstein frame, and run to strong coupling.** 36 automated checks pass.

Reproduce (`code/`): `phase3e0_string.py` → `phase3e0_string.json` (≈6 min), `tests/test_all.py`. Units α′ = 1, so ε = α′|curvature| = |curvature|.

## Part A. The model and its reduction (what is and is not fundamental)

String-frame action S = (1/2κ₄²)∫√−g e^{−2Φ}[R + 4(∂Φ)² + cG], c = α′/8 (heterotic sign), G the Gauss–Bonnet scalar, FLRW with curvature k. **This is a deliberately schematic truncation**, as in the proposal: the NS–NS three-form H, gauge fields, the internal-volume modulus, one-loop and all higher-order α′ corrections are dropped, and the six-dimensional compactification is replaced by a torus with frozen volume. It is not the heterotic ten-dimensional theory, and it is not claimed to be equivalent to it; the proper reduction and field-redefinition analysis remain undone.

Verified: √−g G = 8 d/dt[ȧ³/N³ + 3kȧ/N] exactly; the closed forms of R and G agree with a brute-force curvature computation of the S³-sliced metric (residual 0 and 0). Euler–Lagrange equations for N, a, Φ derived symbolically; the constraint is
2a e^{−2Φ}[ 3a²H² − 6a²HΦ̇ + 2a²Φ̇² − 24cHΦ̇(a²H² + k) + 3k ] = 0.
Check: with c = 0, k = 0 it is the shifted-dilaton constraint φ̄̇² = 3H², and the pre-big-bang power law a ∝ (−t)^{−1/√3} solves all three equations.

## Part B. Exact results

1. **Flat collapse (k = 0): no bounce.** At H = 0 the constraint is 4a³Φ̇² = 0, so Φ̇ = 0, and the a- and Φ-equations then give ä = aΦ̈ and Φ̈ = 3ä/(2a), hence **ä = 0 = Φ̈**: H = 0 is a static point, never a turning point.
2. **Closed collapse (k = +1): no bounce.** All Gauss–Bonnet terms in the constraint carry a factor of H. At H = 0 it reduces to 4a³Φ̇² + 6a = 0, i.e. **Φ̇² = −3/(2a²) < 0**: impossible. This is the geometry class of a collapsing ball of ordinary matter (Oppenheimer–Snyder interior has k ≥ 0).
3. **Open slicing (k = −1)** permits H = 0 with Φ̇² = 3/(2a²), and ä/a = 2(a² − 12c)/(a⁴ − 24a²c + 48c²) > 0 for a² > 12c = 1.5. The leading-order part (c → 0) of this is 2/a² > 0 already.
4. **Scale argument (derived from the constraint).** The ratio of the Gauss–Bonnet term to the Einstein–Hilbert terms is ≈ 8cΦ̇H(1 + k/(a²H²)) = α′HΦ̇(1 + …): the α′ term competes with the leading terms only where **ε ~ O(1)**, i.e. exactly where the truncated expansion is uncontrolled. No scale in this model separates "strong enough to bounce" from "controlled".
5. **Shear-sector ghost condition (exact quadratic action in Bianchi I).** The coefficient of b̈ in the linearised shear equation is a³e^{−2Φ}(32cHΦ̇ − 4): healthy iff **1 − 8cHΦ̇ > 0, i.e. α′HΦ̇ < 1**. The ghost boundary coincides with ε = O(1). (The gradient coefficient and the scalar sector were not derived here.)

## Part C–D. Homogeneous feasibility scan (480 initial data, both Φ̇ roots, contracting, c = α′/8)

| k | runs with real data | string-frame bounces | also without the GB term | Einstein-frame bounces | shear-ghost runs |
|---|---|---|---|---|---|
| 0 | 160 | **0** | – | 0 | 3 |
| +1 | 32 | **0** | – | 0 | 4 |
| −1 | 160 | 143 | **143 of 143** | **0** | 67 |

- **k = ±0, +1:** every run is singular or never turns, as the exact results predict.
- **k = −1:** 143 of 160 runs show an H = 0 turning point in the string frame, with ε at the turning point from 5×10⁻⁴ to 0.52 (114 below 0.1; 45 below 0.1 along the whole path). But **all 143 occur identically with c = 0**, so the α′ (Gauss–Bonnet) correction plays no role; **none is a bounce of the Einstein-frame scale factor a_E = e^{−Φ}a** (0 of 160): the string-frame bounce at H = 0 has Φ̇ ≠ 0, so the Einstein-frame universe is still contracting. The physical gravitational frame is the Einstein frame (4D Planck units fixed); a string-frame turning point with a contracting Einstein-frame scale factor is a frame artifact, and the dilaton runs to strong coupling (g_s up to 9×10¹² along these trajectories), where loop corrections that were dropped dominate. Also, an open k = −1 slicing is not the interior of a collapsing star.

## Part E. Comparison with earlier phases (does it resolve or shift the obstruction?)

| earlier obstruction | status in this model |
|---|---|
| 3C Weyssenhoff spin-fluid c_s² / μ = 0 pathology | not applicable (no fluid); replaced by a different absence of bounce |
| 3D.0 ECD Planck-scale breakdown | replaced by the string scale; the Gauss–Bonnet term matters only at ε ~ 1: same control problem |
| 3D.3 missing physical source | the tree-level model provides a derivable action but no NEC-violating bounce for k ≥ 0; ghost boundary at ε ~ 1 |
| 3D.1/3D.2 prescribed geometry | not reproduced: no regular transition derived |
| formation of an expanding daughter | unresolved |

**The obstruction is shifted, not resolved**: the new scale is α′, but the correction competes with gravity only at ε = O(1) where the truncation fails, and the structural results (Part B 1–2) are independent of that scale.

## Decision gate

| requirement | result |
|---|---|
| action and reduced equations internally consistent | **yes** (identities and limits verified) |
| solution inside a justified validity regime | **no** for any bounce candidate: GB acts only at ε ~ 1; the k = −1 string-frame turning points at ε ≪ 1 do not depend on GB and run to strong coupling |
| regular transition candidate | **none** (flat/closed exact no-go; open case = Einstein-frame no bounce) |
| no immediate fatal instability | not reached; shear ghost at α′HΦ̇ ≥ 1 |
| plausibly embeddable in inhomogeneous collapse | not applicable |
| **Decision** | **NO-GO for this truncated model; do not proceed toward numerical collapse** |

## What this does and does not exclude

It does not exclude string theory, or even string-inspired bounces: dilaton potentials, one-loop-corrected Gauss–Bonnet couplings (the Antoniadis–Rizos–Tamvakis / Kanti–Rizos–Tamvakis line of work), H-flux, moduli, non-perturbative effects and the actual ten-dimensional compactification were not included; the known "graceful exit" difficulty of the lowest-order string cosmology is recalled from the literature, not re-verified here. What it does show is that the lowest-order tree-level heterotic structure with a Gauss–Bonnet term cannot by itself produce a controlled, physical, collapse-to-expansion transition for the geometry relevant to a collapsing star, and that apparent bounces must be checked in the Einstein frame.

## Recommended next step (bounded)

If string-inspired gravity is to continue: add exactly one more ingredient with a derivable origin, the one-loop dilaton-dependent Gauss–Bonnet coupling (with a justified sign and strength), and repeat this scan with the Einstein-frame bounce as the pass criterion, adding the coupling's own consistency bounds (g_s ≲ 1, ε ≲ 0.1). Only if a bounce with k ≥ 0 appears under those constraints is an inhomogeneous collapse calculation justified.

## Limits

Truncated 4D tree-level model, homogeneous only; both roots of the constraint but a random (not exhaustive) sample of 480 initial data; gradient and scalar-sector stability not derived; ghost condition derived only for the shear sector; no loop or higher-order α′ terms; no check against the ten-dimensional equations.
