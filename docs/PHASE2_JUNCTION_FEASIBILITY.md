# Phase 2, step 1: junction feasibility of a pre-existing bridge around the exact cold cycles

Branch `phase2-bridge`. Code `code/bridge_junction.py`; data `code/bridge_junction.json`; check 19 of 19. Exploratory; not part of the CQG manuscript.

## Question

Phase 1/1B found exact classical cold cycles in a closed region (unstable) and no stable cycle with dissipation. Phase 2 asks whether a **pre-existing bridge changes the stability** of that closed region, and whether a self-consistent bridge configuration exists at all with the same interior parameters (identical interior for the bridge and no-bridge cases). Before any two-ended solver, a tractable junction model decides what is worth evolving.

## Model

The closed region (FRW, K = +1, curvature radius a(t), the exact orbits of Phase 1) is cut at polar angle χ₀ **beyond its equator** (π/2 < χ₀ < π, the baby-universe branch of Paper 3, Section 4) and glued across a spherical thin shell to a static Schwarzschild–(A)dS exterior of constant mass parameter m and vacuum density −qV_f (q = 1 AdS, 0 Minkowski, −0.5 lower de Sitter). The shell is comoving with the interior, so the interior stays exactly homogeneous. Units 8πG = 1, R_dS = 1.

Israel conditions, all in closed form: Γ_in = cos χ₀ (exact for the FRW interior), Γ_out = −√F with F = 1 + Ṙ² − m/R + qR², R = a sin χ₀,

σ = 2(Γ_in − Γ_out)/R,  p = −(σR²)˙/(2RṘ) = −(cos χ₀ + √F)/R − (2R̈ + m/R² + 2qR)/(2√F),

with shell energy conservation (σR²)˙ + p(R²)˙ = 0 (no flux across a comoving wall; checked numerically). Requirements over the whole cycle: σ > 0, σ + p ≥ 0 (null energy), constant m > 0. σ > 0 for all t is equivalent to m < m_max = min_t R(sin²χ₀ + Ṙ² + qR²).

## Result 1 (analytic): the comoving bridge is dynamically inert

Because the wall is comoving and the interior is homogeneous, the interior obeys exactly the equations of Phase 1: ä = −(a/6)(ρ + 3P), the scalar equation, the Friedmann constraint. The shell absorbs whatever the exterior requires. Hence **the return map of the interior, and its Floquet multipliers, are identical with and without the bridge**: λ_max^bridge = λ_max^no-bridge for every χ₀, m, q. This is the "different name for the de Sitter neck" outcome in its sharpest form: a bridge that does not exchange energy with the interior cannot change its stability. A bridge can matter only through an energy or momentum flux across the wall (a non-comoving wall, an inhomogeneous interior, or an exchange with the exterior).

## Result 2 (numerical): the required shell is nearly a domain wall and usually violates the null energy condition

For 6 exact orbits (five small-excursion, one wide), 3 exterior depths and 4 cap sizes (χ₀/π = 0.55, 0.7, 0.85, 0.95), 72 combinations, scanning m in 60 steps up to m_max:

* σ > 0 for all m < m_max (by construction), with σ of order 0.03–3 R_dS⁻¹; for m → 0 the shell is positive-tension and p ≈ −σ, a domain-wall-like stress.
* **σ + p ≥ 0 over the whole cycle is satisfied for some m in only 15 of 72 combinations: 1 of 24 with an AdS exterior (q = 1), 6 of 24 with a Minkowski exterior, 8 of 24 with a de Sitter exterior of energy 0.5V_f**, almost all at the smallest cap (χ₀ = 0.55π, barely beyond the equator) and small m.
* For the best m in each combination the worst NEC margin min(σ + p)/max σ has median −0.19 (AdS), −0.17 (Minkowski), −0.09 (de Sitter): typically the shell needs a stress that violates the null energy condition by 9–19 % of its tension at some phase of the cycle.
* The wide cycle with the de Sitter exterior admits no positive-tension shell (m_max < 0).

This agrees with the no-go of the manuscript: a bridge held homogeneous around the cycle needs a wall that cannot be built from the model's canonical scalar and non-negative dust. Where it can (small caps, lower-energy exteriors) the wall barely exceeds the equator, which is the regime least like a bridge.

## What this decides

* **Model A** (pre-existing non-traversable bridge, no energy exchange): inert; no change of stability is possible. It is not worth a two-ended solver for the stability question.
* **Model B** (bridge coupled to an energy reservoir): the only route. The coupling must be a consistent flux with junction conditions and global conservation (for example outgoing radiation through a non-comoving wall into an exterior of increasing mass, Vaidya-type). It forces the interior out of exact homogeneity, so it needs the Einstein–scalar system with dust and radiation in a spherically symmetric setting. The existing double-null solver (regular centre, exterior) is the natural base but lacks dust and null-fluid sectors.
* **Model C** (traversable or quantum-corrected): excluded by the null-energy result unless a stress-energy correction is specified.

## Next step (proposed)

Specify the flux model and test it on the simplest consistent version: an interior of closed FRW plus radiation that can escape through a boundary of area 4πR² into an exterior Vaidya–(A)dS region, with the constraint residual tracked as a diagnostic. Decisive output: the change in the return-map stability (λ_max with the leak against without) for the exact orbits. If radiation leaks without destroying the neck and the multipliers fall below 1, the bridge (as a sink) has a necessary role; if not, that is a second stability obstruction.
