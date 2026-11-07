# Pre-existing bridge → cyclic universe without a big bang: what can be reused

Status: assessment of the existing Zenodo/GitHub material plus one new calculation (`code/closed_bounce.py`). Nothing here is published or released.

## Bottom line

1. **A cyclic universe without a hot big bang already exists in your own work, in a limited regime.** Paper 1 (late-time vacuum decay, §4.1) finds that with the effective bounce at ρ_c ≤ 10⁴ V_c the cycles repeat (growth factor Γ_n ≈ 0.99) with **no hot dense epoch**; at ρ_c = 10⁵ and 10⁶ V_c they grow (Γ_n = 1.19, 2.31). So "cold" cycles are possible if the bounce happens at or below about 10⁴ times the dark-energy density. What fails is not the cycle but the **bounce**: Paper 1 §4.6 says plainly that it is assumed, not derived, and that no established theory predicts a bounce at that density.
2. **The bounce is the only place a pre-existing bridge could help**, and the classical bridge geometry gives a bounce at exactly the right density, but not for our universe (below).
3. **A bridge that carries matter from one universe to another is still excluded** by the null-energy argument of Paper 3 (traversable throat) and by the causal-accessibility result (Section 5.4: no causal curve from the exterior into the core). The workable idea is a bridge in **time**, not in space.

## The new result: the neck of global de Sitter space is a cold, classical, null-energy-respecting bounce

A closed universe with false-vacuum energy V_f and dust of conserved parameter C = ρ_m a³/3 has ȧ² = V_f a²/3 + C/a − 1. On the branch that comes in from large a it **bounces at a_min > 0** (no singularity, no big bang) if and only if

C < C_c = 2R_dS/(3√3) = 0.385 R_dS (the Eddington–Lemaître condition),

with the bounce at a = 1.00 to 0.62 R_dS and density 1.0 to 2.6 V_f as C goes from 0 to C_c. It needs the strong energy condition to fail (the vacuum energy does that) and **not** the null energy condition, so Proposition 1 of Paper 3 is not violated: the minimum is a minimum of the radius along a timelike curve, not of the areal radius along a null ray. The bounce density (about V_f) is far inside Paper 1's cold-cycle regime (≤ 10⁴ V_c).

**But the observed curvature rules it out for our own universe.** Today's matter in a closed universe of curvature Ω_k has C/C_c = Ω_m√Ω_Λ/(0.385 Ω_k^{3/2}):

| Ω_k | 0.001 | 0.004 | 0.01 | 0.1 | 0.5 | 0.8 |
|---|---|---|---|---|---|---|
| C/C_c | 2.1×10⁴ | 2.7×10³ | 6.7×10² | 21 | 1.9 | 0.94 |

A bounce needs C/C_c < 1, i.e. Ω_k ≳ 0.77. The data allow |Ω_k| of order 10⁻³ (Planck 2018 + BAO), so the matter of a nearly flat universe exceeds the critical value by three to four orders of magnitude. The universe must contract by a large factor before reaching the neck, the matter density rises by the cube of that factor, and matter, not vacuum energy, dominates: a singularity, not a bounce. **So the global dS neck cannot be the bounce of our universe.** It can describe an isolated, strongly closed false-vacuum region (the persistent core of Paper 3), which is a different universe from ours.

## What can be reused, piece by piece

| Piece | Where | Use for the new idea | Limits |
|---|---|---|---|
| Held vacuum, release by dilution, turnaround, restoration | Paper 1 §2–3, `matter_coupling.py`, `dilution_background.py` | The whole expansion → recollapse → restoration half of any cycle, independent of the bounce | Derived within the model; needs a coherent seed for the release |
| Cold-bounce cycles, hysteresis law a_max,n = W_n/3, stability criterion PΥ ≳ ln μ | Paper 1 §4.1–4.2, `cyclic_extension.py`, `oscillation.py` | Defines the **cold cycle** regime (bounce at ≲ 10⁴ V_c) that a no-big-bang cycle must live in | Bounce prescribed; cycles above 10⁴ V_c grow |
| Daughter-region integrator with spatial curvature, events (nucleation, restoration, bounce) | Paper 2, `bubble_gravity.py` (class `Daughter`), `bubble_to_cycle.py`, `bubble_handoff.py` | Hand-off of an isolated region into a cycle; add K > 0 and the closed-dS bounce in place of the assumed one | Currently uses the assumed effective bounce |
| Independent constraint-preserving solver; 100-cycle runs without resets; attracting fixed point of the cycle map | Paper 1 §5–6, `independent_solver.py`, `continuous_cycles.py`, `cycle_map.py` | Verification pattern and the cycle-to-cycle map for any new bounce | Homogeneous only |
| Black holes through the bounce, anisotropy tolerance, failure map | Paper 1 §4.3, 4.5, 6.2 | Constraints any bounce (including a bridge bounce) must satisfy: percolation bound, shear dilution | Prescription-dependent |
| Finite-tension junction analysis, bridge branch (Γ_in < 0), critical mass, horizon roots | Paper 3 §4, `thinwall.py` | Level 1 of Table 5: the mathematical existence of the bridge; the global-dS-neck bound C_c is the cosmological analogue of m_c | Thin wall only |
| Double-null solver with regular centre; quadrant diagnostics; constraint and expansion convergence | Paper 3 §5, `dnull.py` | Extend to two-ended data through a bifurcation sphere (stage 2 of `PRIMORDIAL_BRIDGE_PLAN.md`) | Needs new boundary treatment; no dark matter |
| Causal-accessibility argument (u, v non-decreasing) | Paper 3 §5.4 | The test every proposed bridge must pass for Level 3 | Holds for the ingoing-wall solutions only |
| Null-energy proposition | Paper 3 §3 | Rules out traversable spatial throats however they arose; justifies looking for a bounce in time instead | Classical model |
| Voids-first environmental decay, granularity bound | Paper 2 | Where in a cycle the transition starts; not needed for the bridge | Not related to the bounce |

## Model families that could be made to work (ranked)

1. **Region-to-region reproduction without a big bang.** The persistent expanding core of Paper 3 is a cold de Sitter region. In this model a false vacuum needs matter above ρ_crit to be held, so such a core would dilute its matter, release, recollapse and could be the parent of the next one. This reuses the Paper 1 cycle machinery and the Paper 2 integrator. It is reproduction of regions, not a return of our universe; it needs the matter content of the core as an input (and the closed-dS bound C < C_c if the neck is to be the bounce).
2. **Cold cycle with a derived bounce for a closed region.** Replace the assumed bounce of Paper 1 by the closed-dS neck, in the cold regime where Paper 1 already shows stable repeats. Success test: the integrator completes release → recollapse → restoration → neck bounce → held expansion with C < C_c. Failure test (already known for our universe): Ω_k ≳ 0.77 is required.
3. **Black-hole-to-white-hole style bounce.** The causal structure (region II to region IV) is not available classically; it needs the same kind of assumed quantum bounce as Paper 1, with the additional requirement of a local null-energy violation. It adds a geometric picture, not a derivation.

## What not to do

* Do not present the dS neck as the bounce of our flat universe: the numbers above exclude it.
* Do not claim matter transfer through a non-traversable bridge: Table 5, Level 3 is absent in the present solutions.
* Do not treat "undetectable" as support: a hidden bridge needs a prediction that can fail (stability boundary, causal-accessibility criterion, distinction from black-hole collapse).

## Suggested next computation (smallest useful step)

Add a closed-universe option (K > 0, no assumed bounce) to the Paper 2 `Daughter` integrator and scan (C, Ω_k-at-restoration): does any cycle close, and how far is the nearest closing cycle from the observed curvature? This turns the three-to-four-orders-of-magnitude statement above into a map and either yields a viable cold-cycle region or closes the question.
