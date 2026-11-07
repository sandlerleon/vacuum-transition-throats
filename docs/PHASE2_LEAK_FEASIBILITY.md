# Phase 2, Stage A: reduced leak-model feasibility test (branch `phase2-bridge`)

Status: **negative for this reduced model.** No wormhole is claimed, and this is not a two-ended solver. It is not part of the CQG manuscript or its validated results.

Reproduce: `python code/bridge_leak.py` (writes `bridge_leak.json`, 122 runs); `python code/tests/test_all.py` (20 checks).

## What was tested

Question: if a pre-existing bridge carries energy out of a closed cold cycle (so dissipative heat is not stored as radiation), do the exact cold cycles of Phase 1 survive, and is a bridge needed?

Interior: the Phase 1 closed model plus dissipation Γ and a radiation leak Q (units 8πG = 1, R_dS = 1).

| law | Q | meaning |
|---|---|---|
| `diffuse` | ε ρ_r | volume emission; any ordinary radiating boundary |
| `area` | ε ρ_r / a | flux through a fixed-size mouth (shell area / volume); the bridge-specific law |
| `sink` | Γ φ̇², ρ_r ≡ 0 | ideal sink, all heat leaves instantly; the best any radiation leak can do |
| `none` | 0 | zero-leak baseline with the same Γ |

Orbits: μ=20, ρ_crit=0.2 V_f, C=0.3 C_c (period 0.692, free-fall: 10 bounces) and μ=10, ρ_crit=0.4 V_f, C=0.6 C_c (free-fall: 6 bounces). Γ ∈ {0, 0.03, 0.1, 0.3, 1, 3, 10}; ε ∈ {0.1, 1, 10, 100}; 40 periods maximum.

Bookkeeping (stated as approximations):
- Energy received by the exterior, E_leak = ∫Q a³ dt, with exterior mass parameter m_leak = (π/2) E_leak.
- A homogeneous interior cannot lose energy and keep the Friedmann constraint ȧ² = a²ρ/3 − 1. So the constraint is *tracked*, not enforced: K_eff = a²ρ/3 − ȧ², with the exact identity dK_eff/dt = −a²Q/3 (checked in the test). The model is self-consistent only while K_eff ≈ 1. A large drift means the homogeneous reduced model has broken down and the inhomogeneous problem (Stage B) would be needed.
- Shell stresses (σ > 0, σ + p ≥ 0) are evaluated quasi-statically with the growing exterior mass at χ0 = 0.7π, q = 1.

## Results

1. **No stabilisation.** The finite-time exponent (growth of a 10⁻⁹ separation of nearby cold-cycle starts until 10⁻³) is positive in every run: 2.2–8.1 for orbit 1 (zero-leak baseline 4.2–7.4) and 1.3–3.1 for orbit 2 (baseline 1.5–2.4). The longest-lived leak run lasts 5 bounces; free fall already gives 10 and 6. No run is periodic.
2. **Self-consistent leaks do not help.** Among runs that keep K_eff ≥ 0.9 (10 runs), none exceeds the zero-leak baseline with the same Γ.
3. **Leaks that "beat" the baseline are not valid.** 19 runs gain a bounce or two over their same-Γ baseline, but 17 of them are in orbit 1 at Γ ≥ 1 with K_eff from −0.5 to −121 (the interior is no longer closed), and the other two (orbit 2: Γ = 0.03 `area`, ε = 100, K_eff down to 0.19; Γ = 0.1 `diffuse`, ε = 10, K_eff down to 0.78) gain one or two bounces while the constraint has already drifted by 20–80%. They are consequences of leaving the constraint surface, not evidence of recurrence.
4. **The ideal sink is no better than a realistic leak.** The `sink` runs (everything dissipated leaves at once) have at most 4 bounces (orbit 2, Γ = 0.03) and at most 2 (orbit 1), with K_eff drifting far from 1. No leak of radiation can do better than the ideal sink in this reduced model, so the dissipation channel itself (Γ drains the scalar's climbing energy) is the limit, not where the heat goes. This is consistent with the Phase 1B obstruction.
5. **Milestone 3 (bridge versus ordinary boundary) cannot be decided at this level.** `diffuse` and `area` leaks give essentially the same interior behaviour. A reduced homogeneous interior cannot distinguish a bridge from any other energy sink. Whatever is special about a bridge has to come from the exterior geometry, junction conditions and causal structure, which this model does not contain.
6. **Shell stresses worsen with the leak.** The fraction of the cycle with σ > 0 falls from 1.0 to 0.16–0.9 at large Γ·ε as the exterior mass grows. The fraction with σ + p ≥ 0 is 0–0.4 there, against 0.5–0.9 for zero leak, so the null energy condition on the shell is violated more, not less.

## Milestones (status)

| milestone | status |
|---|---|
| 1. energy transfer is possible | yes, as bookkeeping, but only by leaving the constraint surface (K_eff drifts) |
| 2. cold cycle becomes long-lived | **no** in this model: lifetime ≤ free fall; exponent never ≤ 0 |
| 3. bridge necessary versus ordinary boundary | **undecided** at this level (reduced model blind to the difference) |

## Recommendation

Do not proceed to a two-ended solver on the strength of this test. Stage B would need to be justified by a new mechanism, such as a source that resupplies the scalar's climbing energy, not by more of the same leak. Phase 2 stays on its separate branch; nothing here enters the manuscript or its validated results.
