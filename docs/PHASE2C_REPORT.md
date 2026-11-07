# Phase 2C: can a moving, dynamically coupled bridge support a stable nonsingular cycle? (branch `phase2-bridge`)

Status: **no cyclic wormhole is demonstrated, and the bridge is not needed for a stable cold cycle.** The results below are numerical (reproducible, 21 automated checks) except where marked *derived*. Nothing here touches the CQG manuscript.

Reproduce: `bridge_moving_shell.py`, `bridge_reservoir.py <m_chi>`, `bridge_reservoir_followup.py`, `bridge_phase2c.py`, `single_elliptic.py`, `single_elliptic_survival.py`, then `tests/test_all.py` (all in `code/`).

## 1. Summary

| question | answer |
|---|---|
| Can a moving (non-comoving) thin shell be made consistent with the closed cold cycle? | Yes. Equations derived from the Israel conditions; flux law verified to 1e-6; periodic shells exist, some with elliptic return maps. |
| Can a gravitational-only bridge return energy to the scalar field? | **No, structurally** (§3): an exactly homogeneous interior is blind to any shell trajectory. |
| Can a reversible, action-derived exchange with a bridge-associated reservoir give a stable cycle? | A stable (elliptic) orbit with reservoir amplitude exists (m_χ = 20, g = 600–2000). |
| Is the bridge/reservoir what makes it stable? | **No** (§6). The same stability exists with the reservoir decoupled and in the plain single-field model at C/Cc ≈ 0.35. |
| Correction to Phase 1 | Elliptic exact cold cycles exist (§7); Phase 1's "all sampled orbits are hyperbolic" holds only for its sampled parameter sets. |

## 2. Moving thin shell (2C-3, junction consistency)

Interior: exactly homogeneous FRW, ds² = −dt² + a²(dχ² + sin²χ dΩ²), units 8πG = 1, R_dS = 1. Shell worldline (t(τ), χ(τ)) with rapidity β relative to the fluid: dt/dτ = cosh β, a dχ/dτ = sinh β.

- R = a sin χ, R′ = ȧ sin χ cosh β + cos χ sinh β
- Γ_in = cos χ cosh β + ȧ sin χ sinh β (Γ_in² = 1 + R′² − ρR²/3 exactly, by the Friedmann constraint)
- K_in^τ_τ = β′ + H sinh β
- exterior (mass m, vacuum density −qV_f): Γ_out = −√F, F = 1 + R′² − m/R + qR², K_out^τ_τ = −(2R″ + m/R² + 2qR)/(2√F)
- Israel conditions: σ = 2(Γ_in − Γ_out)/R, p = K_out^τ_τ − K_in^τ_τ + (Γ_out − Γ_in)/R

With a shell equation of state p = wσ this is a closed second-order system for β′ (R″ is linear in β′). No energy source or sink is added. **Independent check:** the Gauss–Codazzi flux law (σR²)′ + p(R²)′ = −R²(ρ + P) cosh β sinh β holds to ≤ 1.0×10⁻⁶ along all 86 scan runs that completed a period (4 of the 90 w = −1 runs ended early); the opposite sign disagrees at O(1). The equation of motion never used this law, so it verifies the derivation (including the sign of K_in).

Results (Phase 1 orbit μ = 20, ρ_crit = 0.2 V_f, C = 0.3 C_c; `bridge_moving_shell.json`):
- A moving shell with an ordinary equation of state (w ∈ [−1/2, 1]) keeps σ > 0 and σ + p ≥ 0 for the whole cycle in the non-periodic scan runs; w = −1 violates the null energy condition for 17–49% of the cycle in every completed run. The earlier comoving shell violated it in 57 of 72 cases, so motion relative to the fluid removes most of that obstruction.
- Time-symmetric **periodic** shells (β = 0 at the neck and the turnaround) exist for 4 of 30 settings. Three have elliptic return maps (|μ| = 1, det 1); one is hyperbolic (2.38, 0.42).
- Around the stable reservoir orbit of §4: 2 periodic shells in the same 30 settings (one elliptic, NEC satisfied for 79% of the cycle; one hyperbolic, 53%).

Geometry (item 7): in both reservoir-orbit shells the shell radius is 4.3–7.4 times the exterior horizon radius r_h (r_h = 0.171, 0.049). The exterior side is a static untrapped region (ℓ± R = R′ ± Γ_out has opposite signs). The only minimal sphere is the exterior bifurcation sphere behind a horizon. **This is an ordinary Schwarzschild–(A)dS exterior; no traversable throat, no causal access to a second region.** The null-energy proposition of the manuscript excludes one without added stress-energy.

## 3. Structural obstruction (derived)

The interior equations for (a, φ, ρ_r) contain no shell variable. For an exactly homogeneous interior, the scalar return map and its Floquet multipliers are therefore identical for every shell trajectory, comoving or moving, and for every shell equation of state. A gravitational junction exchanges energy with the shell and exterior only (through volume change and the flux law above); it cannot act on φ. To act on φ a bridge needs a **non-gravitational** coupling: a bulk interaction (§4) or a boundary action that sources inhomogeneous modes (which needs the PDE problem, i.e. the full two-ended solver).

## 4. Reversible exchange from an action (2C-1, 2C-2)

Reservoir χ with the most general renormalisable Z₂ action: L = −½(∂φ)² − V₀ − Aρ_* − ½(∂χ)² − ½m²χ² − ½gφ²χ². Then J_φ = −gφχ², χ̈ + 3Hχ̇ + (m² + gφ²)χ = 0, total energy–momentum is conserved and the Friedmann constraint holds identically (checked: residual ≤ 1e-12). Nothing is tuned to force periodicity: periodic orbits are found by 2-D shooting on (φ₀, χ₀) (φ̇ = χ̇ = 0 at the turnaround), and the 4-D return map (φ, φ̇, χ, χ̇ at the next bounce) gives the multipliers.

Scan (5 reservoir masses × 5 couplings, `bridge_reservoir_*.json`): orbits with χ₀ > 0 mostly have large multipliers (up to 10⁵). One window is linearly stable: **m_χ = 20, g = 600…2000** has all four multipliers on the unit circle (|μ| − 1 ≈ 3×10⁻⁹; phases 2.712 and 0.0448 rad at g = 1000; same for four Jacobian step sizes). Orbit at g = 1000: φ₀/v = 0.23175, χ₀ = 0.04523, T = 0.62113, a ∈ [0.8867, 0.8979], reservoir energy fraction 10.70–10.86%, constraint residual 7×10⁻¹³. Nonlinear check: perturbations up to 10⁻³ survive 400 bounces with deviation proportional to the perturbation; 1% is the edge. The window in m_χ at g = 1000 is [≈ 20.0, ≈ 20.5]; below 19.75 the multiplier rises (1.5 at 19.75, 3.6 at 19.0), and the orbit is lost at 20.75.

The reservoir exchanges only ≈ 0.2% of the total energy over a cycle: it is stored, not shuttled.

## 5. What the exchange does *not* do

The reservoir energy is nearly constant over the cycle (10.70–10.86%): χ is a fast, small-amplitude oscillator (m_χ ≫ 1/T) that behaves like additional dust-like energy. The question is whether this, rather than any bridge property, is what stabilises the orbit.

## 6. Controls (2C-4): the bridge is not necessary

1. **Coupling off.** The stable orbit's initial data evolved with g = 0 gives 399 bounces (400 with g = 1000). The coupling is not what keeps it going.
2. **Decoupled orbit.** At g = 0 with the same χ₀ the φ-sector has an orbit at φ₀/v = 0.2239 with multipliers 1.0006, 1.0001 (near-elliptic; the 4-D state does not close exactly because χ is not commensurate).
3. **No reservoir at all.** The plain single-field model (Phase 1 code `closed_cycle.py`, `closed_stability.py`, `closed_orbits.py`; no χ, no shell, no bridge) has an **elliptic exact cold cycle at C/Cc = 0.350**: φ₀/v = 0.22295, T = 0.59636, return error 1.6×10⁻¹⁴, trace −1.7043, det 1.000000, |μ| = 1, and likewise at C/Cc = 0.355 (trace −1.408) and 0.360 (−1.164). Nonlinear survival (neck map iterated 200 bounces from perturbed starts): all three seeds survive for perturbations up to 3% (neck deviation ≈ 3.8× the perturbation, no growth); at 10% one of three escapes.

So the success criterion of 2C-4 (stable regime disappears when the bridge coupling is removed) is **not met**. The stable regime is a property of the ordinary closed matter-coupled scalar model, with or without a reservoir, and is simply tuned in by the amount of dust-like energy.

## 7. Correction to Phase 1

The orbit branch that Phase 1 found at C/Cc = 0.30 (φ₀/v = 0.071, multipliers 12.7, 0.079) continues smoothly with C: trace −12.7 (0.300), −6.9 (0.315), −3.7 (0.330), −2.51 (0.340), −2.07 (0.345), then −1.70 (0.350). The trace crosses −2 near C/Cc ≈ 0.346: a period-doubling (flip) bifurcation in which the hyperbolic orbit becomes elliptic. Phase 1's orbit searches sampled C/Cc = 0.3, 0.6, 0.9 (the trajectory scan also 0.1) and so never saw the elliptic window. Consequences:
- "The exact cold cycles are hyperbolic" is true only for the sampled parameter sets. Linearly stable (and, numerically, nonlinearly stable over 200 bounces) exact cold cycles exist in the conservative closed model.
- The Phase 1B result is **unchanged**: with any irreversible dissipation the entropy parameter ρ_r a⁴ is non-decreasing, so no exact periodic orbit exists for Γ ≢ 0. What Phase 1B does not say is that the conservative model has no stable recurrence.
- `docs/PHASE1_CLOSED_CYCLE_RESULTS.md` and the v1.3.0 Zenodo description state the hyperbolic-only conclusion; they need a corrected sentence before v1.3.0 is published. This is left for Leon's decision (the instruction was not to modify Phase 1).

## 8. Milestones

| milestone | status |
|---|---|
| 1. energy transfer possible and reversible | yes for the reservoir (conservative, constraint satisfied) but tiny (≈ 0.2% per cycle); impossible for a purely gravitational shell (§3) |
| 2. cold cycle long-lived | **yes in the conservative closed model, with or without a bridge** (elliptic orbit, 200–400 bounces under perturbation) |
| 3. bridge physically necessary | **no** (§6) |
| junction consistency | moving shell consistent; periodic shells exist; NEC holds for 53–79% of the cycle in the reservoir case, 100% for non-periodic ordinary-EOS runs |
| genuine spacetime bridge | **not demonstrated**: ordinary static exterior behind a horizon, no traversable throat |

## 9. What this implies

1. The stability problem of Phase 1 was a parameter-window issue, not an intrinsic obstruction of the conservative model; the remaining obstruction to a cyclic universe is **irreversibility** (Phase 1B): any dissipation accumulates entropy in the cycle, and a reversible mechanism that removes entropy (rather than energy the scalar needs) would be required. No such mechanism has been identified; the reservoir here conserves entropy and does not remove any.
2. A wormhole is not needed for the stable conservative cycle, and the homogeneous reduction cannot say anything wormhole-specific. If a wormhole is to matter, the physics must be in the inhomogeneous, two-ended problem: a boundary action coupling φ to the bridge so that the exterior can absorb entropy.
3. Requirements for a full two-ended solver (not justified by this report): (i) a boundary action for φ at the shell, (ii) double-null evolution on both sides with the shell as a moving timelike boundary and the flux law above, (iii) a null-expansion map on both sides to test traversability, (iv) a conserved-entropy ledger. The reduced results do not show that this is promising.
