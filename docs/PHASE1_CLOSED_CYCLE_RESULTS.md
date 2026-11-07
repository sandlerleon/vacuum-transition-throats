# Phase 1: closed matter-coupled region, no assumed bounce — results

Status: computed (code in `code/closed_cycle.py`, `closed_orbits.py`, `closed_orbit_scan.py`, `closed_stability.py`, `make_closed_figure.py`; data in the matching JSON files; two new automated checks, 17 in total). **Not committed, not released.** Nothing here enters the CQG manuscript.

## What was computed

The Paper 1 scalar (symmetron potential, dark-matter coupling A = 1 + φ²/2M², held vacuum V_f, true vacuum −qV_f) in a **closed** FRW universe (K = +1, curvature radius a), units 8πG = 1, R_dS = 1, with **exact** equations and no bounce prescription:

ä = −(a/6)(ρ + 3P),  φ̈ + 3Hφ̇ + V′_eff = 0,  Friedmann constraint ȧ² = a²ρ/3 − 1 monitored (and imposed exactly on the initial data).

A bounce is whatever minimum of a the equations produce. Start: the classical neck (ȧ = 0, φ = small seed, φ̇ = 0), dust parameter C = ρ_*a³/3 below C_c = 0.385 R_dS.

## Result 1: generic trajectories do not cycle (216 of 216)

216 runs (q = 1; ρ_crit/V_f = 0.02, 0.05, 0.2; μ = 20, 30 H_Λ; C/C_c = 0.1–0.9; seeds 10⁻³ to 10⁻²⁰ v; dissipation Υ = 0, 0.3, 3): every one expands from the neck, holds, releases, **turns around once and then collapses to a singularity**; none has a second bounce. At the turnaround ρ + 3P ≥ 0.16 V_f > 0 (the strong energy condition holds), so the contraction accelerates; the released vacuum energy (scalar kinetic energy up to 1.2 V_f at the turnaround) is not returned to the potential. Dissipation changes nothing here (and would feed radiation, which also defeats the neck).

## Result 2: exact cold cycles exist, as tuned, time-reversal-symmetric orbits

The equations are time-reversal invariant. A trajectory that starts at the neck with φ̇ = 0 and reaches a turnaround with **φ̇ = 0 as well** is mapped to its own time reverse: the contraction retraces the expansion, restoration returns the scalar to the seed, and the universe returns to the neck. That is one condition on one free parameter (the seed), so the orbits are isolated: a **coherent seed**, as in Paper 1. A scan of 18 parameter sets (q = 1; μ = 10, 20; ρ_crit/V_f = 0.05, 0.2, 0.4; C/C_c = 0.3, 0.6, 0.9) found roots in 17 and orbits that return to the start in at least two cases with a large excursion and in five with a small one:

| Kind | Parameters | a_max / a_neck | Return error | Survives (bounces) at relative seed error 0 / 10⁻¹² / 10⁻⁹ / 10⁻⁶ |
|---|---|---|---|---|
| wide cycle (Fig. C1) | μ = 10, ρ_crit = 0.05 V_f, C = 0.9 C_c | 8.6 | 6×10⁻⁸ v | 1 / 1 / 1 / 0 |
| wide cycle | μ = 10, ρ_crit = 0.05 V_f, C = 0.6 C_c | 7.8 | 4×10⁻⁶ v | 1 / 1 / 0 / 0 |
| breather | μ = 20, ρ_crit = 0.2 V_f, C = 0.3 C_c | 1.03 | 3×10⁻¹³ v | 11 / 11 / 8 / 6 |
| breather | μ = 10, ρ_crit = 0.4 V_f, C = 0.6 C_c | 1.25 | 2×10⁻¹² v | 8 / 7 / 5 / 3 |

The wide cycles are the physical ones: the region expands from a = 0.68 R_dS (the neck) to 5.7 R_dS, a held epoch of about two e-folds, a release and restoration around the turnaround (φ makes a brief excursion to about 1.3 v and returns), contraction, and a nonsingular cold bounce at a density of 1–2.2 V_f. The breathers are small-amplitude oscillations of the scalar between the two vacua with a nearly constant scale factor (a_max/a_neck 1.03–1.3), not an epoch of expansion.

## Result 3: every exact cycle found is unstable

Floquet multipliers of the Poincaré map on the neck (area-preserving, det = 1 where well resolved): the five breathers are hyperbolic, with leading multipliers 12.7, 72, 460, 1.2×10³ and 1.6×10⁴ per cycle; none is elliptic. The wide cycles are far more unstable: a relative seed error of 10⁻⁶ (or 10⁻⁹) destroys the second bounce, so the multiplier per cycle is at least of order 10⁶ to 10⁹. In double precision the wide cycles close one full cycle and part of a second (Fig. C1); the longest-surviving breather completes 11 bounces.

## Against the four success criteria

| Criterion | Status |
|---|---|
| Nonsingular bounce (a_min > 0, finite invariants, regular field) | **Met** (exact cycles; neck at 0.68–0.94 R_dS, density 1–2.2 V_f) |
| Repeated turnaround under the actual scalar dynamics | **Met** for the tuned orbits (turnaround recurs) |
| Vacuum restoration without external resets | **Met** for the tuned orbits (φ returns to the seed) |
| Long-term stability | **Not met**: all orbits hyperbolic; wide cycles with multipliers ≳ 10⁶; no attractor, unlike the hot cycle of Paper 1 |

So the answer to the Phase 1 question is: **a classical cold cycle exists without a wormhole, but only as an isolated, tuned, unstable periodic orbit; trajectories near it leave it within one to a few cycles, and generic trajectories collapse after the first turnaround.** That strengthens Paper 1 in one respect (the bounce can be derived, not assumed) and weakens it in another (the cold cycle is not an attractor; Paper 1's stability came from dissipation into radiation, which the neck cannot tolerate).

## Consequences for the wormhole programme

* "A single bounce is not yet a cycle" is borne out: the cycle exists but is not robust.
* Phase 2 now has a sharp, testable form (Result 3 of the proposed outcomes): **does a pre-existing bridge change the stability of these orbits** (the multipliers, or the existence of an elliptic island)? If the multipliers are unchanged, the bridge is an optional interpretation.
* The orbits are regions of size about R_dS = 5.3 Gpc with spatial curvature of order unity; they are not our observed, nearly flat universe (see `PREEXISTING_BRIDGE_CYCLE_REUSE.md`).

## Limits

Homogeneous; scalar plus coupled dust only (no radiation in the orbit search; dissipation tested only on generic trajectories); q = 1; μ = 10–30 H_Λ; the stiff hierarchy of Paper 1 (μ ≫ H) is not reached; a one-parameter tuned family was sampled on a coarse grid (the count of orbits is not a census); instability estimates for the wide cycles are lower bounds set by double precision; entropy and inhomogeneities not treated.

## Next

1. Add radiation/dissipation to the orbit search (does weak dissipation make an attracting cold cycle, or only a hot one?).
2. Refine the stability map over (q, ρ_crit/V_f, μ, C): is there any elliptic island?
3. Phase 2: the bridge, via the two-ended double-null data of `PRIMORDIAL_BRIDGE_PLAN.md`, comparing multipliers with and without it.

---

# Phase 1B: can self-consistent dissipation stabilize the cold cycle? (2026-11-06)

Code: `closed_dissipation.py`, `closed_dissipation_lifetime.py`, extensions of `closed_cycle.py`; data: `closed_dissipation.json`, `closed_dissipation_fine.json`, `closed_dissipation_bestseed.json`, `closed_dissipation_lifetime.json`; check 18 of 18. Still **uncommitted and unreleased**; the CQG manuscript is unchanged.

## Method

Energy is never deleted: φ̈ + (3H + Γ)φ̇ + V′ = 0 and ρ̇_r + 4Hρ_r = Γφ̇², with ρ_r in the Friedmann and acceleration equations. Γ is (1) constant, (2) gated to rapid field motion (|φ̇| above half the peak of the undamped cycle), or (3) thermal, Γ = C_T ρ_r^{3/4} from a radiation seed of 10⁻⁶ V_f. The constraint is satisfied at every start (residuals ≤ 5×10⁻¹³ at the turning points).

## 1. The exact orbits are periodic in the full state (baseline, Γ = 0)

After one period the full state (a, ȧ, φ, φ̇) returns to the start: for the five small-excursion orbits |Δa| ≤ 1.3×10⁻¹³ R_dS, |Δȧ| ≤ 7×10⁻¹⁷, |Δφ|/v ≤ 7×10⁻¹¹ and |Δφ̇|/v ≤ 6×10⁻¹⁰ (limited by the instability amplifying round-off); the wide orbit returns to |Δa| = 7×10⁻⁸ and |Δφ̇|/v = 3×10⁻⁵ for the same reason. So these are periodic orbits, not just T–B–T–B sequences.

## 2. Weak dissipation does not stabilize; it detunes

Stretching of a seed perturbation over one return (the Floquet proxy; Γ = 0 values 67, 132, 1.3×10³, 9.6×10³, 6.6×10⁴) is **unchanged to three digits for Γ ≤ 10⁻³ and rises for Γ ≥ 10⁻² (67 → 72 → 94 at Γ = 0.1, 0.2)**. A tuned orbit loses its later bounces as soon as Γ > 0 (μ = 20, ρ_crit = 0.2 V_f orbit: 11 bounces at Γ = 0, 6 at 10⁻⁶, 3 at 10⁻³): the dissipation is a perturbation of the tuned state and is amplified at the undamped rate, so the lifetime falls only logarithmically with Γ. The wide orbit loses its second bounce at any Γ ≥ 10⁻⁶.

## 3. Dissipation cannot give a periodic orbit: entropy grows at every bounce

The entropy parameter S = ρ_r a⁴ increases at every bounce in every dissipative run (check), e.g. 2.5×10⁻⁴, 5.1×10⁻⁴, 7.7×10⁻⁴ at Γ = 10⁻³; at Γ = 0.3 the radiation at the necks is 14–20 % of V_f. There is no fixed point of the neck-to-neck map once Γ > 0; a Floquet analysis applies only at Γ = 0.

## 4. What dissipation does do (broad seeds, no fine tuning)

300 seeds per Γ over log₁₀(s/v) ∈ [−2.6, −0.3], t_max = 150:

| Family | Γ = 0 | Γ = 0.2 | Γ = 0.5 | Γ = 0.7 | Γ ≥ 1 |
|---|---|---|---|---|---|
| μ = 20, ρ_crit = 0.2 V_f, C = 0.3 C_c: max bounces / fraction of seeds with ≥ 3 | 3 / 0 % | 13 / 2 % | 8 / 6 % | 15 / 8 % | 0–1 / 0 % |
| μ = 10, ρ_crit = 0.4 V_f, C = 0.6 C_c | 8 / 4 % | 3 / 1 % | 1 / 0 % | 1 / 0 % | 1 / 0 % |

In the first family moderate dissipation (Γ of order 0.2–0.7 H_Λ) **raises the share of seeds that make several bounces and the longest lifetime (up to 15 bounces on the grid; a local refinement at Γ = 0.3 found 14, 22, 25 and 35)**, with radiation of order 0.14–0.2 V_f at the necks; in the second family it only reduces them. Gated and thermal dissipation behave like constant dissipation or worse (at most 5 bounces for the thermal case). In **every** run the cycle ends in a singular collapse within about 25 time units (a few tens of bounces): the radiation parameter S keeps growing until the neck fails (e.g. S doubles in the last bounce of the 14-bounce run). No run survives indefinitely, and none approaches a limit cycle with saturating S.

## 5. A structural obstruction

For the undamped orbit with leading multiplier λ, damping can only stabilize it if the phase-space contraction per cycle beats the instability: ΓT ≳ ln λ (for the best breather λ = 12.7, T = 0.69, so Γ ≳ 3.7 H_Λ). But the energy that dissipation takes from the scalar, which is what the restoration needs to climb back to the false vacuum (ΔV = (1 + q)V_f), appears as radiation that the neck cannot tolerate: at Γ ≳ 1 the first return already fails in all 22 runs tried, and at Γ = 3.7 the radiation reaches 0.2 V_f within one period. Within the sampled models the two requirements (damp the instability; keep restoration and a vacuum-dominated neck) are incompatible.

## Verdict on Phase 1B

* A stable, constraint-satisfying, repeated cold cycle was **not found**. Exact periodic orbits exist at Γ = 0 and are hyperbolic; with consistent dissipation there is no periodic orbit, entropy grows each cycle, and the best case is a modest, family-dependent lengthening of the transient (tens of bounces), never a limit cycle.
* The stability obstruction is quantitative: ΓT ≳ ln λ versus radiation below about V_f at the neck.
* Limits: two families (q = 1), homogeneous, three dissipation laws, a coarse seed grid, one-parameter tuned orbits; the thermal and gated laws are illustrative, not derived from a microscopic model; a different energy-transfer channel (for example to a stable dark component that does not gravitate like radiation) is not excluded.

## Consequence for the wormhole programme

The scalar–matter system does not sustain a stable cold cycle on its own, so the milestone for deferring Phase 2 is met. Phase 2 stays a specific hypothesis, not an expectation: whether a pre-existing bridge could change the gravitational boundary conditions or allow a physically consistent energy flux that alters the stability of the closed region. It would need a genuine coupled spacetime solution; adding a damping term and calling it a bridge effect would not test it. The stated outcome, if no stable cycle is found, is a stability obstruction (above) rather than a cyclic-universe mechanism.

---

# What is proven, what is numerical, what is open

**Proven (analytic).** With the equations of Phase 1B, the comoving radiation parameter obeys d(ρ_r a⁴)/dt = Γ a⁴ φ̇² ≥ 0 (from ρ̇_r + 4Hρ_r = Γφ̇²; it follows for any sign of H). If Γφ̇² is not identically zero over a cycle, ρ_r a⁴ increases and the system cannot return to its initial state: **no exact periodic orbit of the full scalar–radiation system exists for Γ ≥ 0, Γ ≢ 0 under these assumptions.** (This is an obstruction to exact recurrence in this dissipative model; it does not exclude other energy-transfer channels or other models.)

**Numerical, within the sampled models (not exhaustive).**
1. The undamped orbits are exact periodic orbits to round-off (full state) and are linearly unstable (multipliers 12.7 to 1.6×10⁴ for the small-excursion family; ≳ 10⁶ for the wide cycles).
2. Moderate dissipation can lengthen the transient (up to 35 bounces, one family) but produced no stable recurrence; every run ended in collapse within about 25 time units.
3. Stabilization by damping needs ΓT ≳ ln λ, which conflicts with a radiation content that the neck can tolerate (radiation ≈ 0.2 V_f within one period at the required Γ).

**Open.** Whether a different geometry or energy-transfer mechanism (a sink that does not gravitate like radiation, an exchange with an exterior region, inhomogeneous dynamics) can stabilize the cycle. Not every classical cold cycle has been shown to be unstable.

# Reproduction (Phase 1 and 1B)

Equations and units: top of `closed_cycle.py` (8πG = 1, R_dS = 1, V_f = 3, K = +1; exact ä = −(a/6)(ρ+3P) with the constraint imposed on the initial data by `neck_a`). Scan: `python closed_cycle.py` (216 generic runs) → `closed_cycle.json`; orbit search `closed_orbits.py`, `closed_orbit_scan.py`; Floquet and return map `closed_stability.py`; dissipation `closed_dissipation.py`, `closed_dissipation_lifetime.py`; figures `make_closed_figure.py`, `make_phase1b_figure.py` (FigC1, FigC2); checks `tests/test_all.py` (18).
