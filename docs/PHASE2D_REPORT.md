# Phase 2D: a stable cyclic region inside a two-ended bridge geometry (branch `phase2-bridge`)

Status: **Level 2 reached numerically for a non-traversable bridge; Level 3 (matter transfer) is excluded in this geometry; no traversable wormhole is claimed.** The cycle is the plain classical single-field cold cycle of Phase 1 (C/Cc ≈ 0.35). No dissipation, leakage or reservoir is used. Nothing here modifies the CQG manuscript or the Phase 1 documents. 22 automated checks pass.

Reproduce (all in `code/`): `phase2d_cycle_validation.py`, `phase2d_winding.py`, `phase2d_bifurcation.py`, `phase2d_window_end.py`, `phase2d_shell_search.py 0` and `1`, `phase2d_bridge.py`, `phase2d_convergence.py`, `make_phase2d_figure.py` (writes `figures/FigD1_phase2d.png`), then `tests/test_all.py`.

## 1. Independent reproduction of the stable cycle (item 1)

`phase2d_cycle_validation.py` does not import any Phase 1 module. Equations are re-coded from the model definition, the clock is cosmic time for the flow and conformal time for a second orbit solve, the integrators are DOP853, LSODA and Radau, and stability is computed three ways.

| quantity (C/Cc = 0.35, μ = 20, ρ_crit = 0.2 V_f, q = 1) | value |
|---|---|
| orbit φ₀/v (φ̇ = 0 at neck and turnaround) | 0.222948413319 (Phase 1: 0.22294841) |
| period T | 0.596356281246 (cosmic time); same from the conformal-time clock |
| full-state return error (a, ȧ, φ, φ̇) | 1.6×10⁻¹⁴ |
| Friedmann constraint at start / end of the period | 1.1×10⁻¹⁶ / 4.4×10⁻¹⁶ |
| neck-map trace, det | −1.704348, 1.0000001 |
| 4-D monodromy eigenvalues (variational equations) | −0.852174 ± 0.523259i, 1, 1 (|μ| = 1 to 1.5×10⁻⁷) |
| rotation number θ/2π | 0.4124 (close to the 2/5 resonance) |

Convergence: φ₀ and T agree to 10⁻¹² between rtol = 10⁻¹⁰ and 10⁻¹³; return error 1.2×10⁻⁹ (rtol 10⁻⁸) → 1.3×10⁻¹⁴ (10⁻¹³). LSODA, Radau and DOP853 at rtol 10⁻¹¹ agree on T to 10 digits. Same agreement at C/Cc = 0.355 (trace −1.40817, rotation 0.3743) and 0.360 (trace −1.16426, rotation 0.349).

The monodromy has two eigenvalues equal to 1 (the flow direction and the constraint-violating direction) and the elliptic pair; the elliptic pair equals the neck-map multipliers, an independent confirmation of the stability classification.

**What "stable" means here.** Linear (elliptic) stability of a time-reversal-symmetric conservative system means bounded perturbations, not attraction. Numerical evidence for nonlinear boundedness (Phase 2C/2D): the neck map iterated 200 bounces from perturbed starts stays bounded for perturbations up to 3% (deviation ≈ 3.8 × the perturbation), one of three starts escapes at 10%. The neck-map portrait (`FigD1`, panel d) shows nested invariant curves around the fixed point and island chains at larger amplitude. This is evidence, not a proof; the twist condition needed for a KAM argument has not been checked, and the rotation number lies near a resonance.

## 2. The bifurcation near C/Cc = 0.346 (item 2): the label "period-doubling" is NOT supported

- The main orbit branch is smooth in C: φ₀/v = 0.1421 (0.330), 0.1993 (0.345), 0.2229 (0.350), 0.2651 (0.3575), 0.2818 (0.360).
- Trace of the neck map: −3.74 (0.330), −2.51 (0.340), −2.07 (0.345), −1.70 (0.350). It equals −2 at **C/Cc* = 0.345844** (root-found). There a double multiplier −1 occurs: hyperbolic-with-reflection on the lower side, elliptic with rotation number just below 1/2 on the upper side.
- A period-doubling would produce a local 2-cycle emerging from the fixed point. **None was found.** The winding number of d(z) = P²(z) − z on circles of radius 0.002–0.05 v around the fixed point equals the index of the fixed point itself (−1 below C*, +1 above) at C/Cc = 0.3418, 0.3448, 0.3468, 0.3498, 0.3558 (`phase2d_winding.json`), so P² has no other zero within |δz| < 0.02–0.05 v. A multi-start search finds symmetric 2-cycles only at large amplitude (|δφ| ≈ 0.38–0.47 v, φ̇ = 0) that do not shrink to the fixed point as C → C*.
- Conclusion: the multiplier passes through −1 and the orbit changes type, but the standard flip scenario is not confirmed. The Phase 2C wording "period-doubling" is withdrawn and corrected in `PHASE2C_REPORT.md`. The correct statement is "the trace crosses −2".
- **End of the elliptic window.** The branch continues to C/Cc = 0.3600 (φ₀/v = 0.2818, trace −1.176, rotation 0.349, still elliptic). Beyond φ₀/v ≈ 0.284 the start state no longer admits a Friedmann neck (`neck_a` undefined), so the family ends at the boundary of the set of cold bounces, not at a resonance. The elliptic window is therefore C/Cc ∈ (0.3458, ≈ 0.360) for this (μ, ρ_crit, q); Phase 1's sampled values (0.3, 0.6, 0.9) lie outside it. The point at C/Cc = 0.3625 found by the extension scan is a different orbit (T = 1.717, trace 10⁶).

## 3. The geometry (items 3–4)

Specified, not inferred. Parent: maximally extended Schwarzschild with f(R) = 1 − m/R, Kruskal regions I (right exterior), II (black hole), III (left exterior), IV (white hole); in III with future-directed Killing time s: U = e^{κ(r* + s)}, V = −e^{κ(r* − s)}. A spherical timelike shell S = {R = R_s(τ)} is placed in III. **Keep** the part of III with R < R_s together with I, II, IV; **discard** the left asymptotic end (R > R_s). Glue the closed FRW cap (polar angle χ₀ > π/2 of the S³ of the cold cycle) to the kept exterior across S. Junction orientation: Γ_in = cos χ cosh β + ȧ sin χ sinh β (< 0 here, ≈ −1.0), Γ_out = −√F (< 0): the areal radius decreases across the shell toward the throat on both sides. Israel conditions and the shell equation of motion are those of `bridge_moving_shell.py` (Phase 2C, with the flux law verified); equation of state p = wσ.

On a spacelike slice through the cap, the shell and the bifurcation sphere, the areal radius rises from 0 to a maximum a (the equator of the cap, a marginally trapped sphere at the necks), decreases through the shell, reaches its minimum r_h = m (the throat, area 4π r_h²) at the bifurcation sphere and then increases without bound into region I. There are two ends (the closed cyclic region and the asymptotically flat region I). Complete initial data on this slice: cap interior (Friedmann constraint satisfied to 10⁻¹⁶), the shell (σ, p from the Israel conditions) and vacuum exterior.

## 4. Results (items 5–8)

Periodic moving shells around the C/Cc = 0.35 cycle (shell return map in (χ, β) over one interior period, `phase2d_shells_q0.json`). Scan: q ∈ {0, 1}, m/m_max ∈ {0.1, 0.3, 0.5, 0.7, 0.85}, w ∈ {−1/2, −1/4, 0, 1/3, 1}. **Only w = −1/2 with a flat exterior (q = 0) gave periodic shells (5 of 50 settings); no AdS exterior (q = 1) shell was found.** All five have return-map multipliers |μ| = 1 (elliptic), σ > 0 and σ + p ≥ 0 over the whole cycle.

| m/m_max | χ₀/π | m = r_h | shell R | R/r_h min | σ range | max abs β | long-run survival |
|---|---|---|---|---|---|---|---|
| 0.10 | 0.8934 | 0.00334 | 0.309–0.384 | 93 | 0.32–0.42 | 0.42 | 60 periods, all perturbations ≤ 1% |
| 0.30 | 0.9321 | 0.00268 | 0.199–0.314 | 74 | 0.16–0.33 | 0.68 | not run |
| 0.50 | 0.9440 | 0.00252 | 0.165–0.302 | 65 | 0.09–0.31 | 0.85 | not run |
| 0.70 | 0.9505 | 0.00245 | 0.146–0.298 | 60 | 0.05–0.31 | 0.98 | not run |
| 0.85 | 0.9538 | 0.00242 | 0.136–0.297 | 56 | 0.02–0.31 | 1.09 | 60 periods to 10⁻⁴; **fails at 10⁻³** (σ → 10⁻⁵, |β| ≈ 8 within 6 periods) |

So the stiffest shell is linearly but not practically stable. The conservative shell (m/m_max = 0.10) survives 60 interior periods (≈ 36 time units) under random perturbations of the interior neck state and of the shell data (χ, β) of 10⁻⁶, 10⁻⁴, 10⁻³ and 10⁻², keeping R/r_h ≥ 89, σ ≥ 0.30 and σ + p ≥ 0.

Throat, horizons, expansions (conservative shell): throat radius r_h = 0.00334 (the cap's equator is a ≈ 0.95), exterior horizon surface gravity κ = 150. The shell is untrapped on both sides at all times (null expansions R′ ± Γ have opposite signs). The only apparent horizons are the cap's equator (marginally trapped at the necks and turnarounds, where ȧ = 0) and the Schwarzschild horizon behind the shell. Ordering a > R_s > r_h holds at all times.

Constraint and consistency diagnostics along the periodic shell (`phase2d_convergence.json`): Friedmann constraint 8×10⁻¹⁶; Γ_in² identity 9×10⁻¹⁶; p − wσ 2×10⁻¹⁵; mass function R³(ρ/3 + q) + σR²Γ_in − σ²R³/4 against the constant exterior mass, relative 9×10⁻¹⁴ (these three follow algebraically from the definitions and mostly test the code); the **independent** check is the flux conservation law (σR²)′ + p(R²)′ = −R²(ρ + P) cosh β sinh β, relative residual 1×10⁻⁶ (finite-difference limited). The shell return error is 10⁻¹⁴ and converges between rtol 10⁻⁹ and 10⁻¹³.

No-bridge baseline (item 7): the interior history in the coupled run equals the interior-alone run to 2×10⁻¹⁰ (integrator tolerance) for every perturbation. This is not a dynamical finding but the structural fact of Phase 2C §3: a homogeneous interior does not feel the shell, so the bridge cannot change the stability of the cycle for good or ill. The shell has its own stability (table above).

**Causal structure (item 8, derived).** The shell worldline lies in region III (U > 0, V < 0) because R_s > r_h at all times. In the kept exterior, a future-directed causal curve has non-decreasing U and V, so every point of the causal future of the cap has U ≥ U_shell > 0; region I is {U < 0}. Therefore **no future-directed causal curve joins the cap to region I, and none joins region I to the cap** (region I is spacelike separated from region III). Numerically (`phase2d_bridge.json`): 24 radial null rays emitted from the cap centre at different phases all reach the shell, after 4.3–4.8 interior periods (the cap is nearly a full S³, so a radial ray needs several cycles), cross it, and fall along U = const toward the horizon V = 0 into region II (FigD1 panel c). Nothing in this geometry carries matter or signals from the cap to the parent universe.

## 5. Levels

| level | status |
|---|---|
| 1. stable cyclic cosmology | **numerically reproduced independently** (elliptic cycle, C/Cc ∈ (0.346, 0.360), bounded 200+ bounces under ≤ 3% perturbations); not a proof of long-term stability |
| 2. stable cyclic region within a genuine bridge geometry | **reached numerically, for a non-traversable Einstein–Rosen-type bridge** (cap + Schwarzschild throat, flat exterior, conservative shell survives 60 periods at ≤ 1% perturbations); the throat is tiny (r_h ≈ 0.0033 against a cycle radius ≈ 0.95) and hidden behind horizons; the exterior is eternal (includes the white-hole region IV; formation not addressed) |
| 3. wormhole-mediated recycling | **not achieved and excluded here**: static vacuum exterior, region I and the cap are spacelike separated, no NEC-violating matter exists in the model |

## 6. Limits and what is open

- The shell equation of state p = wσ is assumed, not derived from a microphysical action; only w = −1/2 produced periodic shells in the scan. AdS exteriors gave none.
- Spherical symmetry and an exactly homogeneous interior throughout; inhomogeneous perturbations and shell–field interactions are untested and require the PDE solver.
- Irreversibility (Phase 1B) is untouched: the conservative cycle is stable, but any dissipation destroys exact recurrence.
- The Kruskal extension is eternal. A bridge that forms dynamically is a different problem.
- A traversable bridge would need NEC-violating stress-energy at the throat (the manuscript's null-energy proposition); none was introduced.
- Next steps if desired: (i) derive an equation of state for the shell from a field-theory action; (ii) AdS and de Sitter exteriors; (iii) test whether the stable cycle persists with small inhomogeneities (full double-null evolution on both sides of the shell); (iv) clarify the correction to the Phase 1 documents and the v1.3.0 release text (see Phase 2C report §7).
