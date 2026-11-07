# Phase 2F: attempts to stabilize the primordial nontraversable cyclic region (branch `phase2-bridge`)

Status: **no parameter region was found in which the homogeneous cycle is also stable against inhomogeneous perturbations.** Every attempted change (weaker coupling, a different hilltop curvature, a different hilltop shape, a fluid with pressure) leaves growth per cycle of at least 8.6 for n = 3 and at least 3.6×10³ for the worst mode. Two analytic criteria for the instability are derived and validated. The shell–bulk exchange could not be derived from an action. 25 automated checks pass. The CQG manuscript was not touched; the Phase 2E manuscript is revised separately (§7).

Reproduce (`code/`): `phase2f_derive.py` (symbolic, ~minutes), `phase2f_model.py` (validated against Phase 2E), `phase2f_continuation.py up|down`, `phase2f_continuation_mu.py up|down`, `phase2f_pressure.py`, `phase2f_potential.py`, `phase2f_worstcase.py`, `phase2f_attribution.py`, `phase2f_analysis.py`, `make_phase2f_figure.py` (→ `figures/FigF1_phase2f.png`), `tests/test_all.py`.

## 1. Method

- **Generalised derivation** (`phase2f_derive.py`): the Phase 2E linear perturbation equations re-derived for a polynomial potential V = 3 + a₂φ² + a₄φ⁴ + a₆φ⁶ and a coupled fluid **with pressure**: particle density n, ρ_f = A(φ)n + Kn^γ/(γ−1), P = Kn^γ (γ = 5/3), A = 1 + φ²/2M². Field equation □φ = V′ + A′n; fluid number conservation and (ρ + P)u^ν∇_ν u_μ = −h_μ^ν(∂_νP + A′n∂_νφ), which follows from total energy–momentum conservation. The traceless equation again gives Ψ = Φ. Background: a″ = −(a/6)(2φ̇² − 2V + An + e + 3P_f), φ″ = −3Hφ̇ − V′ − nφ/M². **Validation:** at K = a₆ = 0 it reproduces the Phase 2E multipliers (n = 3: 14.3398; 5: 88.56; 10: 1061.85; 60: 154 820).
- **Continuation, not scans.** Each family is followed along the elliptic orbit at **fixed return-map trace** (−1.7043; rotation number 0.41), solving for φ₀ (time-symmetric condition) and adjusting C/C_c by secant at every parameter step. At every point the full homogeneous return map and the inhomogeneous growth for n = 3, 6, 12, 30, 60 (and up to 240–600 where stated) are computed.

## 2. Results

| change | range followed | n = 3 | n = 12 | n = 60 | notes |
|---|---|---|---|---|---|
| baseline (μ = 20, ρ_crit = 0.2 V_f) | – | 14.3 | 2.05×10³ | 1.55×10⁵ | |
| weaker/stronger dust–scalar coupling ρ_crit/V_f | 0.0086 → 0.51 (orbit lost beyond both ends) | 14.3–16.4 (15.9 at 0.0086) | 2.0×10³–2.4×10³ | 1.55×10⁵ → 8.6×10⁴ (5.0×10⁵ at 0.0086) | C/C_c retunes 0.02 → 0.68; fifth-force rate at the neck 24.3 → 15.6 only |
| hilltop scale μ | 20 → 307 (branch jumps beyond ~320; the raw continuation file also contains the post-jump points, `phase2f_muscan_up_rc0.2_valid.json` the valid ones) | 14.3 → 8.6 | 2.0×10³ → 183 | 1.55×10⁵ → 3.5×10³ | worst over n: 1.55×10⁵ (μ=20), 4.3×10⁴ (40), 1.4×10⁴ (78), 4.5×10³ (163), 3.6×10³ (229); min over n ≥ 8.7 |
| hilltop shape (sextic, same top curvature and well depth) | v/v₀ 0.98 → 1.74 | 12.3–16.7 | > 10³ | 3.7×10⁴–3.9×10⁵ | orbit lost at v/v₀ = 0.96 and 1.78 |
| fluid pressure (γ = 5/3) | K = 0 → 0.08 (c_s² to 0.038, P/ρ_d 0.024); lost at K = 0.16 | 14.4 | 1.99×10³ | 1.5×10⁵ → 3.6×10⁴ | n = 120: 2.7×10⁵ → 395; n = 240: 3.2×10⁵ → 1.0 |

Coupled fluid with pressure is the only change that removes any modes, and only those with n ≳ 85. The homogeneous elliptic cycle ceases to exist (continuation lost) before c_s² reaches 0.08.

## 3. Why: two mechanisms, two criteria

Both are derived analytically and compared with the Floquet numerics (`phase2f_analysis.py`).

**Fifth force (coupled fluid).** In the quasi-static limit the density contrast obeys D̈ = (g² − c_s²k²/a²)D with g² = 6Cφ₀²/(M²a³(2M² + φ₀²)). Mode n grows when c_s²k²/a² < g², i.e. **n < n_c = g a/c_s**, with exponent S₅(n) = ∫√(g² − c_s²k²/a²) dt. For dust (c_s = 0) every mode grows, at the k-independent exponent ∫g dt = 12.57 (rate 18–23), the large-n limit of the numerics. **Validation:** with pressure, S₅ against ln G: K = 0.04, n = 120: 8.88 vs 9.11; n = 240: 2.68 vs 2.58; K = 0.08, n = 120: 5.64 vs 5.98; n = 240: 0.00 vs 0.00. With g ≈ 18 and a ≈ 0.94 the cut is n_c ≈ 85 at K = 0.08, as found. It fails at low n (it predicts 12.8 where the exponent is 2.7) because the quasi-static assumption breaks when k is small.

**Hilltop (tachyonic) instability.** The scalar sector alone, δφ″ + 3Hδφ′ + (k²/a² + V″(φ₀) + n/M²)δφ = 0, grows when the effective mass is negative, which happens near the top of the potential for k < a μ: S_h(n) = ∫√(max(−(k²/a² + V″ + n/M²), 0)) dt. The scalar-only Floquet multiplier (`phase2f_attribution.json`) is 2.7 (n = 3) and 1 (n ≥ 12) at μ = 20 but 15, 34, 109, 548, 1.2×10³ for n = 3, 6, 12, 30, 60 at μ = 208: at large μ the hilltop mechanism accounts for most of the log-growth (7.1 of 7.9 at n = 60). The product μT stays at 12–16 along the whole μ family and ln G(n = 3) stays at 2.2–2.7, so the hilltop exponent is scale-free: raising μ shortens T but not μT. The sum S₅ + S_h overestimates the true exponent (the mechanisms are not additive), so these are mechanism-level estimates, accurate only in the pressure-cut regime above.

**Consequence.** A stable cycle needs both the fifth force and the hilltop band closed. Pressure can close the first only above n_c ≈ g a/c_s, and c_s² cannot exceed about 0.04 before the cycle itself is destroyed; the hilltop band k < a μ is not affected by anything tried. Large μ weakens the fifth force (g T falls from 10.7 to 0.5) but widens the tachyonic band to n ≲ μ. Hence the **worst-case growth never falls below 3.6×10³ per cycle**, and the low harmonics keep G ≳ 9. A 10⁻⁶ perturbation becomes nonlinear in about two cycles (worst mode) to five cycles (n = 3).

## 4. Does any region satisfy all requirements?

No. Within what was explored there are parameters with nonsingular repeated cold cycles, bounded homogeneous perturbations, consistent constraints, a positive-tension shell satisfying the NEC and a valid non-traversable two-ended geometry (Phase 2D, at μ = 20, ρ_crit = 0.2 V_f), but never with bounded or substantially reduced inhomogeneous perturbations.

## 5. Shell–bulk interaction (item 6): not derived

For a shell that moves through the homogeneous interior, energy–momentum conservation across it requires the flux law (σR²)′ + p(R²)′ = −R²(ρ + P)cosh β sinh β. It vanishes only for β = 0 (comoving shell, which is NEC-violating, Phase 2C) or ρ + P = 0. Decomposing the gross exchange over one cycle of the conservative periodic shell (`phase2f_analysis.json`): **dust 0.494 (92%), scalar field 0.042 (8%)**, nets 3×10⁻¹⁴ and 2×10⁻¹⁵.

The dust part admits a bookkeeping model: the shell absorbs dust that leaves the cap as it moves inward and emits comoving dust when it moves outward, so the exterior stays vacuum. That is a perfect absorber/emitter boundary condition, not an action. The scalar part would require a boundary condition for a homogeneous classical field that exchanges energy with the shell without exciting inhomogeneous modes; no such condition is available in the homogeneous reduction. I could not derive a shell–bulk coupling from an action that reproduces the exchange, and I did not introduce one by hand. The shell remains an effective junction model.

## 6. Limits of the search

Linear, spherical (radial harmonics n ≥ 3 on the full S³; n = 2 not covered; no shell boundary conditions); one fluid model (polytrope γ = 5/3); the return-map trace fixed at −1.7043 (rotation number 0.41) so other elliptic orbits and windows are not covered; q = 1; μ from 20 to 307 and ρ_crit/V_f from 0.0086 to 0.51 (orbits lost beyond); the sextic family only. Nonlinear spherical evolution was **not** run because no promising region exists (item 7). The WKB criteria were validated numerically but are not rigorous bounds; the low-n exponent (ln G ≈ 2.2–2.7, G ≳ 8.7) is an observed invariant, not a derived one.

## 7. Consequences for the manuscript

The new results materially strengthen the negative conclusion and supply the "general instability criterion" route: the manuscript (Phase 2E draft) gains a subsection and an abstract sentence on attempts to stabilize (§5.3), and its discussion no longer leaves the weak-coupling and pressure routes open. The title and construction are unchanged. It is not submitted.
