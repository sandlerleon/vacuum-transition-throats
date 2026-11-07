# Phase 3 (3A–3C): black-hole-mediated cosmological recycling — baseline, theory selection, preliminary feasibility (branch `phase3-recycling`)

Status: **stopped after 3C as instructed. The classical baseline is validated. The most defensible candidate mechanism (spin–torsion bounce, Einstein–Cartan type) gives a homogeneous bounce, but fails two decision gates for proceeding to the full inhomogeneous evolution (3D): it cannot be matched to a vacuum exterior after a horizon forms, and its effective-fluid form is singular and ill-behaved through the bounce.** No daughter universe, recycling or wormhole is claimed. 27 automated checks pass. Phases 1–2 manuscripts, results and public releases are untouched (this branch starts from `phase2-bridge`).

Reproduce (`code/`): `phase3a_baseline.py`, `phase3c_ec_homogeneous.py`, `phase3c_ec_perturb.py`, `phase3c_ec_window.py`, `phase3c_ec_control.py`, `tests/test_all.py`. Units: 8πG = 1 for the scalar solver and Phase 3C; G = 1 for the exact dust solutions in 3A (stated where used).

## 3A. Classical baseline

**Audit of reusable code.**

| component | status for Phase 3 |
|---|---|
| `dnull.py`: spherical Einstein–scalar, double-null, regular centre, trapped/anti-trapped diagnostics, constraint monitor | reusable; scalar field only (no dust, no fluid with pressure); regular initial data on an outgoing cone |
| `collapse_scan.py` results (72 runs), `validation_extra.json` | reusable as the stored baseline for the matter-coupled vacuum-transition model |
| `closed_cycle` / `phase2f_model` / perturbation equations | reusable for homogeneous backgrounds and linear radial perturbations of a fluid with pressure |
| `bridge_junction`, Kruskal/junction code | reusable for exterior matching |
| coupled scalar + dust **inhomogeneous nonlinear** collapse | **does not exist** (gap; needed for 3D) |

**Results** (`phase3a_baseline.json`):
- *Scalar collapse.* Of 72 stored runs, a trapped surface forms in 65. In the 42 runs with a regular, not-expanding initial core, an anti-trapped (expanding) region was created later in **0**; the other 30 started with an expanding false-vacuum core that merely persists. 70 runs ended by non-finite values and 2 by σ > 12. A fresh run (A = 0.05, V_w = 0.8, w = 0.05) forms its first trapped point at u = 0.575 (h = 0.0125) and 0.569 (h = 0.00625), within one grid cell, with no anti-trapped point. (I did not re-validate the solver's Misner–Sharp mass in this session; an attempted extraction was noisy and is not reported.)
- *Oppenheimer–Snyder* (closed dust ball, Schwarzschild exterior, M = 1, χ₀ = 0.9): numerical singularity time 6.536130 against exact 6.536131, horizon time 4.817308 against exact 4.817308, scale-factor error 2×10⁻¹³, Misner–Sharp mass constant to 10⁻¹⁵, Friedmann residual 7×10⁻¹², no re-expansion, density → ∞.
- *Inhomogeneous LTB dust collapse* (96 shells, exact parametric solution): singularity times 3.69 to 8.80 with the centre first; all 96 shells become trapped before their singularity; no shell crossing before collapse; equation-of-motion residual 4×10⁻⁴ and first-integral residual 3×10⁻⁵ (finite-difference limited).

**Gate 3A: passed.** Classical collapse forms horizons and singularities with no daughter region; the expected negative result is reproduced. Limitation: the combined scalar-plus-dust inhomogeneous case is not covered.

## 3B. Candidate-theory comparison and selection

| | Einstein–Cartan spin fluid | LQC-inspired effective model | noncanonical scalar P(X, φ) |
|---|---|---|---|
| defining equations | Einstein–Cartan–Sciama–Kibble action with a Weyssenhoff spin fluid; torsion is algebraic and eliminated | modified Friedmann equation H² = (κ/3)ρ(1 − ρ/ρ_c); no covariant action for inhomogeneous collapse used here | explicit action; bounce needs ρ + P < 0 |
| origin of the bounce | spin–spin contact term: ρ_eff = ρ − (κ/4)s², p_eff = p − (κ/4)s², s = σn (literature form, **not re-derived** here) | the correction term itself (ρ_c is a free Planck-scale parameter) | NEC violation by the kinetic structure |
| conservation | effective stress conserved (Bianchi) | effective, homogeneous | canonical |
| validity | high spin density; spin averaging/polarization assumed | homogeneous isotropic; extensions to inhomogeneous collapse are model-dependent | classical field theory |
| known risks | spin averaging, realistic (Dirac) matter, inhomogeneities, exterior matching | the bounce is built in; no derivation for black-hole interiors | stable NEC violation is reported to be difficult (ghost/gradient instabilities) — a literature statement **not verified** here |
| arbitrary bounce function? | no: one parameter B = κσ²/4 set by the spin | yes (ρ_c) | no, but NEC violation is imposed by the model |

**Selection: Einstein–Cartan**, as the instructions suggested, because the high-density correction comes from a specified action with one physically set parameter, not a prescribed bounce. A result derived in 3C reduces the choice: for dust, the LQC-type modified Friedmann equation is, as a conserved effective fluid, **the same fluid** (ρ_eff = ρ − ρ²/ρ_c, p_eff = −ρ²/ρ_c = −Bn² with B = 1/ρ_c; checked to 10⁻⁹), so every feasibility test below applies to both. The noncanonical scalar was not implemented.

## 3C. Preliminary feasibility (homogeneous and linear)

**Effective fluid.** With n = 3C/a³ the conserved number density: ρ_eff = n − Bn², p_eff = −Bn², B ≥ 0. Internal consistency: the thermodynamic identity p = n dρ/dn − ρ holds to 3×10⁻¹², and continuity ρ̇ + 3H(ρ + p) = 0 and the Friedmann constraint hold along the solution (residuals 5×10⁻⁴ finite-difference, 1×10⁻¹¹).

**Homogeneous bounce: exists.** Closed ball: ȧ² = C/a − 3BC²/a⁴ − 1, bounce and turnaround at the two roots of Ca³ − a⁴ − 3BC² = 0 (they exist iff B ≤ 9C²/256). For C = 0.01, B = 2×10⁻⁷: bounce a = 0.00195, turnaround a = 0.00994, period 0.0317, return error 9×10⁻¹⁶. The bounce occurs at n = 4.0×10⁶ > n* = 1/(2B) = 2.5×10⁶, where ρ_eff + p_eff = n(1 − 2Bn) < 0: the effective NEC is violated only there, as a bounce requires.

**Matching to a vacuum exterior: fails.** Exterior Schwarzschild with F = 1 + Ṙ² − m/R, m = C sin³χ₀ (horizon R = m). At a bounce Ṙ = 0, so a vacuum exterior needs R_bounce ≥ m: **a shell cannot bounce inside a vacuum horizon** (inside it every timelike curve has R decreasing). Over 35 combinations of χ₀ and B, in 29 the boundary bounces at R < m (impossible), in 6 no horizon forms (the spin term wins before trapping, so there is no black hole), and in every case the interior effective mass R³ρ_eff/3 varies over the cycle by 24% to 97% of the exterior mass, so it cannot equal the constant exterior mass. A homogeneous ball also has p_eff ≠ 0 at its surface, which a vacuum exterior cannot have without a shell. So a regular transition after horizon formation needs a non-vacuum region outside the bouncing core (an inhomogeneous problem), not a homogeneous ball.

**The effective fluid is pathological through the bounce.** Sound speed c_s² = dp/dρ = −2Bn/(1 − 2Bn): **negative (imaginary sound speed) for n < n***, diverging at n*, and **greater than 1 (superluminal) for n > n*** (minimum 1.39 above n* on the grid). The perturbation equations of the fluid are singular at n*, where ρ + P = 0: a monodromy through the bounce gives determinants ≈ 1.8×10¹⁰ instead of 1, so numbers through that point are not meaningful. In the regular window n ≤ 0.7 n* the transfer matrices have det = 1 for n = 3 and 6 and show growth e^{10.5}–e^{15.6}. A **plain dust control** (same contraction leg, B = 0) already grows by e^{8}–e^{12}, so ordinary collapse is violently unstable; the spin term adds a k-dependent excess on that leg, ln-growth difference (B > 0 minus B = 0): −0.2 to −0.5 at n = 3, +0.8 to +1.9 at n = 6, +2.4 to +4.6 at n = 12, increasing with n and with the closeness to n*, as an imaginary sound speed predicts. The WKB exponent ∫k|c_s|/a dt predicts a larger excess (about 2× at n = 12) and growth ∝ k beyond, which was not verified numerically for n ≥ 30 (loss of precision). I therefore claim the sign and trend of the extra instability, not its UV behaviour.

## Decision gates for 3D

| gate | criterion | result |
|---|---|---|
| G1 | a homogeneous regular continuation exists | **pass** (bounce, cycle) |
| G2 | derived from an action, not a prescribed bounce | **pass** for EC (one parameter; coefficient taken from the literature, not re-derived) |
| G3 | consistent matching of the transition region to the exterior (vacuum or a derived non-vacuum) | **fail** for a homogeneous ball: R_bounce < m in 29/35, mass mismatch (24–97%) in all |
| G4 | the effective equations are well-posed through the transition (real c_s² ≤ 1, regular perturbation equations) | **fail** for the barotropic effective fluid: c_s² < 0 below n*, singular at n*, superluminal above |
| G5 | classical baseline reproduced | **pass** |

**Verdict: do not proceed to the expensive 3D evolution with this model.** The failures concern the effective perfect-fluid description. They do not show that spin–torsion gravity cannot produce a regular collapse: the Weyssenhoff fluid has spin dynamics and polarization that a barotropic fluid hides, and realistic inhomogeneous collapse leaves the outer layers as dust.

## Demonstrated, assumed, unresolved

- **Demonstrated (numerically or exactly):** the classical baselines (§3A); thermodynamic and conservation consistency of the effective fluid; the equivalence of the LQC-type and spin-type effective fluids for dust; a homogeneous bounce and cycle; the impossibility of a bounce inside a vacuum horizon; the sign structure of c_s²; the singularity of the fluid perturbation equations at n*; the extra k-growing instability at n ≤ 12.
- **Assumed:** the coefficient and form of the spin term (ρ_eff = ρ − κs²/4, s = σn); a barotropic perfect-fluid completion; radial linear perturbations in longitudinal gauge; the effective cosmological constant of the Phase 2 machinery (negligible at a ≤ 0.01).
- **Unresolved:** whether a non-barotropic spin fluid (evolving spin density, Dirac field with the Hehl–Datta contact term) avoids the effective-fluid pathology; whether an inhomogeneous collapse with a bouncing core and an infalling dust envelope can be globally consistent; formation of an expanding region with a verified causal connection to the infalling matter (3E); matter-energy transfer (3F). The literature claims in the 3B table were not checked in this session.

## Next step if you want to continue

Before any 3D code: (i) re-derive the spin-fluid effective stress from the Einstein–Cartan field equations for a Weyssenhoff fluid with an evolving spin tensor and check whether the barotropic completion is justified or whether spin propagation modifies c_s²; (ii) consider a different matter completion (e.g. a Dirac field), or a different candidate (the noncanonical scalar, after checking its stability constraints in the literature); (iii) only then build the spherical horizon-penetrating solver with a dust envelope. Otherwise stop here and report: Phase 3A–3C is a negative-with-structure result that can be written up with the earlier negative results.
