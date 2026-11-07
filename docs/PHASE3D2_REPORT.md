# Phase 3D.2: free-fall transport through the black-bounce benchmark (branch `phase3-recycling`)

Status: **stopped before any collapse simulation, as instructed. Question A (can matter travel through the geometry?) is answered yes, kinematically and with finite tides. Question B (can physically admissible gravity and matter create the geometry?) is not answered and is constrained by the Penrose theorem.** 33 automated checks pass.

Reproduce (`code/`): `phase3d2_freefall.py` → `phase3d2_freefall.json` (≈70 s), `tests/test_all.py`. G = c = 1, κ = 8π, lengths in units of m. Benchmark: Simpson–Visser, ds² = −f dt² + dr²/f + R²dΩ², R² = r² + a², f = 1 − 2m/R; **prescribed, not shown to form** (Phase 3D.1).

## Part A. Geodesic transport

Timelike geodesics: (dr/dτ)² = E² − U(r), U = f(1 + L²/R²), integrated as r̈ = −U′/2 in proper time (regular across horizons and turning points). Conserved E² drifts by < 4×10⁻¹¹.

- **Radial free fall crosses the minimum-radius surface.** Released from rest at R₀ = 20m (E² = 0.9): black bounce a/m = 1 reaches r = 0 at τ = 99.95; wormhole a/m = 3 at τ = 102.4. Inside the horizon U < 0, so there is no turning point whatever E is; outside, U rises monotonically with R, so every radial infall reaches the horizon.
- **Bound fall does not stop at the throat.** It emerges, reaches the mirror radius r = −19.97 (R = 20.00), turns around and returns: the worldline is periodic (turning points at τ ≈ 200, 400, 600, 800 for a/m = 1; 205, 410, 614, 819 for a/m = 3). Unbound (E ≥ 1) fall continues to r → −∞ with no turning point.
- **Crossing criterion** (derived and checked by direct integration at 0.98 and 1.02 of the critical value, all 20 cases correct): crossing iff E² > max_{R∈[a,R₀]} U(R).

## Part B. ISS-style orbital comparison

- **Angular momentum never helps.** ∂U/∂(L²) = f/R² > 0 outside the horizon, so the set of crossing orbits is {L ≤ L_c(E)}: radial infall is the best case. Inside the horizon f < 0 and there are no turning points.
- **Critical L for crossing:** for a/m = 0.5, 1, 1.9, 3 it equals the **Schwarzschild barrier value** (4.0000 at E = 1; 5.3728 at E = 1.2), i.e. the throat does not alter the plunge condition when it lies behind the horizon or behind the usual peak. For the large-throat wormhole a/m = 6 it is throat-limited: L_c = a√(E²/f(a) − 1) = 4.243 (E = 1), 6.462 (E = 1.2).
- **Circular geodesics.** Exterior stable/unstable orbits exist as in Schwarzschild for a/m ≤ 3 (4 stable, 4 unstable among L ∈ {1,...,8}); for a/m = 6 only stable ones (no photon-sphere-like peak outside the throat). **The throat r = 0 is itself a circular geodesic of a static wormhole (a > 2m)**: stable (potential minimum) for all L when a ≤ 3m; for a > 3m stable only for L² < m a²/(a − 3m) (L < 3.46 at a/m = 6) and a barrier (unstable) beyond. So ISS-style orbiting matter can sit **in** the wormhole throat or be blocked by it; it never enters preferentially. Orbital motion should not be confused with traversability.

## Part C. Tidal forces and matter transport

In the freely falling frame of a radial geodesic (symbolic; radial K₁₁ = R_trtr, transverse K_θ from R_θtθt and R_θrθr):

- Both are **finite at the horizon** (K_θ(r_h) = 3/32 − E²/16 for a = m = 1) and **finite at the throat**: max |K₁₁| = 1.000 and max |K_θ| = 1.900 in units 1/m² for a/m = 1, scaling as ≈ m/a³ and ≈ 2m/a³ (Schwarzschild at R = a would give 2m/a³): **a/m = 0.25: 64 and 126; 0.5: 8 and 15.6; 1.5: 0.30 and 0.55; 1.9: 0.15 and 0.26.** Curvature invariants and observer-frame tides are both finite; the strength is set by the throat radius a, not by m.
- **Jacobi fields.** The transverse deviation of neighbouring radial geodesics equals the proper shell size, ξ_θ ∝ R(τ), to 1.2×10⁻¹¹ (validates K_θ and the integration); it is never zero (R ≥ 1.000). The radial deviation grows by a factor 31 in one passage and has **one sign change** (a radial caustic of the family, finite, not a divergence).
- **For astrophysical masses:** the tide at the throat is ≈ c²m/a³ ~ (1/m²)·(c²) for a ~ m: for m = 1 M_sun (1.5 km) this is ~10¹¹ s⁻²: destructive for anything extended, but finite and smaller by (M/M_sun)⁻² for supermassive holes (order 10⁻⁶ s⁻² at 4×10⁶ M_sun... independent of the observer being in free fall).

## Part D. Global causal structure

- Region sequence along the bound radial fall (a/m = 1): exterior (r>0) → **trapped** (τ = 98.5, R decreasing) → **anti-trapped** (τ = 100.0, R increasing; r = 0 is the marginal throat) → exterior (r<0) (τ = 101.4) → turning point → **trapped** (τ = 298.4) → anti-trapped → exterior (r>0) (τ = 301.3): a chain of copies. Inside |r| < r_h = √(4m² − a²), r is timelike and dr/dτ has a fixed sign (U < 0), so **the passage is one-way**: no future-directed curve returns to the first exterior. For r<0 the radial null expansions are both positive (R increasing): this is the anti-trapped, white-hole-like region (a Kantowski–Sachs bounce).
- **Is it a daughter universe?** The geodesic congruence in the anti-trapped part has positive volume expansion Θ (range 0.003–1.50 for E² = f(20m); up to 10.4 for E = 0; 0.004–1.52 for E = 1). But the region leads to **another asymptotically flat exterior** with its own horizons, not to a closed cosmological expanding universe; matter is conserved along the worldlines and no new hot big bang appears. An expanding areal-radius region is not by itself a cosmological daughter.

## Part E. Physical feasibility

- The stress-energy is a property of the metric, **independent of E, L and of the observer**: radial null T_kk = −2a²/(κR⁴) < 0 at every r, ∫T_kk R² dr = −a/4 (G = 1). Free fall changes none of this.
- **Penrose theorem.** The maximal extension contains closed trapped surfaces (r > 0 inside the horizon), is geodesically complete along the worldlines computed, and has a non-compact Cauchy surface. By the Penrose singularity theorem, geodesic completeness with a trapped surface requires NEC violation, consistent with T_kk < 0. **Formation from NEC-satisfying classical collapsing matter is therefore excluded by the same theorem**; any formation channel must violate the NEC (quantum, torsion/ECD at Planck densities, or unknown). The Phase 3D.0 ECD bounce satisfies none of the requirements for this (barotropic mean field, unstable).
- Test-particle transport cannot establish matter recycling: back-reaction of the infalling matter (which satisfies the NEC) was not computed.

## Decision gate

| question | result |
|---|---|
| does radial free fall cross the minimum-radius surface? | **yes** for every radial infall (bound or unbound), E² drift < 4×10⁻¹¹ |
| does angular momentum help or hinder? | **hinders only**: crossing iff L ≤ L_c(E) (Schwarzschild barrier for a/m ≤ 3, throat-limited at a/m = 6); the throat is a stable or blocking circular orbit, not a funnel |
| are tidal forces finite and tolerable? | **finite** (≈ m/a³ to 2m/a³, regular at the horizon and throat); tolerable only for objects far smaller than a; destructive for extended bodies at m ~ M_sun |
| entry to a genuinely expanding causal region? | **anti-trapped region: yes** (Θ > 0), **cosmological daughter universe: no**; the region leads to a second exterior; one-way: yes |
| evidence that the spacetime can form dynamically? | **none**; NEC violation required (Penrose); 3D.1: no physically supported source; 3D.0: ECD bounce unstable |

## Limits

Static, prescribed geometry; test particles only (no back-reaction, no collisions, no pressure); radial tidal analysis and Jacobi fields for L = 0 only; the circular-orbit table samples six angular momenta; the supermassive scaling is an order-of-magnitude statement; global structure shown by region classification along a worldline and by the sign of dr/dτ, not by a full conformal diagram.

## Where this leaves the programme

Question A is answered favorably but conditionally: **if** a Simpson–Visser-like geometry exists, matter can follow freely falling worldlines through it with finite tides, one way, into an anti-trapped region. Question B is the obstacle and the Penrose theorem sharpens it: the source must violate the NEC. The next specific step is not more kinematics but a source: a dynamical model (even spherically symmetric, with a consistent NEC-violating but physically motivated field, e.g. a semiclassical or torsion-based effective stress) in which a trapped surface forms from ordinary matter and the interior evolves toward a bounce, and which then asks whether the NEC violation is generated rather than imposed.
