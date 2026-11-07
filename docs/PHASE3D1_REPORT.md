# Phase 3D.1: quantum-wave propagation through a hypothetical recycling throat (branch `phase3-recycling`)

Status: **stopped before any coupled quantum-gravity or axisymmetric solver, as instructed. Result: waves do propagate through, and are filtered by, the benchmark throat; the benchmark is causally passable; but its required stress-energy is not physically supported, and nothing here identifies a self-consistent recycling mechanism.** 32 automated checks pass. Phase 1–2 manuscripts and public releases untouched.

Reproduce (`code/`): `phase3d1_wave.py` → `phase3d1_wave.json` (≈9 min), `tests/test_all.py`. Units G = c = ħ = 1, κ = 8π.

## Part A. The benchmark geometry

Simpson–Visser "black-bounce" family, ds² = −f dt² + dr²/f + R²dΩ², R² = r² + a², f = 1 − 2m/R, r ∈ (−∞, ∞). One parameter spans the three classes the instructions ask to separate:

| a/m | class | horizon | radial E=1 geodesic crosses r = 0 (checked) |
|---|---|---|---|
| > 2 | traversable two-way wormhole | none (f_min = 1 − 2m/a > 0) | yes (a/m = 3) |
| = 2 | extremal one-way transition | degenerate, r_h = 0 | yes |
| < 2 | regular black hole; interior bounces through a white-hole region into a second asymptotic region (one-way: r is timelike inside) | r_h = √(4m² − a²) | yes (a/m = 1; proper time +50 → −50 finite, 334) |
| 0 | Schwarzschild | – | – |

Null expansions: trapped iff f < 0, i.e. only for a < 2m. Curvature is regular: Kretschmann at the throat = 0.066 (a/m = 3), 16 (a/m = 1), 1472 (a/m = 0.5) in units m = 1, finite for all r. **This is a mathematical benchmark. It is not claimed to solve any physically admissible field equations**; its source is computed in Part C. It is also not a model of collapse: it is an eternal, static geometry.

## Part B. Quantum diffraction and transmission (test scalar field)

(□ − μ² − ξR)Φ = 0 → u″ + (ω² − V_ℓ)u = 0, dr*/dr = 1/f, V_ℓ = f[ℓ(ℓ+1)/R² + μ² + ξ𝓡] + (1/R)d²R/dr*². Transmission T (wormhole) or greybody Γ (black bounce) from a Wronskian-normalized integration. **Flux conservation T + R = 1 holds to 8×10⁻⁷ worst case over the 108 massless runs; T is converged in the box size to 2×10⁻⁷ (ωL = 1000, 2000, 4000) and in tolerance.** The massive runs carry the long-range Coulomb tail, with residuals 10⁻⁵ to 5×10⁻⁴ (tail projection error, stated).

Selected results (ℓ = 0, massless):

| geometry | T(ω=0.05) | T(0.1) | T(0.2) | T(0.4) |
|---|---|---|---|---|
| wormhole a/m = 3 | 0.0070 | 0.15 | 0.996 | 1 |
| wormhole a/m = 6 | 0.064 | 0.57 | 0.996 | 1 |
| Ellis (m = 0, a = 3) | 0.010 | 0.051 | 0.295 | 0.899 |
| black bounce a/m = 1.9 | 0.069 | 0.43 | 0.96 | 1 |
| black bounce a/m = 1 | 0.061 | 0.35 | 0.93 | 1 |
| Schwarzschild | 0.059 | 0.32 | 0.91 | 1 |

- **The throat is a high-pass filter whose cutoff scales with 1/(throat radius): ωa ≈ 0.3–0.6 for the wormhole.** Low-frequency transmission is strongly suppressed relative to the black-hole greybody factor (0.007 vs 0.059 at ω = 0.05 for a/m = 3). Fitted low-ω exponents: ℓ = 0: 3.2 (wormhole), 2.2 (Ellis); ℓ = 1: ≈ 6–7. Higher ℓ is suppressed far more (T ≈ 10⁻⁶ or less at ω = 0.08): **the throat redistributes outgoing flux toward the lowest partial waves**, as the slit analogy suggests, but the "diffraction" is the centrifugal-barrier filter of the geometry, not a force.
- **Compton wavelength vs throat radius:** for μa = 0.3 the field is reflected just above threshold (T = 0.034 at ω = 1.05μ, 0.76 at 1.5μ, 0.997 at 3μ); at μa = 0.9 T ≈ 0.99 already at threshold (gravitational attraction of the throat region).
- **Black bounce vs Schwarzschild:** greybody factors differ by < 30% for a/m ≤ 1.9 (T ratio 1.06 at a/m = 1, ω = 0.1): from outside, the bounce interior hardly changes absorption. This is "transmission" into the horizon, not into a second universe observable from outside; what happens after the horizon is a statement about the interior geometry only.
- Transmission only characterizes propagation through a geometry that is assumed. **It does not create, open or stabilize a throat.**

## Part C. Gravitational backreaction requirements

The Einstein tensor of the benchmark (sympy) gives, in closed form, for the radial null direction **T_kk = −2a²/(κR⁴) < 0 at every r, for every m and every a > 0** (radial null energy condition violated everywhere; checked symbolically and on grids). The throat has ρ = 0.0015, p_r = −0.0044, p_t = 0.0030 (a/m = 3, units m = 1). The volume integral ∫(ρ + p_r)dV is **−7.0** (a/m = 3) and **−16.5** (a/m = 6) in Planck units (negative, growing ∝ a): a net negative-energy requirement of order a/G.

- **Classical test-field waves cannot supply this.** A minimally coupled scalar has T_kk = (k·∇Φ)² ≥ 0; non-minimal coupling would need Φ² ≳ c/(κ|ξ|), a Planckian amplitude, which breaks the test-field approximation by construction.
- **Quantum vacuum.** With a Ford–Roman-type bound (flat-space free massless scalar, Lorentzian sampling, τ₀ ≲ a) the required density needs **N ≳ 5.6 (a/ℓ_P)² field species**: 2×10⁷⁰ for a = 1 m, 2×10⁷⁶ for 1 km, 10⁸⁸ for a solar radius. So free-field vacuum polarization supports a throat only near the Planck scale. This is an order-of-magnitude estimate from a flat-space inequality, not a semiclassical calculation. **⟨T_μν⟩_ren on this background in a physically specified state was not computed**; so I cannot exclude a state-specific mechanism, only state that none has been found, and that the burden is O((a/ℓ_P)²) species or Planck-scale throats.
- Known sources of this specific metric (phantom scalar plus non-linear electrodynamics) are imposed, not derived from a quantum state.

## Part D. Relation to Einstein–Cartan–Dirac (3D.0)

The wave analysis uses a prescribed background; the ECD bounce does not enter it, and nothing computed here modifies the high-density collapse. The Dirac contact interaction acts only at k_F ~ M_Pl; the transmitted-wave stresses are far below that. The first-order Hartree–Fock model was not extrapolated. A Dirac-field scattering calculation on the benchmark (spinor radial equations) was not done.

## Decision gate

| question | result |
|---|---|
| does quantum-wave transmission occur in the benchmark? | **yes**: T → 1 above ωa ≈ 0.5, high-pass filtered, higher ℓ suppressed; flux conserved to 8×10⁻⁷ |
| is the geometry causally passable? | **yes** (two-way for a > 2m; one-way through a horizon for a < 2m; geodesic checks) |
| is the required gravitational stress-energy physically defensible? | **no, as specified**: radial NEC violated everywhere, ∫(ρ+p_r)dV < 0 ∝ a; free-field quantum support needs ~5.6 (a/ℓ_P)² species; no derived quantum state |
| has a self-consistent mechanism for matter recycling been identified? | **no** |

## Honest limits

Static, eternal benchmark (no collapse, no formation process); test scalar field only (no Dirac or gauge fields); ξ = 0 in the scattering table; the quantum-inequality estimate is flat-space and order-of-magnitude; no ⟨T_μν⟩_ren or semiclassical backreaction solved; massive-field residuals up to 5×10⁻⁴.

## What would be the next specific step

(i) A state-specific ⟨T_μν⟩_ren for a conformally coupled field on the wormhole branch (a > 2m), e.g. in the Page/DeWitt-Schwinger approximation, to see whether the sign and the magnitude of T_kk can be matched at any a ≫ ℓ_P; (ii) a time-dependent version (formation of a throat from collapse) rather than the static benchmark.
