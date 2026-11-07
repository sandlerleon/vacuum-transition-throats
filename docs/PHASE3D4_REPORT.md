# Phase 3D.4: proper-time and instability-window diagnostic (branch `phase3-recycling`)

Status: **stopped before any coupled collapse evolution, as instructed. Verdict: a kinematic proper-time window exists on the prescribed benchmark, but it does not carry over to the instability that matters. Fast transit shortens proper time, not the Killing time over which a global mode grows, which has a light-crossing floor; no growth rate was derived for the only source the benchmark requires (a ghost); gravitational time dilation does not weaken local growth, remove the NEC violation or supply the source.** 35 automated checks pass.

Reproduce (`code/`): `phase3d4_propertime.py` → `phase3d4_propertime.json` (≈20 s), `tests/test_all.py`. G = c = 1, lengths and times in units of m unless stated. Benchmark: Simpson–Visser; **prescribed, not shown to form**.

## Parts A–B. Invariant high-curvature region and proper-time transit

Curvature from the orthonormal Riemann components (validated: reproduces the 3D.1 values K(0) = 0.0658, 16, 1472 for a/m = 3, 1, 0.5 to 10⁻⁹ and the Schwarzschild limit 48m²/r⁶). K peaks at the throat for every a/m tested; the throat curvature length K^{−1/4} is 0.16, 0.50, 1.26, 1.97, 3.6 m for a/m = 0.5, 1, 1.9, 3, 6. The region is defined invariantly as **K ≥ ε·K(0)** (not by a coordinate radius).

Proper time through it, radial, E = 1 / E = 5:

| a/m | ε = 0.5 | ε = 0.1 | ε = 0.01 |
|---|---|---|---|
| 0.5 | 0.17 / 0.06 | 0.38 / 0.14 | 1.2 / 0.39 |
| 1 | 0.51 / 0.14 | 1.69 / 0.44 | 4.0 / 0.91 |
| 1.9 | 2.0 / 0.41 | 5.5 / 1.0 | 11 / 1.8 |
| 3 | 4.1 / 0.67 | 9.5 / 1.4 | 19 / 2.5 |
| 6 | 10 / 1.2 | 22 / 2.4 | 42 / 4.2 |

- **Δτ grows with the throat radius, ≈ a (not with m), and falls roughly as 1/E at large E.** Horizon-regular: it is a proper-time integral dτ = dr/√(E² − U).
- **Angular momentum is not monotonic in its effect on transit time.** For the black bounce (a/m = 1, ε = 0.1) L = 0.9 L_c *shortens* the interior transit (1.69 → 0.78, factor 0.46) because inside the horizon f < 0 and the L²f/R² term is negative; for the wormhole (a/m = 3) it lengthens it (9.5 → 20, ×2.1). It still never helps crossing from outside (3D.2).
- **Clock facts:** from rest at 20m, a/m = 1: proper time 98.5 m to the horizon and 99.95 m to the throat, regular. The Killing time to reach r_h + δ is 123, 134, 145 m for δ = 10⁻², 10⁻⁴, 10⁻⁶: logarithmic divergence, a coordinate artifact of t, not a physical freezing.

## Part C. Redshift from four-momenta (exact)

For a radial infaller with Killing energy E, a radially ingoing photon with Killing frequency ω is seen at ω_obs = ω/(E + √(E² − f)); verified against explicit four-vectors (null residual 0, frequency match 10⁻¹²). This is finite and regular across the horizon (a/m = 1: 0.795 at r = 30m, 0.514 at r = 2m, 0.414 at the throat). The signal sent *out* by the infaller and received at infinity is redshifted by f/(E + √(E² − f)) → 0 at the horizon (0.742 at r = 30m, 0.054 at r = 2m): **the "freezing" is a property of signals received outside, not of the infaller's clock.** Static observers (clock ratio √f) exist only where f > 0; there is none inside |r| < r_h. Nothing in the geometry defines a ratio between a clock in the first exterior and one in the second; only proper time along matter worldlines is meaningful (the chain of regions takes ≈ 3 m of proper time between horizon crossings).

## Part D. Growth rates: what is derived and what is not

- **No growth rate is derived for the only source this metric requires** (the forced phantom scalar, 3D.3). The literature reports an unstable radial mode for static Ellis–Bronnikov-type wormholes (secondary sources, not recomputed). I do not invent one.
- **Critical rate (derived, model-independent):** matter survives a perturbation seeded at entry if γ·Δ ≤ A_max = ln(δ_nl/δ₀), where Δ is the relevant elapsed time. For proper time and ε = 0.1, radial: **γ_c = A_max/Δτ**: a/m = 1: 4.1 m⁻¹ (E = 1) to 16 m⁻¹ (E = 5) for A_max = 6.9 (δ₀ = 10⁻³); in seconds, ×1/(4.9 μs·M/M_sun), i.e. 8×10⁵ s⁻¹ for 1 M_sun. A throat-scale instability with γ ~ 1/a would therefore give A = O(1)–O(few): neither automatically safe nor automatically fatal.
- **One computed instability model (a toy, not the source's instability):** a non-minimally coupled test scalar (ξR coupling) on the wormhole branch has a bound state with Killing growth rate γ_t = √(−Λ₀) above a critical coupling, **ξ_c (ℓ = 0) = 2.40 (a/m = 3), 0.79 (a/m = 6), 0.51 (Ellis)**; e.g. γ_t = 0.215/m (a/m = 3) and 0.399/m (a/m = 6) at ξ = 12, and none at ξ ≤ 1 for a/m = 3 (including conformal ξ = 1/6). Potential validated against the symbolic form to 10⁻¹².
- **The decisive subtlety: which time.** A global unstable mode grows in Killing time t, and the packet's Killing-time residence in the region cannot drop below the light-crossing time: **16.9 m (a/m = 3), 17.1 m (a/m = 6)** for ε = 0.1, saturating already at E ≈ 5 (residence 22.3 → 17.0 → 16.9 m for E = 1, 5, 20 at a/m = 3). So A = γ·t_res has the **floor A_floor = γ t_light = 3.62 (a/m = 3) and 6.83 (a/m = 6)** against A_max = 6.9: fast transit reduces proper time ∝ 1/E but not the exposure of a Killing-time mode. For a/m = 3 (ξ = 12) A stays below A_max for all E (4.8 at E = 1, 3.6 at E = 20); for a/m = 6 it passes only for E ≳ 5 and then by 1%. Where γ ≥ A_max/t_light ≈ 0.4/m even a light-speed transit does not help.
- For the black-bounce branch (a < 2m) t is not a global coordinate across the horizons, so a Killing-time mode criterion does not apply; a different (retarded-time or interior-evolution) analysis would be needed, which was not done.

## Part E. Verdict

| question | result |
|---|---|
| is the proper-time transit finite and horizon-regular? | **yes**: ≈ a-scale, from 0.4 m (a/m = 0.5) to 22 m (a/m = 6), ∝ 1/E at large E |
| is a finite parameter region consistent with a perturbation tolerance and causal passage? | **only in the toy**, and only if the instability is a Killing-time mode with γ < A_max/t_light; a/m = 3 (ξ = 12) passes at every E, a/m = 6 only marginally for E ≳ 5 |
| does time dilation weaken local instability? | **no**: proper time shortens with E, but a global mode's exposure has a light-crossing floor; the redshift/freezing is a feature of received signals |
| does it remove the NEC violation or supply the source? | **no** (T_kk = −2a²/(κR⁴) is observer-independent) |
| does it help the ghost source? | **not derived**; a ghost with an unbounded spectrum is not rescued by a finite transit |
| local finite-transit stability vs global stability, ghost-freedom, formation | **different questions**: only the first is addressed; the other three are unchanged (3D.3: no admissible source) |

## Limits

Prescribed benchmark; radial transit maps for three energies (L maps for a few values); the instability model is a test-field toy with a hand-chosen coupling (ξ up to 12 is large and not motivated by a theory); ε = 0.1 region definition (other ε give the same ordering in the table); interior (a < 2m) growth not analysed; no derived growth rate for the forced ghost source; no coupled evolution.

## Where this leaves the programme

Phase 3D.4 clarifies that a short proper time is not the same as short exposure to an instability: the right comparison is the instability's own time variable (Killing time for a static global mode, proper time only for a locally generated one). It does not change the hierarchy of obstacles: the source remains the problem. The open, bounded next step stays the one from 3D.3: one concrete derivable model (limiting-curvature/mimetic, or a Horndeski bounce) with its own derived perturbation spectrum, to which the amplification integral A = ∫γ dτ can then be applied honestly.
