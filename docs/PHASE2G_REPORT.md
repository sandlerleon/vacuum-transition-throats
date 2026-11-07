# Phase 2G-A: perturbation spectrum on the actual cap, with consistent shell boundary conditions (branch `phase2-bridge`)

Status: **the finite cap with consistent no-flux boundary conditions does not remove the inhomogeneous instability, and no choice of local boundary condition at the shell can.** Phase 2G-B (a scalar-dependent shell tension) was analysed, not implemented, because it is blocked by a field-continuity inconsistency of the construction that Phases 2D–2F did not check. 26 automated checks pass. The CQG manuscript is untouched; the Phase 2E manuscript carries a short added subsection.

Reproduce (`code/`): `phase2g_comoving_nec.py`, `phase2g_cap_spectrum.py`, `phase2g_nu_curve.py`, `phase2g_nu_low.py <nu>`, `phase2g_shell_data.py`, `tests/test_all.py`.

## 1. The boundary condition that Phases 2D–2F left undetermined

The scalar field has to obey some condition at the shell. For a homogeneous interior and a shell moving with rapidity β relative to the fluid, ∂_nφ = sinh β · φ̇ and the dust crosses the shell. A no-flux (Neumann) condition ∂_nφ = 0 with zero normal dust velocity is therefore satisfied by the background **only if β = 0, i.e. for a comoving shell**. The moving periodic shells of Phase 2D need an absorbing/emitting condition (the flux law, Phase 2C) that no field theory supplied; they are not compatible with a field confined to the cap. The consistent passive model is the comoving shell with a field confined to the cap (Neumann), for which the flux law gives zero exchange (this removes the "97% energy exchange" problem of Phase 2E §4 for that model).

**Do comoving shells with σ > 0 and σ + p ≥ 0 exist around the stable cycle?** Scan over exterior vacuum density −qV_f, q ∈ {1, 0, −0.5, −0.9}, χ₀ ∈ (0.505π, 0.995π) and the exterior mass (9900 combinations, `phase2g_comoving_nec.json`): 3486 satisfy both conditions over the whole cycle (this is far more than for the larger-excursion orbits of Phase 2 step 1, because the stable cycle changes a by only 1%). But for the exterior that matches the model, q = 1 (V(±v) = −V_f, AdS), they exist **only for χ₀ ≲ 0.6π**, a cap not much larger than a hemisphere, with large tension σ ≈ 2.3–2.8 and a strongly time-dependent p/σ in [−0.99, −0.54] (`phase2g_shell_data.json`). The shell matter would need an equation of state that changes by a factor of two within a cycle.

## 2. The cap spectrum

Radial modes regular at the pole, Q_ν = sin(νχ)/sin χ, k² = ν² − 1. The no-flux conditions (∂_χδφ = 0 and dust normal velocity U = vQ′ = 0 at χ₀) both reduce to Q_ν′(χ₀) = 0, i.e. **tan(νχ₀) = ν tan χ₀**. (ν = 1 is the homogeneous mode.) At χ₀ = π the roots are the integers (the three-sphere of Phase 2E). The Floquet multipliers use the Phase 2F equations with non-integer k².

| χ₀/π | ν₂, ν₃, ν₄ (lowest) | G of the lowest modes | worst G (ν ≤ 40) |
|---|---|---|---|
| 0.55 | 2.76, 4.57, 6.38 | 10.4, 64, 210 | 7.5×10⁴ |
| 0.60 | 2.57, 4.21, 5.86 | 7.8, 48, 155 | 7.6×10⁴ |
| 0.70 | 2.28, 3.66, 5.07 | 4.6, 30, 93 | 7.6×10⁴ |
| 0.80 | 2.11, 3.28, 4.49 | 2.9, 20, 61 | 7.6×10⁴ |
| 0.90 | 2.02, 3.06, 4.12 | 2.1, 15, 45 | 7.7×10⁴ |
| 0.95 | 2.00, 3.01, 4.02 | 1.9, 14.5, 41 | 7.7×10⁴ |

The shell-consistent caps (χ₀ ≲ 0.6π) have lowest-mode growth 7.8–10.4 per cycle. The lowest mode of a larger cap is closer to the gauge-degenerate point ν = 2 (k² = 3), but the growth is a smooth function of ν there (`phase2g_nu_curve.json`): G(2.002) = 1.90, G(2.15) = 3.3, G(3) = 14.34 (the S³ value), G(5) = 88.6.

**General statement.** The growth G(ν) (the same function for any boundary condition, because the instability is a bulk property) satisfies G(ν) > 1 for every ν > 1 and tends to 1 only as ν → 1 (G = 1.08, 1.27, 1.44 at ν = 1.02, 1.2, 1.7; `phase2g_nu_low.json`), then rises monotonically for ν > 2 up to about 10⁵. A boundary condition at the shell, whatever it is, only selects which ν values occur. Any local condition leaves infinitely many modes with ν spaced by about π/χ₀, so there are always modes with ν ≳ 3, where G > 14. Hence **no boundary condition can make the cap stable**; the neutral point is the homogeneous mode alone. Numerical statement, not a proof.

## 3. Phase 2G-B: the scalar-dependent shell tension

For S_shell = −∫d³ξ √−h σ(φ) with a continuous bulk field φ, variation gives the scalar junction condition [∂_nφ] = σ′(φ_s), a shell stress S_ab = −σ(φ_s)h_ab, i.e. **p = −σ exactly (w = −1)**, and a linear Robin condition for the perturbations. Two things stop it from being a "next experiment":

1. **Field continuity across the shell.** The construction keeps a static exterior with constant field. For a vacuum exterior of density −qV_f the field would sit at its minimum, φ = v. But a continuous field requires φ_in(shell) = φ_out, while the cycle has φ(t) oscillating between 0.223 v and 0.612 v (`phase2g_shell_data.json`). A thin transition between φ_in(t) and v is a domain wall whose tension is fixed by the potential: the thin-wall estimate ∫_{φ_in}^{v}√(2(V − V_vac))dφ is 0.38 down to 0.11, against the σ ≈ 2.3–2.8 that the NEC-satisfying comoving shells need: a mismatch of a factor 6 to 20. Alternatively the field exists in the exterior with φ(t) matched at the shell, which makes the exterior non-static and not Schwarzschild. Phases 2D–2F silently assumed the field lives only in the cap; the junction then has no field action to derive a boundary term from.
2. **The stress history.** A scalar-tension shell has p/σ = −1 at all times, whereas the consistent shells of §1 need p/σ between −0.99 and −0.54, and the required history cannot be reproduced.

The Robin condition would at most move the roots ν_j (§2); it cannot lower G(ν) for ν ≳ 3. So Phase 2G-B, as posed, cannot change the conclusion even if the first obstruction were overcome. Phase 2G-C (nonlinear evolution) is not warranted.

## 4. Corrections and consequences for earlier phases

- The Phase 2D/2E moving shell (w = −1/2) is an effective junction model with an unexplained flux. The comoving NEC shells of §1 are self-consistent with a confined field, but only for χ₀ ≲ 0.6π with q = 1 and with a non-constant p/σ.
- Phase 2E's statement that the field is a bulk field with a static vacuum exterior is inconsistent for the oscillating field; the consistent reading is a field confined to the cap.
- The finite-cap analysis does not change the headline number: the worst-case growth is 7.7×10⁴ per cycle and the lowest mode grows by at least 1.9 (7.8–10.4 for shell-consistent caps).

## 5. Limits

Linear; radial (ℓ = 0) modes only; the stable cycle at C/C_c = 0.35, μ = 20, ρ_crit = 0.2 V_f only; modes ν ≤ 40; the gauge-degenerate neighbourhood of ν = 2 is smooth in the numerics but the n = 2 gauge mode itself is not treated; the boundary condition is derived from the no-flux requirement, not from a shell action; the junction perturbation conditions for δσ, δp are not imposed (they determine the shell's matter response and do not constrain the interior modes once the shell stays comoving).
