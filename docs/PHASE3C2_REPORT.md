# Phase 3C.2 and 3C.3: rotation, spin–torsion and magnetic fields — feasibility (branch `phase3-recycling`)

Status: **stopped before any axisymmetric evolution, as instructed. In the exact controlled model I could build, neither rotation nor a magnetic field helps the spin–torsion bounce: both raise the critical spin strength, and the magnetic field at the bounce is far too weak to repair the effective fluid's pathologies. Two genuine leads remain open. Nothing here establishes a daughter universe or recycling.** 28 automated checks pass. Phases 1–3C and all public releases are untouched.

Reproduce (`code/`): `phase3c2_kinematics.py`, `phase3c2_bianchi.py`, `phase3c3_mhd.py`, `phase3c3_bounce_fields.py`, `tests/test_all.py`. Units 8πG = κ = 1 unless stated.

## 0. Inputs checked against the literature (this session)

- **Pasmatsiou, Tsagas & Barrow, "Kinematics of Einstein–Cartan universes", arXiv:1611.07878** (read via the web; the last ~10% of the page was truncated). Verified: the Einstein–Cartan Raychaudhuri equation for a Weyssenhoff fluid (their Eq. 73) has the GR form; the Riemannian form (Eq. 76) is Θ̃′ = −Θ̃²/3 − (κ/2)(ρ + 3p) − 2(σ̃² − ω̃²) + … + (κ²/2)s² − κ s_ab ω̃^ab; the effective fluid is ρ_eff = ρ − κs²/4, p_eff = p − κs²/4 (s² = s_ab s^ab/2); the kinematic relation is ω_ab = ω̃_ab − κs_ab/2; the paper assumes random spin orientation (⟨s_ab⟩ = 0, ⟨s²⟩ ≠ 0), which removes the cross term. **This removes the main assumption flagged in Phase 3C: the effective-fluid coefficient is now confirmed, not just remembered.** Not found in the visible text: a Hamiltonian constraint with torsion, and the Riemannian-form shear equation with spin terms.
- **Shear-free conjecture** (a rotating, expanding barotropic perfect fluid cannot be shear-free): proven for dust (Ellis), incoherent radiation and several special cases, **open in general** (checked via search; sources: arXiv:1709.03547, arXiv:1510.05798, arXiv:gr-qc/9702035). It needs ρ + p ≠ 0, which fails at the spin fluid's n*.
- Pasmatsiou–Tsagas–Barrow also state that spin can "inhibit (perhaps even reverse)" collapse; Popławski's bounce work assumes closed, homogeneous, isotropic backgrounds, with the shear statement given only in slides (not verified as a derivation).

## 1. Derived and verified identities (`phase3c2_kinematics.json`)

1. **Effective-fluid link.** −(κ/2)(ρ_eff + 3p_eff) − [−(κ/2)(ρ + 3p) + κ²s²/2] = 0 (symbolic): the spin term of Eq. 76 is exactly the effective fluid's.
2. **Spin and vorticity enter as a perfect square.** For aligned spin and rotation, 2ω̃² − κ s_ab ω̃^ab + κ²s²/2 = 2(ω̃ − κs/2)² (symbolic zero difference), consistent with ω_ab = ω̃_ab − κs_ab/2. For random orientation the cross term vanishes: 2ω̃² + κ²⟨s²⟩/2 ≥ 0, both repulsive. Corotating, polarized spin with ω̃ = κs/2 cancels the repulsion in the Raychaudhuri equation: total spin/rotation effects are **not additive**.
3. **Circulation law** from the verified vorticity equation with axisymmetric shear: ω̇/ω = −2H_⊥, i.e. ω ∝ a_⊥⁻².
4. **Scaling of every term in the constraint** θ²/3 = κρ + σ² − ω² − R³/2 + Λ in an isotropic contraction (V → 0): dust V⁻¹, curvature V⁻²ᐟ³, vorticity and magnetic V⁻⁴ᐟ³, orbital kinetic V⁻⁵ᐟ³, **shear and spin both V⁻²**. Hence as V → 0 only spin (negative) and shear (positive) matter; a bounce requires the net V⁻² coefficient to be negative, B N₀² > Σ_asymptotic. Rotation and magnetic fields can matter only by *changing Σ*.

## 2. An exact controlled model: axisymmetric Bianchi I (`phase3c2_bianchi.py`)

Metric ds² = −dt² + e^{2α}dx² + e^{2β}(dy² + dz²). Einstein equations: 2α̇β̇ + β̇² = ρ; 2β̈ + 3β̇² = −p_x; α̈ + β̈ + α̇² + β̇² + α̇β̇ = −p_y. Matter, each separately conserved (checked): dust ρ = n = N₀/V; random-spin fluid ρ_s = p_s = −Bn²; **rotation proxy**: azimuthally averaged orbital motion of the dust in the equatorial plane with conserved specific angular momentum (v ∝ a_⊥⁻¹, from the circulation law), ρ_r = ½nv², p_y = p_z = ρ_r, p_x = 0; **magnetic field** along x with flux conservation, ρ_B = b₀²e^{−4β}/2, p_x = −ρ_B, p_y = p_z = +ρ_B.

**What the proxy is, and is not.** It reproduces the circulation law and the anisotropic centrifugal stress (flattening), and it is an exact solution class of the Einstein equations. It is **not** a rotating spacetime: there is no frame dragging and no global vorticity. A true rotating homogeneous model (tilted Bianchi, Gödel-type) could behave differently; the conclusions below are for this proxy.

**Checks:** constraint residual 3×10⁻⁹; conservation 4×10⁻⁶; GR limit (dust, shear) is Kasner-like with σ²V² constant to 10⁻¹⁰ and no bounce; the analytic bounce volume V_b = BN₀ − Σ matches the numerics (0.017535 both).

**Results.** Critical spin strength B_c (smallest B giving a bounce; initial state V = 1, n = 1, contracting; B_c is in units of N₀² and equals the initial spin-to-dust energy ratio):

| model (matched start) | B_c | at 10 B_c: V_bounce, max σ² |
|---|---|---|
| A: nonrotating EC, isotropic | < 10⁻⁸ (any B > 0 works) | 10⁻³, 2×10⁻⁹ |
| A′: nonrotating EC, initial shear 0.05 | 0.0025 (= Σ₀) | 0.023, 4.8 |
| B: rotating or magnetized GR, **no torsion** | **no bounce in any case** (5 of 5 singular) | – |
| C: rotating EC, v₀² = 0.01 / 0.05 / 0.2 | 5.5×10⁻⁴ / 5.0×10⁻³ / 0.031 | 0.0046 / 0.041 / 0.24 |
| D: magnetized, ρ_B0 = 0.01 / 0.05 / 0.2 | 5.3×10⁻⁴ / 7.0×10⁻³ / 0.053 | 0.0046 / 0.060 / 0.42 |
| E: rotating + magnetized, 0.05/0.05 and 0.2/0.2 | 0.016 and 0.102 | 0.13 and 0.76 |

Fits: B_c ≈ 0.27 (v₀²)^1.34 (rotation) and B_c ≈ 0.67 (ρ_B0)^1.54 (magnetic). **Both rotation and the magnetic field act as sources of shear that the spin term must overcome; the critical spin strength rises with each and more with both.** Without torsion nothing prevents the singularity. This answers the instructions' main question for this model: the rotation proxy's centrifugal stress is anisotropic and generates shear that *worsens* the conditions for a spin-supported bounce. Whether true vorticity (which the proxy does not contain) would act differently is untested, and the +2ω² term of the Raychaudhuri equation is not represented by it.

## 3. Magnetized spin fluid: local wave analysis (`phase3c3_mhd.py`)

Relativistic ideal MHD, locally flat, background field b, effective EOS ρ(n) = n − Bn². Derived with sympy and factorized: **(ω² − k²v_A²cos²ψ) · [ω⁴ − ω²k²(v_A² + c_s²(1 − v_A²sin²ψ)) + c_s²v_A²k⁴cos²ψ] = 0**, v_A² = b²/(w + b²) (the factorization is verified symbolically). Consequences, with c_s² = −2Bn/(1 − 2Bn) < 0 below n*:

- the product of the magnetosonic roots is c_s²v_A²cos²ψ < 0, so **one root has ω² < 0 for every angle except exactly perpendicular propagation: a magnetic field cannot cure the imaginary sound speed** (79 of 90 scanned cases unstable);
- exactly perpendicular propagation is stable only if v_A² > −c_s²/(1 − c_s²);
- the inertia becomes w + b² instead of w, and the perpendicular fast speed² at n* (w = 0) is (b² − 2Bn²)/b²: finite, real and subluminal iff b² > 2Bn². So a strong enough field **regularizes the singular fluid momentum equation at n*, but only for perpendicular modes.**

**Does the field present at the bounce reach that threshold?** No (`phase3c3_bounce_fields.json`). In the magnetized bounce solutions, at 1.5–10 B_c, the bounce sits at n/n* = 2–7, ρ_B/|p_spin| = 2×10⁻⁴ to 0.15 (threshold: ≥ 1) and w + b² < 0 (v_A² < 0): the fluid inertia is still negative. The field would have to exceed the spin pressure, which in this model drives B_c up by the scaling above.

## 4. Comparison of the three models and decision gates

| | A: nonrotating EC | B: rotating GR, no torsion | C: rotating EC (proxy) |
|---|---|---|---|
| regular transition in the Bianchi-I class | yes if Σ < BN₀² | **no** | yes, but B_c rises as (v₀²)^1.34 (5.5×10⁻⁴ → 0.031 over v₀² = 0.01–0.2) |
| effect of rotation/field on the condition | – | none | worsens |
| effective-fluid pathologies (c_s² < 0, singular at n*, c_s² > 1) | present | not applicable | present (unchanged) |

| gate | question | result |
|---|---|---|
| G6 | does rotation reduce the critical spin strength? | **no** (raises it, exponent 1.34) |
| G7 | does a magnetic field reduce it? | **no** (raises it, exponent 1.54) |
| G8 | does a field regularize the effective fluid at the bounce? | only perpendicular modes, and the field present is 10⁻⁴–10⁻¹ of the threshold |
| G9 | does any improvement survive shear/perturbations? | none found |
| G10 | is axisymmetric evolution justified? | **not for the barotropic spin fluid with random spin orientation** |

## 5. Leads that remain open (not tested)

1. **Polarized spin and the perfect-square structure.** For aligned spin the Raychaudhuri terms combine into 2(ω̃ − κs/2)², so corotating spin could cancel the vorticity effects. Whether it also cancels the vorticity-induced *shear* depends on the Riemannian-form shear equation with spin terms (the truncated Eqs. 77–78 of the paper), which I could not read. This is the most promising untested direction: it would address exactly the mechanism that raises B_c here.
2. **A true rotating homogeneous model** (tilted Bianchi/Gödel-type) to test the rotation conclusion beyond the orbital-motion proxy, including frame dragging and the non-integrable rest space.
3. **Non-barotropic spin dynamics** (evolving spin density), which would change the sound-speed structure that fails in Phase 3C.
4. **Perpendicular-field geometries** in which only perpendicular modes matter, if a configuration with ρ_B ≳ |p_spin| at the bounce can be reconciled with the shear it induces.

## 6. Limits

Bianchi I has no spatial curvature (irrelevant near the bounce but not at the start); the rotation is a proxy; spin is random-oriented and barotropic; the local MHD analysis is a WKB statement in a flat patch; no inhomogeneous or axisymmetric evolution was run; the shear-free statement is a conjecture in general; Popławski's anisotropy claim was seen only in slides.
