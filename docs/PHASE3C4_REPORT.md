# Phase 3C.4: polarized spin and shear in Einstein–Cartan collapse (branch `phase3-recycling`)

Status: **stopped after 3C.4 as instructed. Decision: FAIL for starting axisymmetric numerical evolution, with a precise, constructive reason: the polarized-spin cross terms have exactly the structure that could help, but the reduced homogeneous system that contains them is not mathematically closed, and the missing piece (the spin–orbit transport law) is not supplied by the paper.** The best-case benefit of the mechanism is bounded: it can at most remove the penalty that rotation adds, not create a bounce. 29 automated checks pass. Nothing here establishes a daughter universe, recycling or a wormhole.

Reproduce (`code/`): `phase3c4_polarized.py`, `tests/test_all.py`. Units κ = 8πG as in the paper unless stated.

## Step 1–2. The complete paper, and what it actually contains

The full text (24 pages) of Pasmatsiou, Tsagas & Barrow, *Kinematics of Einstein–Cartan universes* (Cambridge repository bitstream, extracted locally) was read. Corrections to my earlier extract:
- **Eq. (77)** is the *effective stress tensor* of a Weyssenhoff fluid in Riemannian form, T_ab = (ρ − κs²/4)u_a u_b + (p − κs²/4)h_ab − ½h^{cd}∇̃_d(s_ca u_b + s_cb u_a), and **Eq. (78)** is the shear-free, irrotational GR-form Raychaudhuri equation with spin. They are not shear or vorticity propagation equations.
- The **shear and vorticity equations with torsion are the general Eqs. (63) and (64)** (in terms of the torsionful vorticity ω_ab = ω̃_ab − κs_ab/2). **No Riemannian-form shear or vorticity propagation equation for the Weyssenhoff fluid is given**, no Hamiltonian constraint with torsion is given, and no spin transport law is given.
- The authors state that the **Weyssenhoff fluid is incompatible with the Cosmological Principle** (the spacelike spin tensor defines a preferred direction) and "of limited use from the field theoretical perspective", and that spin aligned along a shear eigenvector, and the combined spin–vorticity action, are left as future work (their §9).
- Their stated conclusion for collapse: worldlines do not focus if κs² > ρ + 3p (their Eq. 76 with Λ = 0), i.e. the quadratic spin term alone, with random orientation.

## Step 3. Derivation of the polarized-spin cross terms (new, from Eq. 77)

For spin s_ab, vorticity ω_ab and an axisymmetric shear σ_ab all about one axis, the part of Eq. (77) that exists only for polarized spin is evaluated symbolically (sympy). The sign convention of ω is **fixed by requiring the paper's own Raychaudhuri cross term, −κ s_ab ω̃^ab (Eq. 76)**; this fixes D_d u_b = (Θ/3)h_db + σ_db + ω_bd. Results (`phase3c4_polarized.json`):

- cross-term energy density T_uu^{(c)} = **2 s ω**; spatial stress **p_⊥ = s ω in the equatorial plane, p_∥ = 0**; trace T^{(c)} = 0; the shear–spin commutator [s, σ] vanishes for axisymmetric alignment;
- Raychaudhuri: 2ω̃² − 2κ s ω̃ + κ²s²/2 = ½(κs − 2ω̃)² (a perfect square, as in the paper); co-aligned spin reduces the repulsion, anti-aligned increases it (up to 2(ω̃ + κ|s|/2)²);
- **shear source** (standard 1+3 shear equation with the effective anisotropic stress): −ω̃_c⟨aω̃_b⟩^c + (κ/2)π_ab ∝ **ω̃(ω̃ + κs/2)** · diag(1/3, 1/3, −2/3). It **vanishes at s = −2ω̃/κ** (anti-aligned spin of magnitude 2ω̃/κ).

So, in this derivation, **anti-aligned polarized spin is the configuration that both increases the repulsive Raychaudhuri term and removes the rotation-induced flattening source; co-aligned spin does the opposite on both counts.** Both the Raychaudhuri square and the shear source have their zero at |s| = 2ω̃/κ, with opposite alignment.

## Steps 4–6. Aligned, anti-aligned, partially polarized; closure and conservation

| configuration | Raychaudhuri repulsion | rotation-induced shear source |
|---|---|---|
| random orientation | 2ω̃² + κ²⟨s²⟩/2 | ω̃² (unchanged) |
| co-aligned s > 0 | reduced; zero at s = 2ω̃/κ | increased by κsω̃/2 |
| anti-aligned s < 0 | enhanced | reduced; zero at s = −2ω̃/κ |
| partial polarization | interpolates (linear in s) | interpolates |

**The reduced system is not closed.** In the homogeneous axisymmetric (Bianchi I) class, the cross term (ρ_c = 2sω, p_⊥ = sω, p_∥ = 0) is conserved only if **d ln(sω)/dt = −(α̇ + 3β̇)**. With the standard separately conserved transports, s ∝ n (fixed spin per particle) and ω̃ ∝ a_⊥⁻² (circulation law from the vorticity equation, Phase 3C.2), one obtains **−(α̇ + 4β̇): the mismatch is exactly β̇**. So either the transports must exchange angular momentum between spin and orbit, or such a homogeneous aligned rotating configuration does not exist (consistent with the known absence of rotating Bianchi I). A reading of the paper's torsionful vorticity equation (Eq. 64, Frenkel condition, s ∝ n) indeed contains an exchange term, ω̃̇ = −2H_⊥ω̃ + c κΘ s with c = 5/6 in my transcription (index placement of the last term is a reading of the paper, not independently verified), which would conflict with the conservation requirement at generic points. The constraint equations (65)–(67) (not extracted) may restrict the allowed configurations further. I did not impose a permanent alignment or an ad hoc transport law to force consistency.

## Step 7. Comparison against the matched unpolarized baseline: a bound, not a simulation

Because the reduced system is not closed, a polarized run would not be a faithful solution. A **what-if bound** is possible: suppose a fraction f of the rotation-induced anisotropic stress is cancelled (by the cross term at |s| = 2ω̃/κ). From the Phase 3C.2 Bianchi-I model with rotation proxy v₀² = 0.2, the critical spin strength B_c is: f = 0: **0.0305**; f = 0.5: 0.0124; f = 0.9: 1.4×10⁻³; f = 0.99: 5.7×10⁻⁵. **The bound: even complete cancellation restores the unrotated baseline (B_c → ~10⁻⁸), it cannot do better.** Reasons: rotation without torsion never bounces (Phase 3C.2), and the cross energy 2sω scales as V^{−5/3} (s ∝ V⁻¹, ω̃ ∝ V^{−2/3}), subdominant to the V⁻² spin and shear terms that decide the bounce.

## Step 8. Perturbative instabilities

Unchanged in this analysis: the cross terms scale like V^{−5/3} and do not enter the isotropic effective fluid (ρ_eff = n − Bn², c_s² = −2Bn/(1 − 2Bn)); the c_s² < 0 / singular-at-n* / superluminal structure of Phase 3C and the extra k-growing instability stand. A polarized medium additionally breaks isotropy, which would add anisotropic perturbation structure not analysed.

## Step 9. Magnetic fields

Treated as secondary and not combined: no consistent polarized model exists to attach them to.

## Gate decision for 3C.4

| gate | criterion | result |
|---|---|---|
| P1 | does the polarized-spin structure point in the right direction? | **yes (structure)**: anti-aligned spin raises repulsion and lowers the shear source; zero of the shear source at s = −2ω̃/κ |
| P2 | is there a closed, consistent reduced system containing it? | **no**: conservation mismatch β̇; spin–orbit transport law missing |
| P3 | can it by itself create a bounce? | **no**: best case removes the rotation penalty (bound above) |
| P4 | does it remove the effective-fluid pathologies? | **no** |
| **Decision** | start axisymmetric numerical evolution? | **FAIL (do not start)** |

## What the cancellation would physically require

At the cancellation point κ|s| = 2ω̃ with s = σ₀n, so ω̃ = κσ₀n/2. At the spin-fluid's n* = 2/(κσ₀²) this gives ω̃ = 1/σ₀, i.e. in natural units a vorticity of order the Planck frequency for σ₀ = ħ/2. Before that epoch κσ₀n/(2ω̃) ∝ a⁻¹ is tiny for any macroscopic rotation: the spin cross term is utterly negligible until the very last stage of the collapse, where the needed orbital vorticity is itself Planckian. Whether the compression can amplify an initial rotation to that value (ω̃ ∝ a_⊥⁻², ratio of scales ~ (n*/n₀)^{2/3}) while keeping the shear it induces below the spin term is exactly the unresolved question.

## What would make this lead testable (specific, not "more of the same")

1. **The spin transport law.** Derive from the Einstein–Cartan field equations and the Frenkel condition the evolution of s_ab for a *rotating* Weyssenhoff fluid, including the torsion–vorticity exchange; check it against stress-energy conservation (the missing closure).
2. **The constraint equations (65)–(67)** of the paper for a rotating, aligned configuration, to learn whether a homogeneous rotating aligned class exists at all.
3. **A true rotating homogeneous class with torsion** (tilted Bianchi type VII/VIII/IX or Gödel-type), replacing the orbital-motion proxy.
4. Only then, the matched comparison against the Phase 3C.2/3C.3 baseline and the perturbative analysis.

## Limits

Cross terms derived from Eq. (77) as transcribed in the text (Greek symbols are partly lost in the PDF text extraction; the κ factors were matched by requiring the paper's own Eq. 76 cross term); homogeneous axisymmetric alignment only; Eq. (64) used as read; no spin transport law; no constraint equations; the Weyssenhoff fluid is a classical model of limited field-theoretic standing (the paper's own caveat).
