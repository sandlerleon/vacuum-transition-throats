# Phase 3C.5: covariant spin transport and constraint closure (branch `phase3-recycling`)

Status: **stopped before axisymmetric evolution, as instructed. Verdict: the anti-aligned cancellation s = −2ω̃/κ is NOT dynamically maintained, and the homogeneous polarized system does not close. The mechanism is excluded in the reduced (homogeneous) class and is not supported as a derived result in any wider class.** No phenomenological exchange term was introduced. 30 automated checks pass. Phase 1–2 manuscripts and public releases are untouched.

Reproduce (`code/`): `phase3c5_closure.py` → `phase3c5_closure.json`, `tests/test_all.py`. Conventions as in Pasmatsiou–Tsagas–Barrow (PTB): signature (−,+,+,+), κ = 8πG, s² = s_ab s^ab/2, Frenkel condition s_ab u^b = 0, ω̃ = Riemannian vorticity.

## Priority check (the falsification)

On the cancellation surface the derived transports give

**d/dτ (s + 2ω̃/κ) = −α̇ s**  (α̇ = axial expansion rate; s ∝ n and ω̃ ∝ a_⊥⁻²).

It is nonzero unless the axis does not evolve. d ln(s/ω̃)/dτ = −α̇ exactly, so s/ω̃ ∝ 1/a_∥ and a collapse that shrinks the axis drives the system **through** the surface, not along it. The surface is invariant only for α̇ = 0, a cylindrical collapse whose shear fraction σ²/(Θ²/3) is fixed at **1/4**, the opposite of isotropization. **Anti-alignment is not preserved; this falsifies the mechanism without a new solver.**

## Steps 1–2. Action, field equations (what was used, what was verified)

The Einstein–Cartan action (1/2κ)∫R(Γ) + S_Weyssenhoff gives (i) the Cartan equation (torsion algebraically sourced by the spin density), (ii) the Einstein equation with the canonical energy-momentum tensor, (iii) angular-momentum conservation (the antisymmetric part of the canonical tensor equals the divergence of the spin density). For a Weyssenhoff fluid with Frenkel s_ab u^b = 0 and s_abc = s_ab u_c, the torsion vector vanishes (s_a = s_ba u^b = 0), so number conservation keeps its Riemannian form, and the system is equivalent to GR with the effective stress of PTB Eq. 77 plus a spin transport law. I used these standard results (Hehl, Obukhov–Korotky, PTB), checking the algebra that I could check symbolically; I did not rederive the action variation from scratch.

## Step 3. Spin transport (derived)

- **Compatibility with the Frenkel condition (A).** Aligned polarized spin s^{ab} = s ε^{abc}e_c is orthogonal to u for every s (residual 0, s_ab s^ab/2 = s²). The surface s = −2ω̃/κ is therefore algebraically allowed, but it is a **relation between two evolving fields**, not a free parameter: the spin density is s = nσ₀ with fixed spin per particle (see below), so it ties ω̃ to n.
- **Transport law (B).** Angular-momentum balance with Frenkel gives Ṡ^{ab} = −ΘS^{ab} + λ(u^b S^{ac}A_c − u^a S^{bc}A_c). Requiring the Frenkel condition to be preserved, d/dτ(S^{ab}u_b) = Ṡ^{ab}u_b + S^{ab}A_b gives the residual (1 − λ)(A_y s, −A_x s) in the (1,2) plane and **forces λ = 1** (Fermi–Walker transport). Consequence: **d ln s/dτ = −Θ exactly, for any acceleration A**: s = nσ₀. The direction is gyroscopically fixed (no spin–orbit exchange appears because the momentum correction is p^a = s^{ab}A_b, which vanishes for A = 0, and it cannot change |s| in any case).
- **Vorticity.** The Riemannian circulation law for geodesic flow, ω̃̇ = −(2/3)Θω̃ + σ_zz ω̃ = −2β̇ ω̃ (ω̃ ∝ a_⊥⁻²), is a kinematic (Helmholtz) identity independent of the field equations.

## Step 4. The β̇ mismatch of 3C.4 is a genuine closure failure

Energy balance of the cross term in the homogeneous class requires d ln(sω̃)/dt = −(α̇ + 3β̇) (re-derived from ρ̇ + Θ(ρ + p̄) + σ^{ab}π_ab = 0 with ρ_c = 2sω̃, p_⊥ = sω̃, π_zz = −(2/3)sω̃). The derived transports give −(α̇ + 4β̇). The three laws (s ∝ n, Helmholtz, energy balance) for the two rates (α̇, β̇) give the residual **β̇ = 0**: **the only closed solution is static**. It is not a wrong scaling assumption and not a missing term I can choose: both transports are theorems, and the energy balance is the Bianchi identity.

## Step 5. Constraints and sign checks (honest scope)

- The garbled extraction of PTB Eqs. (65)–(67) does not allow a reliable index-by-index verification of the torsionful constraints, so I did **not** claim one. Instead the Riemannian-effective constraints were derived/checked independently: on an explicit metric (flat space, rigidly rotating observers, Θ = σ = Ricci = 0, A ≠ 0, ω̃ ≠ 0) the vector constraint is **curl ω̃ + 2ε A ω̃ = 0 with k = 2** (all components), fixing the sign structure that PTB Eq. 66 must reduce to.
- **Commutator obstruction (new, decisive for step 7).** On the same metric, **D_[a D_b] f = −ω_ab ḟ** for generic f (ratio −1.000 in all components). With ω̃ ≠ 0 and ṅ = −Θn ≠ 0 this forces D_a n ≠ 0. **A rotating, time-dependent configuration cannot be uniform in the rest space.** Rotating collapse requires gradients (pressure forces, energy flux q_a, divergence terms) of the same order as the terms in the energy balance, so the 3C.4 reduced homogeneous class (and the rotation proxy of 3C.2) is not an exact class for polarized rotating collapse; Gödel-like cases are consistent only because ṅ = 0. Step 7 (a consistent rotating homogeneous background with torsion) is therefore **not attainable in the uniform class**; it would need a tilted/inhomogeneous construction, not a homogeneous one, and that is no longer the cheap controlled model.

## Step 6. Is the cancellation preserved, fine-tuned or forbidden?

| level | result |
|---|---|
| algebra (Frenkel, spin-tensor) | allowed |
| derived transports | not preserved: d(s + 2ω̃/κ)/dτ = −α̇ s |
| homogeneous class | forbidden for any moving configuration (β̇ = 0 forced) |
| enthalpy-corrected Kelvin form | d/dτ(s + 2ω̃/κ) = s(ḣ/h − α̇); the invariant locus α̇/β̇ = 4Bn/(1 − 4Bn) is a fine-tuned trajectory, and with ḣ ≠ 0 it is itself excluded by the commutator obstruction; formal only |
| numerical illustration (EC Bianchi I bounce, shear 0.05, B = 0.01 / 0.03 / 0.1) | s/ω̃ at the bounce is 15.9 / 5.8 / 2.8 times its initial value; the surface is crossed |

## Step 8. Cancellation of the source term is not cancellation of the shear

Even where the anisotropic source ω̃(ω̃ + κs/2) vanishes, shear is not removed: it only stops being generated. In GR dust with initial shear 0.05, σ²V² stays constant to 3.5×10⁻¹² through the whole singular collapse, so σ²/ρ grows as V⁻¹. Shear is fed by the Weyl term and the initial anisotropy; removal requires a V⁻² term exceeding it (the spin term of 3C.2). In the EC runs σ²V² is also conserved (ratio 1.00) through the bounce window. Non-singular transition and instabilities are unchanged from 3C/3C.4 (effective-fluid pathologies stand).

## Gate decision for 3C.5

| gate | criterion | result |
|---|---|---|
| Q1 | is s = −2ω̃/κ algebraically allowed (spin tensor, Frenkel)? | **yes** |
| Q2 | does the derived spin transport preserve it? | **no**: d(s + 2ω̃/κ)/dτ = −α̇ s |
| Q3 | does the homogeneous polarized system close? | **no**: β̇ = 0 forced (static only) |
| Q4 | does a rotating homogeneous background with torsion exist? | **not in rest-space-uniform form** ([D_a,D_b]f = −ω_ab ḟ) |
| Q5 | does cancelling the source cancel the shear? | **no** |
| **Verdict** | mathematically consistent? physically plausible? | **excluded in the reduced class; not demonstrated beyond it** |

## What this does and does not exclude

It excludes the specific mechanism "polarized anti-aligned Weyssenhoff spin cancels rotation-induced anisotropy in homogeneous collapse". It does **not** exclude a spin–orbit-coupled inhomogeneous rotating collapse in Einstein–Cartan gravity (which would need a real axisymmetric solver) or other torsion theories. Per the instructions I did not add exchange terms and did not start an axisymmetric evolution.

## Recommended next step

The reduced Weyssenhoff route is closed. Of the stated alternatives: **Einstein–Cartan–Dirac** (fermion spin density from a spinor field with a Nieh–Yan/Holst extension, where the axial current gives a four-fermion interaction and, unlike the Weyssenhoff fluid, there is no Frenkel-imposed fixed spin per particle) is the next well-posed theory test; it would need its own homogeneous baseline (3C-style) before any rotation. A tilted rotating Bianchi (VIII/IX/VII_h) baseline with ordinary matter is the prerequisite for any rotating claim.

## Limits

Action variation not rederived (standard results used); PTB Eqs. 65–67 not index-verified (extraction garbled; checks done on the Riemannian-effective constraints); the energy-balance step assumes A = q = D·π = 0 (homogeneous class); Kelvin-enthalpy variant is formal; Weyssenhoff fluid is a classical model of limited field-theoretic standing.
