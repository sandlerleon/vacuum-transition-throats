# Phase 3D.0: Einstein–Cartan–Dirac feasibility (branch `phase3-recycling`)

Status: **stopped after 3D.0 as instructed. Verdict: NO-GO for a stable, mean-field (barotropic-fluid) ECD bounce; the microscopic derivation is clean, but it reproduces the Phase 3C structure exactly and a general theorem shows the pathology is unavoidable for any barotropic fluid. No daughter-universe claim is made or supported.** 31 automated checks pass.

Reproduce (`code/`): `phase3d0_ecd.py` → `phase3d0_ecd.json`, `tests/test_all.py`. Units ħ = c = 1, κ = 8πG; numbers quoted in Planck units (G = 1).

## Step 1. Fundamental equations (derived and checked)

Action S = ∫e [R(ω,e)/2κ + L_D(e,ω,ψ)], independent spin connection ω = ω̃ + K, minimal coupling.

| item | result | how checked |
|---|---|---|
| γ^{abc} | = i ε^{abcd} γ_d γ⁵ for all 24 index triples | explicit Dirac matrices, signature (+,−,−,−) |
| coupling of the contortion | **only the totally antisymmetric part couples**, L_K = (i/4)K_{abc}ψ̄γ^{abc}ψ = c_K K_{abc}ε^{abcd}A_d, **\|c_K\| = 1/4**; vector and tensor torsion decouple (residual 2×10⁻¹⁶ for a random contortion) | numerics |
| gravity quadratic term | K_{abc} = ε_{abcd}a^d ⇒ R(Γ) = R̃ + **6 a·a**; pairing K_{abc}ε^{abcd}A_d = −6 a·A | sympy, signature (−,+,+,+) |
| Cartan equation + elimination | a^d = −κA^d/4·(±), **L_eff = −(3κ/16) A_a A^a = +λ[(A⁰)² − \|A\|²], λ = 3κ/16 = 3πG/2** | sympy; agrees with the literature value (Freidel–Minic–Takeuchi, Perez–Rovelli) in the Einstein–Cartan limit |

The coefficient is negative-definite whatever the sign of c_K. Torsion is algebraic: no new propagating degrees of freedom, hence no ghost sector from torsion itself. The interaction is non-renormalizable (dimension 6).

**Weyssenhoff is not the Dirac spin tensor.** The Dirac canonical spin tensor is totally antisymmetric, s_{abc} = ε_{abcd}S^d; the Weyssenhoff form s_{abc} = s_{ab}u_c is its u-projection only. The energy of the contact term differs accordingly (below), so the Phase 3C/3C.5 spin-fluid model is not simply the fermion theory.

## Step 2. Matter states (what is a state, what is an assumption)

A **classical c-number Dirac spinor** violates Pauli exclusion and cannot be isotropic alone (its spatial vectors A^i, V^i cannot both vanish); a **fermion ensemble** needs degeneracy pressure; the **quantum expectation** of the contact term has Hartree and Fock (exchange) parts. I took the physically defensible case: a degenerate gas, first-order Hartree–Fock, local density matrix D, with ε_int = −λ Σ_μ η_μμ [Tr(M_μD)² − Tr(M_μDM_μD)], M_μ = γ⁰γ^μγ⁵.

| state | ε_int / (λ n²) | meaning |
|---|---|---|
| unpolarized, massless (g = 2) | **−1/2** | B = λ/2 |
| single-helicity Weyl gas (g = 1) | −2 | B = 2λ (Hartree −1, spatial exchange adds −1) |
| unpolarized, nonrelativistic | −3/2 | B = 3λ/2 |
| partially polarized NR (fraction P) | −(3/2)(1 − P²) | polarization **removes** the effect |
| fully polarized NR | **0** | Pauli: one internal state, contact term cancels exactly |

Every state has the same form ε_int = p_int = −B n² (Weyssenhoff form, but with its own B and no spin-fluid assumption). Two consequences: (i) **polarization does not help**: the Phase 3C.4 idea of using polarized spin has no Dirac counterpart (a polarized gas feels nothing); (ii) the NR coefficient is 4.5× the Weyssenhoff value κ/16 = λ/3. For realistic masses m ≪ M_Pl the gas is ultra-relativistic at the relevant densities (k_F ~ M_Pl ≫ m), so the massless rows apply. Anisotropic stresses vanish for unpolarized/isotropic D; the stress tensor is a perfect fluid with ρ = ε_kin − Bn², p = p_kin − Bn² (admissible: no negative-definite pieces at low density). Caveats: antiparticles/Dirac sea ignored; chirality is conserved only for m → 0; HF is first order.

## Step 3. Homogeneous baseline (flat FLRW, exact)

ρ(n) = C n^{4/3} − B n², p = nρ′ − ρ. For the unpolarized massless gas (C = 2.32, B = 3π/4): bounce at **n_b = 0.977** (Planck units, k_F = 3.07 M_Pl), analytic and numerical agree, Friedmann constraint residual 6×10⁻¹²; p_b = −1.5, **R̃_b = 113 (finite)**. For the Weyl gas n_b = 0.173, k_F = 2.17, R̃_b = 14. The geometry is regular (finite curvature and torsion), the bounce is genuine within this model (ρ + 3p < 0 at ρ = 0). Characteristic densities (any C, B): chemical potential μ = dρ/dn vanishes at **n* = (2/3)^{3/2} n_b = 0.544 n_b** (ρ + p = 0: the stress is vacuum-like and the fluid velocity is undefined, an effective-fluid artifact), and **dp/dn < 0 for n > (2/9)^{3/2} n_b = 0.105 n_b**. The spinodal is reached at ε_int/ε_kin = 2/9 (so first-order HF is still semi-quantitative there), at t = 1.12 vs bounce at t = 1.73 (from n = 0.02 n_b); the bounce itself sits at ε_int/ε_kin = 1 (first-order HF unreliable) and at k_F ≈ 2–3 M_Pl (beyond any EFT cutoff).

## Step 4. Stability: a theorem, not a case study

**Theorem (barotropic bounce).** Let F(n) = ρ(n) − c n^{2/3} (c = 0 flat, c > 0 closed), H² = κF/3, a bounce at F(n_b) = 0 with ä > 0 ⇔ F′(n_b) < 0. If the chemical potential is positive at some lower density (F′ > 0, as for any ordinary gas), then (i) by the intermediate value theorem there is n* < n_b with F′(n*) = 0 (for c = 0: ρ + p = nμ = 0), and (ii) by the mean value theorem F″ < 0 somewhere below n_b, i.e. **dp/dn = nρ″ < −(2c/9)n^{−1/3} ≤ 0: negative compressibility (a thermodynamic spinodal) precedes every such bounce.** Verified on 4000 random EOS: all 312 bouncing cases show both. This generalizes the Phase 3C c_s² result to any microphysics expressible as a barotropic fluid, including ECD gases.

What the theorem does **not** say: the dynamic growth rate. I did not derive the collective-mode (RPA/Landau) dispersion, so the effective c_s² is **not** interpreted as a propagation speed; the static (k → 0) instability of the uniform state at dp/dn < 0 is rigorous Fermi-liquid thermodynamics (Pomeranchuk-type), but the growth rate relative to the Hubble rate near the bounce is **not computed**. Ghosts/hyperbolicity: torsion is algebraic, the Dirac equation with a local nonlinearity stays first-order symmetric hyperbolic, so no new ghost or hyperbolicity failure follows from the structure (argued, not proven for the full coupled system). Sensitivity: the 2/9 and (2/3)^{3/2} numbers are independent of C and B, so the no-go is not a tuning issue.

## Step 5. Relevance to recycling

- **Embedding.** The Misner–Sharp mass of a homogeneous ball is m = (κ/6)ρR³ (flat) = 0 at the bounce, or ∝ ρ/n (closed, equal to the energy per particle, which the interaction changes by O(1)). For a closed ball 2m/R = sin²χ₀(1 + ȧ²) ≤ 1 at ȧ = 0: **the boundary is never trapped at the bounce**; a bounce after trapped-surface formation requires the Misner–Sharp mass to drop, i.e. energy shed outward across the boundary, impossible into a vacuum exterior (this is the Phase 3C G3 result, now shown to be independent of the EOS).
- So any bounce inside a horizon must be a **core** bounce with a collapsing, non-bouncing envelope: an inhomogeneous spherical (LTB-like) problem. With a barotropic mean-field fluid that problem is ill-posed (dp/dn < 0 ⇒ elliptic regions), so it cannot be posed with this EOS.
- Not done and necessary: the collective-mode dispersion of the interacting gas (growth vs Hubble), a non-barotropic microscopic treatment (e.g. quantum field theory with the contact term beyond HF, trans-Planckian k_F ≈ 1.5–3 M_Pl so likely outside the EFT), and an inhomogeneous core–envelope calculation. No junction or causal claim is made.

## Comparison with Phase 3C.5

ECD removes the transport-law problem (spin dynamics comes from the Dirac equation, not an imposed transport law), clarifies that polarized spin does not help (Pauli), and fixes B from first principles (B = λ/2 … 2λ). It does **not** remove the Phase 3C pathology: the stress is still of the form ρ = ε − Bn², p = p − Bn², and the theorem above shows that structure forces μ = 0 and negative compressibility before the bounce.

## Gate decision for 3D.0

| criterion | result |
|---|---|
| consistent microscopic spin dynamics | **yes** (Dirac equation; only the axial torsion couples; λ = 3κ/16) |
| regular candidate transition | **yes in the homogeneous mean-field model** (finite R̃, nonsingular, exact baseline) |
| controlled constraints and curvature | **yes** (residual ~10⁻¹¹; curvature finite but Planckian, k_F ≈ 2–3 M_Pl) |
| no immediate fatal instability | **no**: negative compressibility from 0.105 n_b (theorem); dynamic rate not computed |
| plausible route to inhomogeneous collapse | **not with this model** (barotropic fluid ill-posed; mean-field invalid at the bounce) |
| **Decision** | **NO-GO for the stable mean-field ECD bounce. Failure origin: the matter state (barotropic mean field) and the theory's contact interaction (attractive in the density channel), not a symmetry assumption.** |

## What remains open (and is specific)

1. The RPA/Landau collective-mode dispersion of the degenerate gas with the Fierz-decomposed contact interaction: does the spinodal instability grow faster than the Hubble rate, and is it cut off by a finite unstable band?
2. A genuinely non-barotropic state (e.g. the Dirac field's own quantum dynamics beyond HF) at k_F ≈ M_Pl, where the EFT itself is not controlled.
3. Torsion theories with a dynamical torsion sector (Holst/Nieh–Yan, Poincaré gauge) where the contact interaction is not the whole story.

## Limits

First-order HF; massless/NR limits only (no intermediate mass scan); antiparticles ignored; the sign of the Fock part for the single-helicity gas follows the matrix algebra, not an independent literature value; no collective-mode growth rate; the closed-ball embedding fact is geometric, not a full junction solution.
