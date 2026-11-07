# Phase 2E: validation of the nontraversable cyclic-bridge construction (branch `phase2-bridge`)

Status: **the construction is locally and globally consistent as an idealized exact solution, but it is NOT dynamically robust in this model: linear inhomogeneous perturbations grow by factors of 14 to 10⁵ per cycle.** No traversable wormhole, matter transfer or dynamical formation is claimed. 23 automated checks pass. The CQG manuscript and the Phase 1 documents are untouched.

Reproduce (`code/`): `phase2e_junction_global.py`, `phase2e_perturb_derive.py`, `phase2e_perturb_reduce.py`, `phase2e_perturb_floquet.py`, `phase2e_perturb_mechanism.py`, `phase2e_window.py`, `phase2e_coupling_scan.py <rc> <mu>`, `make_phase2e_figure.py` (writes `figures/FigE1_phase2e.png`), `tests/test_all.py`.

## Demonstrated versus unresolved

| claim | status |
|---|---|
| the cyclic region is an exact, independently reproduced, homogeneously stable (elliptic) solution | demonstrated numerically |
| local junction conditions (induced metric, normals, extrinsic curvatures, stress tensor, orientation, time orientation) | demonstrated from raw Christoffel symbols, residuals ≤ 5×10⁻¹⁵ |
| global gluing: kept exterior = I ∪ II ∪ IV ∪ (III with R < R_s), one timelike boundary, shell stays outside the horizon | demonstrated (construction + checks) |
| cap and parent universe are causally disconnected; the bridge is non-traversable | demonstrated (exact Kruskal argument + rays) |
| shell satisfies σ > 0 and the NEC, periodic, elliptic return map, survives 60 periods at ≤ 1% perturbation (homogeneous cap, spherical shell perturbations) | demonstrated numerically |
| stability under inhomogeneous perturbations | **fails** (linear analysis, §3) |
| shell matter from an action | EOS p = −σ/2 follows from an isotropic string gas; **the required energy exchange with the bulk is unexplained** (§4) |
| independence of the eternal white-hole region | the future development of data on a slice through the throat does not contain region IV; the construction still needs a pre-existing throat (§5) |
| dynamical formation, traversability, matter transfer | not demonstrated; the last two are excluded in this geometry |

## 1. Test A: global junction consistency

`phase2e_junction_global.py` recomputes everything from each side's own coordinates: interior g = diag(−1, a², a² sin²χ, …) on (t, χ), exterior f = 1 − m/R on (s, R) in region III with future-directed Killing time s. Unit normals point from the cap into the kept exterior; extrinsic curvatures K_ab = −Γ-terms contracted with n (K_ττ = −n_ν A^ν, K_θθ = −Γ^μ_θθ n_μ), junction S^a_b = −([K^a_b] − δ^a_b[K]). Along one full period of the conservative periodic shell (m/m_max = 0.10, 3001 points):

- induced metric continuity: h_ττ = −1 on both sides to 6×10⁻¹⁶, h_θθ equal to 6×10⁻¹⁷; normals unit and orthogonal to u to 7×10⁻¹⁶;
- orientation: n^χ > 0 inside, n^R < 0 outside (the kept side contains the throat), dt/dτ > 0, ds/dτ > 0, R > r_h at all times;
- σ and p recovered from raw curvatures equal those of the equation of motion to 5×10⁻¹⁵, and p = wσ with w = −1/2 to 4×10⁻¹⁵ (an independent test of the equation of motion);
- the kept exterior contains no boundary other than the shell: regions I, II, IV are complete Kruskal regions and the part of III with R < R_s is bounded by the shell and by the horizons that continue into II and IV;
- curvature: exterior Kretschmann scalar 0.17–0.61 at the shell, **3.8×10¹¹ at the throat** (r_h = 0.0033 in units R_dS = 1), interior 15–24.

The throat is therefore extremely small and strongly curved compared with the cycle scale.

## 2. Test E (item 5): the stability window

Independent code (Phase 2D). Window C/Cc ∈ (0.3458, ≈ 0.360). Nonlinear survival of the neck map (150 iterates, 2 seeds, perturbations 10⁻³, 10⁻², 3×10⁻²): survives at C/Cc = 0.3470 (near the multiplier = −1 point), 0.3500, **0.351358 (rotation number exactly 2/5)** and **0.354887 (3/8)**, with neck deviation ≈ 4–6 × the perturbation and no growth. Near the end of the window (C/Cc = 0.3598) the basin shrinks sharply: at 10⁻² one of two starts is lost immediately, at 3×10⁻² both. A scan of the dust–scalar coupling (ρ_crit/V_f = 0.2, 0.4, 1, 2, 4) in steps of 0.05 in C/Cc found an elliptic cycle only at ρ_crit/V_f = 0.2 (and none at 2 and 4, where no cold cycle was found); the elliptic window is only 0.014 wide in C/Cc, so this scan cannot exclude narrow windows elsewhere.

## 3. Test B: inhomogeneous perturbations

**Derivation.** Spherically symmetric linear perturbations of the closed cycle with scalar field and dust coupled by the mass function A(φ) = 1 + φ²/2M², in longitudinal gauge (Ψ = Φ follows from the traceless equation). The Einstein tensor, the scalar equation □φ = V′ + A′ρ_*, dust number conservation and the dust force law u^ν∇_ν u_μ = −h_μ^ν ∂_ν ln A are expanded to first order with sympy (`phase2e_perturb_derive.py`) and projected on radial harmonics Q_n = sin(nχ)/sin χ, k² = n² − 1. The constraints (00) and (0χ) give Φ algebraically with no division by H, so the equations are regular at the bounces. n = 1 is the homogeneous mode; **n = 2 is not covered** (longitudinal gauge degenerates, k² = 3K). Checks: constraint propagation (d/dt of the algebraic Φ against Φ′ from (0χ)): 6×10⁻⁹ and 3×10⁻⁸ relative (n = 5, 12); the unused trace equation: 2×10⁻⁶ and 1×10⁻⁶.

**Result** (`phase2e_perturb_floquet.json`): for every harmonic n = 3…60 the 4×4 monodromy over one cycle has a multiplier far above 1 (and its inverse): n = 3: 14; 5: 89; 10: 1.1×10³; 20: 1.2×10⁴; 60: 1.5×10⁵ (det = 1). **A 10⁻⁶ perturbation reaches order unity within about 1.2–3 cycles for n ≥ 6 and within about 5 cycles at n = 3.**

**Mechanism.** At large n the growth is k-independent: in the quasi-static limit the dust contrast obeys D̈ = g²D with g² = 6Cφ₀²/(M²a³(2M² + φ₀²)), the attractive scalar-mediated (fifth) force on the coupled dust. Integrated over the cycle, g gives a log-growth of 12.6 (rate 18–23 per unit time) against ln(multiplier) = 11.95 at n = 60 and 11.3 at n = 40 (the asymptote of the full result). With the dust–scalar coupling switched off **in the perturbations only**, the large-n multipliers fall to 1.34 (n = 30, 60), while n = 3 and 10 remain unstable (334 and 98): the field sits near the top of its potential, where long-wavelength perturbations are tachyonic. So two mechanisms act: a coupling-driven instability at all scales and a hilltop (spinodal) instability at long wavelength.

**Limits of this result.** Linear, on the full S³ (the cap has a boundary at the shell whose boundary conditions were not imposed), dust treated as pressureless with no velocity dispersion, one parameter set (μ = 20, ρ_crit = 0.2 V_f, q = 1). A pressure or finite-range force could weaken the coupling-driven part; nothing here shows it. The nonlinear inhomogeneous evolution (double-null with a dust fluid and a moving boundary) was not performed; the linear growth rates make a robust survival implausible without a new stabilizing mechanism.

## 4. Test C: shell matter from an action

An isotropic gas of non-interacting, slowly moving Nambu–Goto strings confined to the two-dimensional shell has energy density σ = μn (n = length per area) and tension μ along each string, so S^{ij} = −μn ⟨t̂^i t̂^j⟩ = −(σ/2)δ^{ij}: **p = −σ/2 exactly** (in two spatial dimensions ⟨t̂^i t̂^j⟩ = δ^{ij}/2). Without exchange the string length would be conserved, giving (σR)′ = 0. But the junction conditions require, for any shell moving through the homogeneous interior, (σR²)′ + p(R²)′ = −R²(ρ + P) cosh β sinh β (verified to 10⁻⁶ in Phase 2C): the shell must absorb and re-emit bulk dust and scalar energy. Numerically (`phase2e_junction_global.json`) the gross energy exchanged with the bulk per cycle is 0.54, which is **97% of the shell's own mean energy** (0.39–0.77) and 1% of the cap's energy; the net over a period is zero. A string gas with fixed tension μ cannot supply this unless strings are created and destroyed by absorbing dust and scalar quanta, and no such interaction has been specified. The effective equation of state is therefore physically motivated, but the shell as a self-contained string network is not a complete matter model; the energy bookkeeping is satisfied only by an unmodelled interaction.

## 5. Test D: dependence on the eternal region

Take the slice Σ₀ = {s = 0} in region III, continued through the bifurcation sphere into region I, with the cap data at the time the shell crosses s = 0. Σ₀ carries regular initial data (the cap satisfies the Friedmann constraint; the shell satisfies the Israel conditions; the exterior is time-symmetric vacuum). Its future domain of dependence contains I, II and the later part of III with the cap and shell, and **does not contain region IV**. So the evolution needs no white-hole region. What it does need is the initial throat: Σ₀ has a minimal sphere and two ends. This is a **primordial, pre-existing non-traversable bridge**; formation from regular collapse data (which has no left end and no region IV) is a different, unaddressed problem. The eternal extension to the past is only the maximal Cauchy development.

## 6. Recommendation

Do not present the cycle as a robust solution. A defensible statement: "an exactly homogeneous, independently reproduced cold cycle admits a consistent non-traversable two-ended embedding, but it is linearly unstable to inhomogeneous perturbations in this model". The paper draft (`manuscript/Nontraversable_ER_Bridge_Cyclic_Region.docx`, `.pdf`; built by `code/phase2e_paper_build.py` from the result files) is written on that basis.
