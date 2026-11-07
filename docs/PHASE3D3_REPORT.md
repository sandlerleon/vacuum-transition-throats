# Phase 3D.3: Goldilocks feasibility and source search (branch `phase3-recycling`)

Status: **stopped before any parameter scan, as the instructions require when no admissible source can be established. Verdict: NO-GO. No candidate source satisfies the hard constraints; the one competition that can be closed quantitatively (free quantum vacuum, supply versus effective-field-theory validity) has an empty window for every throat size, with a scale-free margin of at least 5.6.** 34 automated checks pass.

Reproduce (`code/`): `phase3d3_goldilocks.py` → `phase3d3_goldilocks.json` (4 s), `tests/test_all.py`. Planck units, κ = 8π.

Method note. The attached superconductivity paper was not available to me in this session; I used only the method as described in your message (two competing effects, robustness against disorder, falsification boundaries, hard constraints rather than a weighted score). No result of that paper is used.

## Step 1. Hard constraints (all must hold; no score can offset a failure)

H1 satisfies the field equations; H2 any NEC-violating sector has a derivable action; H3 no ghost (Hamiltonian bounded below); H4 perturbatively stable; H5 inside effective-field-theory (EFT) validity; H6 arises from collapse of ordinary low-density matter; H7 future-directed continuation into an expanding region, not merely geodesic crossing.

## Step 2. What the benchmark forces on its source (derived here)

For the Simpson–Visser metric, G^r_r − G^t_t = −2f a²/R⁴ exactly (sympy residual 0). A minimally coupled scalar with ε = ±1 contributes T^r_r − T^t_t = ε f φ′², and a purely radial magnetic nonlinear-electrodynamics (NED) field contributes zero. Hence κ ε φ′² = −2a²/R⁴: **ε = −1 is forced (a phantom/ghost scalar) with φ = √(2/κ) arctan(r/a)** (residual 2×10⁻¹⁶). A canonical scalar cannot do it. This reproduces the Bronnikov–Walia result that the benchmark is sourced by a phantom scalar with potential plus NED (Phys. Rev. D 105, 044039; search-verified, arXiv:2112.13198) and shows the ghost is not a choice for this metric: it is forced within scalar + radial-NED sources. The "trapped ghost" variant (ghost only in a strong-field region) applies to other metrics, not to this one, where T_kk < 0 at every r.

## Step 3. Collapse baseline

Taken from Phase 3A (not recomputed): with ordinary matter the NEC holds, so by the Penrose theorem trapped surfaces end in singularities (65/72 double-null runs trapped, 0/42 expanding regions created, 96 LTB shells all trapped). A regular transition needs an NEC-violating sector that is generated, not imposed (Phase 3D.2).

## Step 4–5. Census of candidate sources and their competing effects

| source | H1 | H2 | H3 | H4 | H5 | H6 | H7 | decisive reason |
|---|---|---|---|---|---|---|---|---|
| S1 ordinary matter / canonical scalar | pass | pass | pass | pass | pass | pass | **FAIL** | NEC ⇒ singularity (Penrose); 3A |
| S2 phantom scalar + NED (forced source of the benchmark) | pass | pass | **FAIL** | **FAIL** | open | **FAIL** | pass | ghost forced; static Ellis–Bronnikov-type wormholes carry an unstable radial mode (literature, secondary sources) |
| S3 NED regular black hole (dS/Minkowski core) | pass | pass | pass | open | pass | open | **FAIL** | no bounce into an expanding region; inner-horizon (mass-inflation) instability reported generically (Carballo-Rubio et al., 2021) |
| S4 Einstein–Cartan–Dirac, barotropic mean field | pass | pass | pass | **FAIL** | **FAIL** | open | open | 3D.0: negative compressibility from 0.105 n_b; spinodal at k_F = 1.45 M_P (1.02 for a single-helicity gas), bounce at 2.2–3.1 M_P |
| S5 free-field quantum vacuum | open | pass | pass | open | **FAIL** | open | open | Part B below |
| S6 limiting-curvature / mimetic / Horndeski | open | pass | open | open | open | open | open | homogeneous dust case reproduces the 3B effective fluid (already pathological); inhomogeneous collapse and stability not evaluated here |

Exactly one entry is wholly unevaluated (S6): this is a limitation of this phase, not a proof of failure, and it is the natural place for a continued source search. **No source passes all of H1–H7.**

### Part B. The Goldilocks competition for the quantum source (the one computable window)

Two competing requirements, both scaling as a²:
- **Supply.** A throat of radius a needs |ρ + p_r| ≈ c/(κa²). A Ford–Roman-type bound on free fields gives a supply ≤ 3N/(32π²a⁴), so **N ≥ 5.59 a²** (a/m = 6 coefficient from 3D.1; 8.38 a² as a/m → ∞; ≈1000 a² for a typical state with the Page/DeWitt–Schwinger magnitude N/(1920π²a⁴)).
- **Validity.** With N species the gravitational cutoff drops to Λ ≈ M_P/√N (Dvali 2007; Dvali–Redi 2008; conjecture-level, search-verified), and the throat curvature 1/a must lie below it: **N ≤ a²** (G = 1).

The window needs 5.59 a² ≤ N ≤ a²: **empty for every a, by a factor 5.59, independent of the throat size.** On a grid of 361 301 points (a/ℓ_P from 1 to 10⁴⁰, N from 1 to 10⁹⁰) the number of admissible points is **0** for all six constant choices (margins 5.6, 8.4, 1000, and 140, 211, 25 000 with the reduced Planck mass). Robustness: allowing each O(1) constant to vary log-uniformly over a factor 3 either way, only 2.4% (a/m = 6) and 0.06% (a/m → ∞) of draws open a window; the prior is arbitrary and this only measures the slack that would be needed. Caveats: flat-space free-field bound with a heuristic sampling time τ₀ ≲ a; the species bound is a conjecture supported by black-hole arguments with an O(1) (and ln N) uncertainty; no actual state's ⟨T_μν⟩_ren was computed. The bound could fail for strongly interacting or non-free fields, which the census leaves open.

## Step 6. Robustness and falsification

- No fine-tuned or narrow window exists in the computable competition; the failure is scale-free, so it cannot be cured by choosing a throat size a/m, a rotation or a magnetic field (the derived source requirement T_kk = −2a²/(κR⁴) is independent of m and of the observer).
- Falsification conditions met: required matter source not derivable except as a ghost (S2); equations used beyond the EFT cutoff (S4, S5); formation from ordinary matter excluded (S1).
- Not met/needed: curvature is regular and geodesic continuation exists (3D.2) in the benchmark, but that is the benchmark's property, not evidence for a source.

## Decision gate

| deliverable | status |
|---|---|
| derived source equations | the forced ghost sector (Step 2); other sources from earlier phases |
| parameter-space feasibility map | quantum-vacuum map (Part B): empty; others are not parametrically scannable because they fail a hard constraint first |
| explicit Goldilocks boundaries | N = 5.59 a² (supply) and N = a² (validity): they do not intersect |
| robustness | slack factor ≥ 5.6 needed; 2.4% of arbitrary-prior draws |
| causal classification of surviving solutions | none survive; benchmark's causal structure in 3D.2 |
| **go/no-go** | **NO-GO** for a physically admissible Goldilocks window with the sources examined; stop before a large scan |

## What this does and does not exclude

It does not prove recycling impossible: S6 (modified gravity with a derivable action), strongly interacting fields and fully dynamical torsion are unevaluated, and the species-bound argument is heuristic. It excludes (with the stated caveats) the three simplest routes: ghost-supported throats, free-field vacuum support, and barotropic ECD bounce.

## Recommended next step

If the programme continues, the only open source class with a derivable action and no immediate failure is S6. A bounded next phase: take one concrete limiting-curvature or mimetic model, derive its spherically symmetric field equations and the collapse of ordinary dust inside a horizon, and check H3–H7 (ghost, stability, formation, expanding continuation) before any parameter scan.
