# Wormhole Formation in Matter-Coupled Vacuum Transitions: Gravitational Collapse, Horizon Constraints, and the Limits of Cosmological Rebirth

Leon Sandler, Independent Researcher — leon.sandler@alumni.swinburne.edu — ORCID 0009-0007-4584-808X

Code, data and manuscript for the paper of the same title, prepared for **Classical and Quantum Gravity** (IOP Publishing; Paper).
Manuscript: `manuscript/CQG_Wormhole_Formation_Vacuum_Transitions.{docx,pdf}`, supplementary information `manuscript/CQG_Supplementary_Information.{docx,pdf}`, cover letter `manuscript/CQG_Cover_Letter.docx`.
Archived on Zenodo: paper 10.5281/zenodo.24020378, code 10.5281/zenodo.24020381 (version 1.2.0; all versions: 10.5281/zenodo.23995978 / 10.5281/zenodo.23995982). Earlier versions: 1.1.2 (10.5281/zenodo.24020303 / 10.5281/zenodo.24020304); 1.1.1 (10.5281/zenodo.24012315 / 10.5281/zenodo.24012317); 1.1.0 (10.5281/zenodo.24011178 / 10.5281/zenodo.24011179); 1.0.0, first release with the earlier title (10.5281/zenodo.23995979 / 10.5281/zenodo.23995983).

## What the paper does

It asks whether a density-triggered vacuum transition can form a wormhole, a spacetime throat to a new expanding region, in a model whose only matter is a canonical scalar field coupled to dark
matter (false vacuum = observed dark-energy density, true vacuum below it). It shows that (i) the null energy condition holds, so no null ray has a minimum of the areal radius and no traversable throat
forms; (ii) a finite-tension thin-wall (Israel) analysis gives a critical exterior mass m_c(σ) that equals (1+q) R_dS up to corrections of order σ², where the exterior horizon and the interior de Sitter horizon have
equal area (a geometric coincidence of the small-tension limit, not a thermodynamic criterion); (iii) a spherically symmetric Einstein–scalar solver in double-null coordinates with a regular centre, validated against
exact solutions and by the convergence of its constraints, identifies trapped surfaces in 65 of 72 runs at the base resolution and in six more after refinement (one remains unresolved); (iv) a persistent expanding false-vacuum core appears only where the initial data already contain a false-vacuum core larger than its own
de Sitter radius, and in none of the 42 runs that start with a smaller core does one appear, so the runs show persistence, not creation, of the expanding core; (v) at the observed dark-energy scale that region has a radius
of 5.3 Gpc and a vacuum-energy deficit of 1.1×10²³ solar masses, 4.4 times the matter inside such a sphere today, so matter is not transferred.

The numerical walls are far thicker and heavier than the physical wall (the physical wall is thinner by about 10⁻³⁰ and lighter by about 10⁻⁶⁰ in units of the de Sitter radius), the evolution is spherical and
classical, the dark matter is not simulated, no spacelike slice through the bridge is built from the numerical solutions, and the code does not follow the solution beyond the singularity. The results are conditional
statements about this model.

## Files and results

| Result | Paper | Script |
|---|---|---|
| Symbolic derivation of the Einstein–scalar equations in double-null gauge | 3, 5.1 | `derive_equations.py` |
| Double-null solver with a regular centre (second-order diamond scheme, constraint-integrated sigma near the centre, secant solve for the centre value) | 5.1 | `dnull.py` |
| Verification against exact de Sitter and anti-de Sitter solutions | 5.2, Table 3 | `test_dnull.py` |
| Collapse scan (72 runs, q = 1), termination causes, anti-trapped cores | 5.3, 5.4, Table S1 | `collapse_scan.py` |
| Scans with other exterior vacua (q = −0.5, 0, 2) | 5.5, Table S3 | `collapse_scan.py` → `collapse_scan_q.json` |
| Convergence, scan statistics, null-ray minima, outcome classes | 5.2–5.4 | `validate.py` |
| Constraint residuals, null-expansion convergence, classification and refinement of the unresolved runs, creation test, other vacua | 5.2–5.5, Tables S2, S4 | `validate_extra.py` |
| Finite-tension critical mass, exterior horizon roots, entropy-area coincidence, mass and size at the observed scale | 4, 6, 7, Tables 2, 4 | `thinwall.py` |
| Figures 1–6 and S1 | | `make_figures.py` |
| Normalization scales of the constraint residuals | 5.2, Table S4 | `update_constraint_scale.py` |
| Fifteen automated checks | | `tests/test_all.py` |
| Manuscript, supplement, cover letter | | `paper_build.py`, `build_cover_letter.py` |
| Reference verification (Crossref) | | `crossref_check.py` |

## Reproduction

Python 3.11 with the packages of `requirements.txt` (numpy, scipy, numba, matplotlib, sympy, python-docx); LibreOffice for the PDF conversion. From `code/`:

```
python derive_equations.py            # prints the Einstein tensor components (about a minute)
python collapse_scan.py               # 72 + 42 runs, about 25 seconds (writes collapse_scan.json, collapse_scan_q.json)
python thinwall.py                    # writes thinwall.json
python validate.py                    # writes validation.json (about 2 minutes)
python validate_extra.py              # writes validation_extra.json (about 10 minutes)
python make_figures.py                # writes ../figures/Fig1-Fig6, FigS1
python tests/test_all.py              # 15 checks
python paper_build.py && python build_cover_letter.py
```

Every number in the manuscript is read from the JSON files these scripts write.

## Version 1.2.0

Adds Section 9.1 (primordial bridges and observationally hidden continuations, with Table 5 of success levels) and a causal-accessibility statement for the numerical solutions (no future-directed causal curve leads from the exterior into the expanding core). No change to the simulations. `docs/PRIMORDIAL_BRIDGE_PLAN.md` is a work plan for the separate investigation of pre-existing bridges (plan only, nothing computed).

## Version 1.1.2

Manuscript wording and one added validation field: the critical mass is described as dominated by the vacuum-energy difference (matter is 23% of it, and is not needed); constraint residuals are given in absolute and relative form (`update_constraint_scale.py` adds the normalization scales to `validation_extra.json`); the conclusion qualifies the threshold; the abstract states that the numerical evidence is complete only for anti-de Sitter exteriors.

## Version 1.1.1

Manuscript wording only (no change to code or results): the 71/72 result is stated as 65 at the base resolution plus six after refinement; the null-energy no-go is limited to ordinary traversable throats in the classical model; the numerical result is called a persistent expanding false-vacuum core; the limits of the q = 0 and q = −0.5 comparisons are stated.

## Revision 1.1.0

Retitled; the critical mass is the finite-tension maximum of the junction potential (its small-tension limit is (1+q) R_dS) and the horizon-area equality is stated as a coincidence of that limit; persistence of an
expanding core is separated from creation (0 of 42 runs create one); exteriors of other vacuum depth are tested; the constraints and null expansions of the dynamical solutions are shown to converge; the seven runs
without a trapped surface are classified and six resolved by refinement; the thin-wall bridge is stated precisely (the minimal sphere is the exterior's bifurcation sphere) and no throat is claimed from the numerics;
the omission of dark matter from the numerical experiments is stated.

## Citation

Please cite the manuscript and the code by their Zenodo DOIs above. Licence: MIT (code), CC-BY 4.0 (manuscript).
