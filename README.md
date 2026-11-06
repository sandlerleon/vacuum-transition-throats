# Vacuum-Transition-Induced Spacetime Throats: Gravitational Collapse, Matter Transfer, and Conditional Cosmological Rebirth

Leon Sandler, Independent Researcher — leon.sandler@alumni.swinburne.edu — ORCID 0009-0007-4584-808X

Code, data and manuscript for the paper of the same title, prepared for **Classical and Quantum Gravity** (IOP Publishing; Paper).
Manuscript: `manuscript/CQG_Vacuum_Transition_Throats.{docx,pdf}`, supplementary information `manuscript/CQG_Supplementary_Information.{docx,pdf}`, cover letter `manuscript/CQG_Cover_Letter.docx`.
Archived on Zenodo: paper 10.5281/zenodo.23995979, code 10.5281/zenodo.23995983 (version 1.0.0; all versions: 10.5281/zenodo.23995978 / 10.5281/zenodo.23995982).

## What the paper does

It asks whether a density-triggered vacuum transition can produce a spacetime throat to a new expanding region, in a model whose only matter is a canonical scalar field coupled to dark
matter (false vacuum = observed dark-energy density, true vacuum below it). It shows that (i) the null energy condition holds, so no null ray has a minimum of the areal radius and no traversable throat
forms; (ii) an exact thin-wall (Israel) criterion gives a critical exterior mass m_c = (1+q) R_dS, at which the exterior horizon and the interior de Sitter horizon have equal area; (iii) a spherically
symmetric Einstein–scalar solver in double-null coordinates with a regular centre, validated against exact solutions, forms trapped surfaces and an expanding core only when the false-vacuum region
extends beyond its own de Sitter radius; (iv) at the observed dark-energy scale that region has a radius of 5.3 Gpc and a vacuum-energy deficit of 1.1×10²³ solar masses, 4.4 times the matter inside
such a sphere today, so localized collapse cannot produce a daughter and matter is not transferred to it.

The numerical walls are far thicker and heavier than the physical wall (the physical wall is thinner by about 10⁻³⁰ and lighter by about 10⁻⁶⁰ in units of the de Sitter radius), the evolution is spherical and
classical, and the code does not follow the solution beyond the singularity. The results are conditional statements about this model.

## Files and results

| Result | Paper | Script |
|---|---|---|
| Symbolic derivation of the Einstein–scalar equations in double-null gauge | 3, 5.1 | `derive_equations.py` |
| Double-null solver with a regular centre (second-order diamond scheme, constraint-integrated sigma near the centre, secant solve for the centre value) | 5.1 | `dnull.py` |
| Verification against exact de Sitter and anti-de Sitter solutions | 5.2, Table 2 | `test_dnull.py` |
| Collapse scan (72 runs), horizon radius against the analytic value, anti-trapped cores | 5.3, 5.4, Table S1 | `collapse_scan.py` |
| Convergence, scan statistics, null-ray minima, outcome classes | 5.2–5.4 | `validate.py` |
| Exact thin-wall criterion, critical mass, entropy equality, mass and size at the observed scale | 4, 6, 7, Table 3 | `thinwall.py` |
| Figures 1–5 | | `make_figures.py` |
| Nine automated checks | | `tests/test_all.py` |
| Manuscript, supplement, cover letter | | `paper_build.py`, `build_cover_letter.py` |
| Reference verification (Crossref) | | `crossref_check.py` |

## Reproduction

Python 3.11 with the packages of `requirements.txt` (numpy, scipy, numba, matplotlib, sympy, python-docx); LibreOffice for the PDF conversion. From `code/`:

```
python derive_equations.py            # prints the Einstein tensor components (about a minute)
python collapse_scan.py               # 72 runs, about 10 seconds (writes collapse_scan.json)
python thinwall.py                    # writes thinwall.json
python validate.py                    # writes validation.json (about 2 minutes)
python make_figures.py                # writes ../figures/Fig1-Fig5
python tests/test_all.py              # 9 checks
python paper_build.py && python build_cover_letter.py
```

Every number in the manuscript is read from the JSON files these scripts write.

## Citation

Please cite the manuscript and the code by their Zenodo DOIs above. Licence: MIT (code), CC-BY 4.0 (manuscript).
