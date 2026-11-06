# Plan for a separate investigation: pre-existing (primordial) non-traversable bridges

Status: **plan only** — nothing in this document has been computed. The manuscript of this repository (Section 9.1) states the hypothesis and the criteria; this file turns them into a work programme.

## Hypothesis (stated so that it can fail)

A pre-existing, non-traversable spacetime bridge connects a late-time vacuum-transition region to an expanding cosmological region, with the connection hidden from observers in the parent universe.
It needs no null-energy violation and no spontaneous formation of a throat at the end of the cosmic evolution. Being unobservable is **not** a test; the tests are the predictions below.

## What the present paper already gives

* A null-energy argument: no ordinary traversable throat in the classical canonical-scalar + non-negative dust model, however the throat arose.
* A finite-tension thin-wall junction analysis, including the bridge branch (interior beyond the de Sitter equator, Γ_in < 0).
* A double-null solver with a regular centre and a causal-accessibility statement for its solutions: u and v do not decrease along future-directed causal curves, the core is at v < V_w and the exterior at v > V_w, so no causal curve enters the core from the exterior.
* A persistence/creation test (0 of 42 normal-core runs create an expanding region).

## Work programme

| Stage | Question | Method | A result that would count as failure |
|---|---|---|---|
| 1 (Level 1) | What is the global causal structure of the thin-wall bridge branch? | Glue the maximally extended Schwarzschild–(A)dS exterior to the de Sitter cap along the exact wall trajectory R(τ) of Section 4; propagate radial null geodesics across the wall by junction matching; map the Penrose diagram | A future-directed causal curve from each asymptotic region to the cap does not exist (then the cap is causally inaccessible from the parent: Level 3 fails for this bridge) |
| 2 (Model B) | Does a pre-existing bridge survive a density-triggered vacuum transition? | Double-null evolution from two-ended characteristic data through the bifurcation sphere (no regular centre: the centre treatment of dnull.py is replaced by data on two crossing null rays); a false-vacuum region with a wall; same diagnostics (trapped, anti-trapped regions, constraint and expansion convergence) | The bridge pinches off or a trapped surface closes the connection for every tested parameter; or no anti-trapped region survives |
| 3 (Model A control) | What do regular localized data do? | The present paper | (done: collapse without creation of an expanding region) |
| 4 (Model C) | Can a stress-energy correction change the causal behaviour? | Only after stages 1–2: an explicitly defined correction (semiclassical or modified gravity) with stated magnitude, duration and consistency checks | The required negative energy exceeds what the correction can supply, or the geometry is unstable |

## Levels of success (from Section 9.1, Table 5)

1. A mathematically consistent hidden bridge exists — stage 1.
2. The bridge survives a density-triggered vacuum transition — stage 2.
3. A causal path carries energy or matter into the expanding region — stage 1 criterion, then stage 2.
4. The expanding region develops a hot cosmology — needs a reheating mechanism (new physics).
5. The process repeats under the same laws — needs a global construction.

Levels 3–5 require new physics or a global spacetime analysis; levels 1–2 are within reach of the existing tools.

## Predictions that could fail (to be derived, not assumed)

* A stability boundary in (wall tension, exterior mass, core radius) for the bridge.
* A causal-accessibility criterion: the sign of a null-geodesic crossing condition at the wall.
* A quantitative distinction from ordinary black-hole collapse (for example a different late-time r_v, r_u structure in the region behind the horizon).

## Decisions for the author

1. Whether to write this as a second paper (recommended) or fold stage 1 into a revision of the present one (not recommended: the present paper is a finished classical no-go).
2. Whether stage 2 should use the same scalar potential family (recommended) or a different false-vacuum model.
3. Whether to include stage 4 at all: it introduces assumptions that the null-energy argument cannot test.
