# Phase 3G.0: damped oscillations after a high-energy transition ("nine diminishing bounces") (branch `phase3-recycling`)

Status: **stopped after 3G.0. Verdict: Direction B (post-transition cosmological ringdown) is a sound, calculable phenomenon, but it does not contain bounces; Direction A (repeated universe-wide bounces) can reach nine or more consecutive bounces in the Phase 1B toy model, but they are not diminishing, they end singular, and they inherit the unsolved source problem.** The pasted proposal left its section 6 ("Recommended Phase 3G.0") empty, so I took the natural reading of the surrounding text: arithmetic of a damped sequence, a flat-universe ringdown, a repeated-bounce test on the Phase 1B model, and scale constraints. 39 automated checks pass.

Reproduce (`code/`): `phase3g0_ringdown.py` → `phase3g0_ringdown.json` (≈40 s), `tests/test_all.py`. κ = 8πG = 1.

## A. The damped sequence

For successive peaks decreasing by q = 0.7: γ/ω_d = −ln q/(2π) = 0.0568, the ninth peak is q⁸ = 5.76% of the first (the 5.8% quoted). The number of peaks before falling below a threshold is 1 + ln(1/θ)/(−ln q): 7.5 (10%), 9.4 (5%), 13.9 (1%), 20.4 (0.1%): **nine is a consequence of the threshold, not a preferred number.**

## B. Direction B: ringdown in an expanding flat universe

Model: GR + massive scalar (m) with decay width Γ into radiation, φ̈ + (3H + Γ)φ̇ + m²φ = 0, ρ̇_r + 4Hρ_r = Γφ̇², H² = ρ/3, starting at H₀/m = 0.0757 (the value at which Hubble friction alone gives 0.7 between successive extrema of |φ|).

1. **No bounce is possible.** Ḣ = −(1/2)(ρ+p) ≤ 0 with ρ > 0: H decreases monotonically and never changes sign (max Ḣ = −1.2×10⁻¹¹ over the run; zero H-crossings in all three cases). The oscillation is a ripple in ä/a. The cycle-by-cycle change in H is large at first (43% in the first cycle) but secular, never a reversal.
2. **"Bounces" can be counted only as acceleration episodes (ä > 0).** Hubble friction only: 28 episodes in 14 cycles (two per cycle, because ρ + 3p changes sign twice per cycle); Γ/m = 0.03: 9 episodes; Γ/m = 0.2267: **0 episodes** (the decay products make radiation dominate, w → 1/3, ä < 0 for ever). So dissipation both damps and extinguishes the ripples in acceleration.
3. **The amplitude sequence is not geometric with Hubble friction.** Successive extremum ratios: 0.79, 0.83, 0.85, 0.87, 0.89, 0.90 … 0.966 (→ 1): amplitude ∝ a^{−3/2} ∝ 1/t, A₉/A₁ = 0.32, not 0.058. A geometric sequence needs a constant damping rate: with Γ/m = 0.2267 (extremum ratio e^{−πΓ/2m} = 0.7) the ratios rise from 0.55 toward 0.68 as the Hubble friction falls (A₉/A₁ = 0.020; the asymptotic ratio 0.70 is reached only when H ≪ Γ); with Γ/m = 0.03: ratios 0.76 → 0.93, A₉/A₁ = 0.23.

So the "ball bouncing nine times" corresponds, in an expanding universe, to a scalar ringing with a damping that is Hubble friction plus a decay width, and its energy goes to radiation (87% radiation at the end for Γ/m = 0.03; 100% for Γ/m = 0.2267), which is covariantly conserved with entropy produced: the standard picture of reheating.

**Extra-dimensional stabilization.** The oscillating modulus case is the same mathematics (an effective scalar with a potential); no separate calculation was done, and oscillating compact dimensions do not make the three large dimensions bounce. In the Phase 3E.0 truncation the volume modulus was frozen.

## C. Direction A: repeated universe-wide bounces (Phase 1B closed scalar model, Γ = 0.3)

Re-running the seed scan of the Phase 1B family (μ = 20, ρ_c = 0.2V_f, C = 0.3C_c), the best seed (log₁₀ = −0.915) gives **12 consecutive bounces**, so "nine or more" is attainable transiently. But:
- turnaround radii a_max: 0.96, 0.94, 0.93, 0.90, 0.88, 0.92, 0.91, 0.88, 0.91, 0.91, 0.90, 0.94 (no decrease, log-slope +0.005 per bounce);
- neck radii a_min: 0.918 → 0.831 (decreasing, log-slope −0.0085 per bounce);
- oscillation amplitude a_max − a_min: 0.045, 0.027, 0.038, 0.022, 0.011, 0.005, 0.052, 0.001, 0.062, 0.065, 0.061, 0.104: **irregular and growing at the end (log-slope +0.065 per bounce), the opposite of a damped ball**;
- entropy parameter S = ρ_r a⁴ increases monotonically every bounce, from 0.07 to 0.35 (a factor 5), the Tolman effect: dissipation raises, not lowers, the amplitude;
- it ends **singular at t = 14.3**. The stored lifetime scan shows maxima of 13–15 bounces for Γ = 0.2–0.7 and 0–1 for Γ ≥ 1, all ending singular.

The first-principles reason is that radiation accumulates (it is covariantly conserved and cannot be deleted), so each cycle starts with more entropy; there is no damping toward a smaller bounce, and the system tends toward collapse. This repeats the Phase 1B conclusion (no stable cold cycle with dissipation) in the "nine bounces" language. Repeated bounces also need the closed-universe, tuned-seed, bounce-source assumptions of Phases 1–2 and the unsolved gravitational source (3D.3–3F.0).

## D. Scale and observational constraints

Standard cosmology requires the ringing field to have decayed before big-bang nucleosynthesis: H(T = 1 MeV) = 4.5×10⁻²⁵ GeV (t ≈ 0.7 s), so Γ ≳ H_BBN. The number of oscillations before decay is ≈ (m/πΓ)·ln(1/0.05): for m = 1 TeV and Γ ~ H_BBN this is ~10²⁷, and for a field decaying at Γ ~ 0.1m about 10: the number of oscillations is set by m/Γ, not by a universal nine. Unobserved oscillations in the expansion history after BBN or in the perturbation spectrum would constrain any model; this was not computed.

## Decision

| question | answer |
|---|---|
| does a post-transition ringdown contain genuine bounces? | **no**: H never reverses; the signal is ripples in ä/a, at most two acceleration episodes per cycle, none once radiation dominates |
| is the decay a geometric 0.7ⁿ sequence? | **no with Hubble friction alone** (power law, ratios → 1); geometric only for H ≪ Γ |
| can repeated universe-wide bounces diminish like a ball? | **no** in the Phase 1B model: 12 bounces reached, but amplitude not decreasing, entropy ×5, necks shrinking, ends singular |
| does this solve black-hole matter recycling? | **no**: it neither provides the gravitational source nor the causal continuation |
| what is it useful for? | **Direction B** is a well-defined early-universe calculation (reheating / modulus ringdown) once a transition exists |

## Limits

Flat, homogeneous, minimally coupled GR with a quadratic potential and a phenomenological decay width; one starting H₀/m; no perturbations; the Phase 1B model is a closed scalar-potential toy whose bounces are classical (not derived from a fundamental theory); no CMB/BBN data were fitted; the extra-dimensional case was not calculated.

## Recommended next step

If the cosmological ringdown idea is pursued, it should be as a separate early-universe study with an explicit field content (inflaton or modulus, decay channels, reheating) and observational checks, not as a source for black-hole recycling. For the recycling programme the open item is unchanged: the S6 model of the 3D.3 census.
