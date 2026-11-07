# -*- coding: utf-8 -*-
"""One-time editing script for revision 2 of the CQG manuscript sources (content_a/b/c.py, pcommon.py, build_cover_letter.py).  Every anchor must match exactly once, so that no edit is silent.   python apply_revision2.py"""
import io
import os

HERE = os.path.dirname(os.path.abspath(__file__))
BS = chr(92)


def rd(f):
    return io.open(os.path.join(HERE, f), encoding="utf-8").read()


def wr(f, s):
    io.open(os.path.join(HERE, f), "w", encoding="utf-8", newline="\n").write(s)


def sub(s, old, new, count=1):
    n = s.count(old)
    assert n == count, "anchor found %d times (expected %d): %s" % (n, count, old[:90])
    return s.replace(old, new)


# ------------------------------------------------------------------------------------------------ pcommon
s = rd("pcommon.py")
s = sub(s, 'VERSION = "1.1.0"', 'VERSION = "2.0.0"')
s = sub(s, 'DOIS = json.load(open(os.path.join(HERE, "dois.json"), encoding="utf-8"))',
        'DOIS = json.load(open(os.path.join(HERE, "dois.json"), encoding="utf-8"))\nAUD = json.load(open(os.path.join(HERE, "independent_audit.json"), encoding="utf-8"))\nWALLX = json.load(open(os.path.join(HERE, "wall_extrapolation.json"), encoding="utf-8"))\nCAUS = json.load(open(os.path.join(HERE, "causal_diagnostic.json"), encoding="utf-8"))')
wr("pcommon.py", s)

# ------------------------------------------------------------------------------------------------ cover letter email
s = rd("build_cover_letter.py")
s = s.replace("leon.sandler@alumni.swinburne.edu", "sandler.leon@gmail.com")
wr("build_cover_letter.py", s)

# ------------------------------------------------------------------------------------------------ content_a
s = rd("content_a.py")
s = sub(s, "Corresponding author: leon.sandler@alumni.swinburne.edu", "Corresponding author: sandler.leon@gmail.com")
i = s.index("ABSTRACT = (")
j = s.index('add(("p0", ABSTRACT))')
new_abs = '''ABSTRACT = (
    "Rebirth scenarios invoke a wormhole carrying matter into a new expanding region. We ask whether a density-triggered vacuum transition can form one in a scalar-only gravitational benchmark: a canonical scalar whose false vacuum has the observed "
    "dark-energy density, motivated by a coupling to dark matter that the numerical experiments do not evolve. In this classical model the null energy condition holds, so no ordinary traversable throat can form: no null ray has a "
    "minimum of the areal radius. Non-traversable bridges and quantum or baby-universe continuations are not excluded. For a thin spherical wall and regular initial data the Israel junction conditions give a finite-tension critical exterior mass parameter, "
    "equal to %s up to corrections of order $\\\\hat\\\\sigma^2$, below which the configurations on the branch analysed collapse; the equality of horizon areas there is a geometric coincidence of the limit. "
    "An independent derivation reproduces the junction branches, the horizon assignments and the dimensional restoration. A spherical Einstein–scalar code checked against exact solutions identified trapped surfaces in %d of %d runs "
    "at the base resolution and in %d more after refinement; one stayed unresolved. An expanding (anti-trapped) false-vacuum core appears only when the initial data already contain a core larger than its own de Sitter "
    "radius; none develops in the %d runs that start smaller, so the simulations show persistence, not creation. At the observed dark-energy scale that core has a radius of %s Gpc and the exterior mass parameter is "
    "%s solar masses (for $q=1$ the quasi-local mass at the wall is half of it), dominated by the vacuum-energy difference. In the evolved domain no causal curve carries exterior matter into the core; the global structure is not computed. "
    "The numerical evidence is complete only for anti-de Sitter exteriors and for walls much thicker than the physical wall."
    % ("$(1+q)R_{dS}$", SC["n_trapped"], SC["n_runs"], resolved, CRE["n_initially_normal_core"], n(R_dS, 1), "%s × 10$^{%d}$" % (n(M1["M_c_Msun"] / 10 ** 23, 1), 23)))
'''
s = s[:i] + new_abs + s[j:]
# terminology paragraph after the 'Three outcomes' paragraph
anchor = 'add(("p", "The mechanism is not new in its geometry.'
term = '''add(("p", "Four terms are used with fixed meanings. An anti-trapped (expanding) region is a region in which both null expansions are positive, $r_u>0$ and $r_v>0$: a local property that the numerical solutions can test. "
          "The expanding core is the false-vacuum interior of the wall when it contains such a region. A bridge is the non-traversable two-ended geometry of the thin-wall solution (Section 4); no bridge is constructed numerically. "
          "A daughter universe is reserved for a global expanding region that the collapse has created and that has a causal continuation; no run produces one. Persistence is the survival of an anti-trapped region that is present in the initial data, "
          "and creation is its appearance from initial data that have none."))
'''
s = sub(s, anchor, term + anchor)
# novelty paragraph and contributions
old = ("What is new here is the application to a specific density-triggered model with a fixed vacuum scale. The vacuum energy is the observed dark-energy density, so the critical sizes and masses \"\n"
       "          \"are numbers rather than free parameters, and they can be compared with the matter that the universe actually contains. The paper makes four contributions.\"))")
new = ("Known from that literature are the Israel junction analysis of a false-vacuum bubble, the existence of a critical mass below which regular data collapse, and the obstruction to reaching the baby-universe branch from regular initial data. "
       "What this paper adds is quantitative and specific to the model: a null-energy proposition for the scalar-plus-dust model; a finite-tension critical mass with its small-tension limit and a coincidence of horizon areas there, re-derived independently (Supplement S8); "
       "a regular-centre Einstein–scalar experiment that separates persistence from creation, with its dependence on wall thickness quantified (Section 5.6); and dimensioned numbers at the observed dark-energy scale, "
       "with the exterior mass parameter distinguished from the quasi-local mass (Section 4). The vacuum energy is the observed dark-energy density, so the critical sizes and masses \"\n"
       "          \"are numbers rather than free parameters, and they can be compared with the matter that the universe actually contains. The paper does not add a new bubble solution, a traversable throat or a global spacetime for the numerical solutions. It makes five contributions.\"))")
s = sub(s, old, new)
s = sub(s, 'add(("bullet", "An accounting of the size, mass and energy budget that the criterion requires at the observed dark-energy scale, showing that the critical mass is dominated by the vacuum-energy difference and not by matter (Section 7), and falsification criteria (Section 8)."))',
        'add(("bullet", "An accounting of the size, mass and energy budget that the criterion requires at the observed dark-energy scale, showing that the critical mass is dominated by the vacuum-energy difference and not by matter (Section 7), and falsification criteria (Section 8)."))\n'
        'add(("bullet", "Audits of the evidence: an independent derivation of the junction branches, horizons and masses (Supplement S8), a quantitative wall-thickness extrapolation (Section 5.6) and a causal-classification diagnostic that separates the evolved domain from conjecture (Section 5.4)."))')
# 2.1 scope
old = ("One limit of scope follows. The dark matter is not \"\n"
       "          \"simulated: the numerical experiments evolve the scalar field alone in a prescribed double-well potential, so they study the gravitational problem of a false-vacuum region bounded by a wall and not the \"\n"
       "          \"complete matter-triggered process. The dark-matter coupling enters only through the existence of the false vacuum, and the claim that the whole process has been simulated is not made.\"))")
new = ("One limit of scope is stated at once, because it is the largest gap between the motivating model and the computed one: the numerical experiments are a scalar-only gravitational benchmark. \"\n"
       "          \"They evolve the scalar field alone in a prescribed double-well potential, with the false vacuum held by the potential and not by dark matter, so they study the gravitational problem of a false-vacuum region bounded by a wall and not the \"\n"
       "          \"matter-triggered process. What depends on the dark matter and what does not is as follows. The null-energy proposition (Section 3) is proved for the scalar together with non-negative-density dust and does not use the numerical solutions. \"\n"
       "          \"The thin-wall criterion (Section 4) depends on the matter only through the exterior mass parameter, to which dust contributes its conserved mass, and through the vacuum energies. The numerical persistence and creation results (Section 5) are results for the scalar-only system: \"\n"
       "          \"whether positive-energy dust inside the core or the wall can change the classification of any run was not tested, and a controlled subset with dust, which needs a coupled solver, is the most consequential extension. The claim that the whole process has been simulated is not made.\"))")
s = sub(s, old, new)
# section 4 scope
s = sub(s, "outcome A for every regular configuration of subcritical mass, with no assumption on the initial velocity.",
        "outcome A for every configuration of subcritical mass on the branch analysed here (spherical, thin wall, vacuum exterior of the stated form, regular initial slice), with no assumption on the initial velocity.")
# mass definition paragraph at the end of section 4
anchor = 'add(("p", "Two statements follow. First, $m>m_c(\\\\hat\\\\sigma)$ is necessary for any regular configuration to reach $\\\\Gamma_{in}<0$.'
mass_par = '''add(("p", "Which mass is meant needs a statement. The $m$ of Eqs. (6)–(9) is the parameter $m=2GM/c^2$ of the exterior metric function $f=1-m/r-\\\\rho_{out}r^2/3$: for $q>0$ the Schwarzschild–anti-de Sitter mass, for $q=0$ the ADM mass, and for $q<0$ the "
          "Schwarzschild–de Sitter parameter, which is not an asymptotic mass because there is no asymptotically flat region. The quasi-local Misner–Sharp mass at the wall is $m_{\\\\rm MS}(R)=m+\\\\rho_{out}R^3/3$, so the two differ by the exterior vacuum energy. "
          "At $R=R_{dS}$ and vanishing tension the interior and exterior quasi-local masses coincide and equal $R_{dS}$ in these units for every $q$, whereas the parameter is $(1+q)R_{dS}$ (Supplementary Table S6). "
          "For $q=1$ the critical mass quoted in this paper, %s $M_\\\\odot$, is therefore the exterior parameter, twice the quasi-local mass at the wall (%s $M_\\\\odot$); statements about the mass of the region mean the parameter unless ‘quasi-local’ is written. "
          "An independent derivation (Supplement S8) reproduces Eqs. (7)–(9) from the extrinsic curvatures, the critical mass with 120-digit arithmetic (for the model's tension $m_c-2=%s$), the horizon roots at $m_c$ and the conversion to solar masses to $10^{-16}$." %
          (sci(AUD["5_dimensional_restoration"]["q=1"]["M_c_Msun"], 2), sci(AUD["5_dimensional_restoration"]["q=1"]["M_misner_sharp_at_wall_Msun"], 2), sci(float(AUD["2_physical"]["m_c_minus_limit"]), 2))))
'''
s = sub(s, anchor, mass_par + anchor)
wr("content_a.py", s)

# ------------------------------------------------------------------------------------------------ content_b
s = rd("content_b.py")
s = sub(s, "Only the scalar field is evolved; the dark matter of the full model is not included (Section 2.1).", "Only the scalar field is evolved: this is a scalar-only gravitational benchmark, and the dark matter of the full model is not included (Section 2.1).")
s = sub(s, "For $R_w=2R_{dS}$ (panel b) a region of positive expansions in both directions fills the core between the de Sitter horizon and the wall; beyond the wall there is first a ",
        "For $R_w=2R_{dS}$ (panel b) a region of positive expansions in both directions fills the core between the de Sitter horizon and the wall and, in the earliest rows only, extends a little beyond it (Fig. 5); beyond the wall there is first a ")
anchor = 'add(("p", "Two limits of this evidence should be explicit. First,'
c1 = CAUS = None
cause = '''_CE = CAUS["run_expanding_core"]
add(("p", "To separate what the numerical solutions establish from what is conjecture, a causal-classification diagnostic was run on the two runs of Fig. 4 (causal_diagnostic.py). In the evolved domain of the expanding-core run (%d rows, $u\\\\le%s$) the fractions of the lattice in the four sign quadrants are "
          "%s normal, %s trapped, %s anti-trapped and %s with $r_u>0>r_v$. The wall, located by the half-way value of $\\\\phi$ on every row, stays within %s of $v=V_w$, so that the core and the exterior are the same regions of the chart at all $u$. "
          "One statement follows from the chart and one from the numbers. From the chart: a future-directed causal curve has non-decreasing $u$ and $v$, so no exterior point lies in the causal past of a core point; a lattice test finds no exterior cell in the causal past of any of the %s anti-trapped cells inside the core, "
          "and %s of the core cells have an exterior cell in their causal future, so signals can leave the core and cannot enter it. From the numbers: the anti-trapped quadrant is not confined to the core. In the earliest rows, up to $u=%s$, %s of the %s anti-trapped cells (%s) lie beyond the wall, in the initial data of the exterior, "
          "before that region becomes trapped, and those cells do have exterior cells in their causal past. The statement that the exterior cannot send matter into the expanding region is therefore a statement about the false-vacuum core, and in these runs the exterior contains no matter. "
          "Fig. 5 marks what is established (the evolved domain, panel a) and what is not (panel b): the continuation beyond the last row, whether the anti-trapped region reaches a second asymptotic region or closes on itself, and any spacelike slice through a throat. None of these is computed, and the solutions are not labelled wormholes or bridges." %
          (_CE["rows_evolved"], n(_CE["u_last"], 2), pc(_CE["fractions"]["normal"], 0), pc(_CE["fractions"]["trapped"], 0), pc(_CE["fractions"]["anti_trapped"], 0), pc(_CE["fractions"]["r_u>0>r_v"], 0), n(_CE["wall_drift_max"], 2),
           n(_CE["anti_trapped_in_core"]), pc(_CE["fraction_of_core_cells_whose_future_reaches_exterior"], 0), n(_CE["u_last_row_with_anti_trapped_beyond_wall"], 2), n(_CE["anti_trapped_beyond_wall"]), n(_CE["anti_trapped_cells"]),
           pc(_CE["anti_trapped_beyond_wall"] / _CE["anti_trapped_cells"], 0))))
add(("fig", "Fig5_causal.png", "Fig. 5. Causal classification for $R_w=2R_{dS}$. (a) Sign quadrants of the null expansions on the evolved domain, with the wall (dashed), the last evolved row (dotted) and the region that is not computed. "
                               "(b) What is established (core and exterior fixed regions of the chart; no exterior-to-core causal curve) and what is conjecture and not computed (the global continuation).", 6.6))
'''
s = sub(s, anchor, cause + anchor)
# renumber old Fig. 5 -> Fig. 6
s = sub(s, '(Fig. 5, Supplementary Table S3)', '(Fig. 6, Supplementary Table S3)')
s = sub(s, 'add(("fig", "Fig5_other_vacua.png", "Fig. 5. Outcomes for exteriors', 'add(("fig", "Fig5_other_vacua.png", "Fig. 6. Outcomes for exteriors')
# new section 5.6 appended at the end
s = s.rstrip("\n") + '''

add(("h2", "5.6 Extrapolation to the physical wall"))
_BA, _BW, _SW, _SA = WALLX["by_amplitude"], WALLX["by_width"], WALLX["deviation_slope_vs_w"], WALLX["deviation_slope_vs_A"]
add(("p", "The experiments use walls that are thicker than the physical wall by a factor %s and heavier by %s (Table 1); resolving the physical wall on the grid used here would take about %s points per de Sitter radius. Which conclusions can be carried to the physical wall? Three classes can be separated, and the third cannot be answered by this paper." %
          (sci(WALLX["extrapolation_factors"]["thickness_ratio_sim_to_phys"], 1), sci(WALLX["extrapolation_factors"]["tension_ratio_sim_to_phys"], 1), sci(WALLX["extrapolation_factors"]["grid_points_needed_to_resolve_physical_wall"], 1))))
add(("p", "First, conclusions fixed for any thickness by argument. Proposition 1 uses only the null energy condition. The existence of anti-trapped points in the initial core depends on whether the core radius exceeds $R_{dS}$, with an offset of the order of the wall width (Section 5.4). "
          "The thin-wall mass depends on the tension only through $\\\\hat\\\\sigma^2$: at the simulated tensions ($\\\\hat\\\\sigma=%s$–$%s$) the critical mass is %s to %s below its limit, and at the physical tension the shift is $6\\\\times10^{-61}$ (Fig. 7b; Supplement S8). "
          "Evaluating the criterion at the finite simulated tension instead of the limit does not change the classification of any run: the smallest $m_{\\\\rm ext}/m_c$ among the expanding cores is %s with the finite-tension mass and %s with the limit." %
          (n(_BA["0.05"]["sigma_hat_sim"], 2), n(_BA["0.2"]["sigma_hat_sim"], 2), pc(abs(_BA["0.05"]["m_c_shift"]), 1), pc(abs(_BA["0.2"]["m_c_shift"]), 1),
           n(min(v["min_m_over_mc_expanding_core_finite_tension"] for v in _BA.values()), 2), n(min(v["min_m_over_mc_expanding_core_thin_limit"] for v in _BA.values()), 2))))
add(("p", "Second, conclusions that are insensitive within the resolution. The formation of a trapped surface at subcritical mass holds in every run. The radius of the first trapped point agrees with the horizon radius of the measured mass to a median of %s–%s (Section 5.3), against a resolution floor $h/r_h$ of about %s. "
          "The deviation shows no trend with the amplitude $A$, which scales thickness and tension together (95%% interval of the slope %s to %s per unit $A$, consistent with zero), but it does depend on the wall width: the median signed deviation is %s, %s and %s for $w=0.05$, 0.1 and 0.2, "
          "a slope of %s per unit $w/R_{dS}$ (95%% interval %s to %s). A linear extrapolation to $w\\\\to0$ gives %s (95%% interval %s to %s), consistent with the thin-wall value within the resolution floor (Fig. 7a). "
          "The thin-wall prediction is therefore recovered to about one per cent and no better, and the finite-width effect that is visible is of that order." %
          (pc(min(v["median_abs_dev"] for v in _BA.values()), 2), pc(max(v["median_abs_dev"] for v in _BA.values()), 2), pc(_SW["resolution_floor_median"], 1), n(_SA["ci95"][0], 2), n(_SA["ci95"][1], 2),
           pc(_BW["0.05"]["median_signed"], 2), pc(_BW["0.1"]["median_signed"], 2), pc(_BW["0.2"]["median_signed"], 2), n(_SW["slope_median"], 3), n(_SW["ci95"][0], 3), n(_SW["ci95"][1], 3),
           pc(_SW["intercept_at_w0_median"], 2), pc(_SW["intercept_ci95"][0], 2), pc(_SW["intercept_ci95"][1], 2))))
add(("p", "Third, what cannot be extrapolated: the dynamics of a wall $10^{-30}R_{dS}$ thick, which is neither resolved nor tested; any quantity at the level of one per cent or finer in the horizon radius; the stability of the wall and its energy exchange with the exterior; the fate of matter after the code stops at the singularity; "
          "and non-spherical or quantum wall dynamics. The simulations therefore do not validate the physical wall regime. They support the necessary-condition statement $m>m_c$, the null-energy result and the dependence of the expanding core on a core larger than its own de Sitter radius, "
          "and they show that the thin-wall structure is present at finite thickness to the stated accuracy."))
add(("fig", "Fig7_wall_extrapolation.png", "Fig. 7. (a) Radius of the first trapped point relative to the horizon radius of the measured exterior mass against the wall width, with a linear fit extrapolated to zero width and the resolution floor. "
                                           "(b) Relative shift of the critical mass with the wall tension (line: Eq. (9); points: the exact maximum at the three simulated tensions); the physical tension lies far to the left and its correction is $6\\\\times10^{-61}$.", 6.6))
'''
wr("content_b.py", s)

# ------------------------------------------------------------------------------------------------ content_c
s = rd("content_c.py")
s = sub(s, "(Fig. 6b)", "(Fig. 8b)")
s = sub(s, 'add(("fig", "Fig6_assembly.png", "Fig. 6. (a) Critical mass', 'add(("fig", "Fig6_assembly.png", "Fig. 8. (a) Critical mass')
s = sub(s, "A false-vacuum region smaller than $R_{dS}$, whatever its mass, is outcome A: this is what the numerical runs with $R_w\\\\lesssim R_{dS}$ show and, \"\n          \"for $m<m_c$, what the thin-wall analysis proves.",
        "In the scanned initial data a false-vacuum region smaller than $R_{dS}$ is outcome A whatever its mass, and for $m<m_c$ on the thin-wall branch of Section 4 \"\n          \"this is what the analysis proves.")
s = sub(s, "that is a mass $M_c$ = %s $M_\\\\odot$, for a region of radius $R_{dS}$ = %s Gpc. ", "that is an exterior mass parameter $M_c$ = %s $M_\\\\odot$ (the quasi-local mass at the wall is half of it for $q=1$, Section 4), for a region of radius $R_{dS}$ = %s Gpc. ")
# section 8
s = sub(s, 'add(("p", "Collecting the results, the model has three regimes for a false-vacuum region of core radius $R_w$ and exterior mass $m$, in units of $R_{dS}$ and with $q=1$:"))',
        'add(("p", "Collecting the results, the model has three regimes for a false-vacuum region of core radius $R_w$ and exterior mass $m$, in units of $R_{dS}$ and with $q=1$. They are limited to the thin-wall branch of Section 4 and to the initial data that were scanned (%d runs for $q=1$ and 14 for each other $q$); they are not theorems about all configurations:" % SC["n_runs"]))')
s = sub(s, "the wall is confined and collapses; a trapped surface forms (outcome A). This includes every region that is smaller than $R_{dS}$ and has no more than its own vacuum energy.\"))",
        "for regular initial data on the spherical thin-wall branch, the wall is confined and collapses; a trapped surface forms (outcome A). In the scan this holds for every run with $R_w\\\\lesssim R_{dS}$ and $m<m_c$; no theorem extends it beyond the branch and the data studied.\"))")
# section 9 first paragraph
s = sub(s, "below a finite-tension critical mass \"\n          \"every regular configuration collapses;", "on the thin-wall branch analysed, below a finite-tension critical mass \"\n          \"every regular configuration collapses;")
s = sub(s, "The dark matter that motivates the transition is not simulated: the experiments evolve the scalar field in a prescribed double-well potential, so the complete matter-triggered process has not been \"\n          \"simulated.",
        "The dark matter that motivates the transition is not simulated: the experiments are a scalar-only gravitational benchmark, with the scalar field in a prescribed double-well potential, so the complete matter-triggered process has not been \"\n          \"simulated and the effect of dust on the classification of the runs is untested.")
# 9.2 and renumber Table 5 -> Table 6 not needed (no new table in main text)
anchor = 'add(("h1", "10. Conclusions"))'
sec92 = '''add(("h2", "9.2 Externally seeded vacuum transitions and critical growth"))
_RC = 2.0 * SIGMA / (3.0 * 2.0)
add(("p", "A further speculative extension is whether an externally prepared, localized high-energy configuration could initiate a vacuum transition whose subsequent evolution is driven by the vacuum-energy difference rather than by continued external energy input. "
          "This is analogous, at the level of threshold dynamics, to a self-propagating phase transition, but not to a nuclear chain reaction in its microscopic mechanism. The initiating configuration could in principle arise through an engineered intervention by a sufficiently advanced technological civilization; "
          "no such intervention is assumed or modelled here."))
add(("p", "In the non-gravitating thin-wall approximation a bubble separating two vacuum phases has the energy"))
add(("eq", "E(R)=4\\\\pi\\\\sigma R^2-\\\\frac{4\\\\pi}{3}\\\\Delta V\\\\,R^3,\\\\qquad R_c=\\\\frac{2\\\\sigma}{\\\\Delta V}", "12"))
add(("p", "where $\\\\Delta V>0$ is the energy-density advantage of the interior phase. Above $R_c$ expansion is energetically favourable within the approximation. In the units of this paper $\\\\Delta V=(1+q)V_f=3(1+q)/R_{dS}^2$, so that $R_c/R_{dS}=2\\\\hat\\\\sigma/[3(1+q)]$, which for the model's tension and $q=1$ is %s, about %s mm: "
          "of the order of the wall's own thickness (%s mm), so the estimate sits at the edge of the validity of the thin-wall picture. Gravity can change the conclusion substantially, and the junction conditions, the causal structure, the null energy condition and the stability of the wall decide whether such growth is physically possible. "
          "A growing vacuum bubble is not necessarily a wormhole, a nonsingular daughter universe or a mechanism for restarting the parent cosmology." % (sci(_RC, 1), n(_RC * TW["R_dS_m"] * 1e3, 2), n(TW["wall_thickness_mm"], 2))))
add(("p", "The present results impose an important constraint on this proposal: localized collapse did not create an expanding daughter region in the numerical experiments, and the expanding cores that persisted were already present in the initial data. "
          "The gravitational threshold of Section 4 is also not lowered by the energy of a seed, because $m_c$ is fixed by the vacuum-energy difference at the scale $R_{dS}$. Consequently, an engineered seed would require a physically justified mechanism that overcomes the identified formation obstruction rather than merely increasing the initial energy."))
add(("p", "A future investigation could compare localized perturbations, externally prepared supercritical bubbles and initial data containing a pre-existing bridge. The decisive tests would be whether regular initial data evolve into a newly formed expanding region, whether the solution remains stable, and whether a future-directed causal continuation exists. "
          "Such a study would distinguish vacuum-transition propagation from genuine cosmological reproduction."))
'''
s = sub(s, anchor, sec92 + anchor)
# conclusions: mass parameter wording
s = sub(s, "At the observed dark-energy scale that core has a radius of %s Gpc and a vacuum-energy deficit of %s $M_\\\\odot$,", "At the observed dark-energy scale that core has a radius of %s Gpc and an exterior mass parameter, equal to the vacuum-energy deficit, of %s $M_\\\\odot$,")
# supplement: new sections appended
s = s.rstrip("\n") + '''

adds(("h1", "S8. Independent audit of the junction branches, critical mass, horizons and masses"))
adds(("p", "independent_audit.py re-derives the junction conditions symbolically from the extrinsic curvatures and does not use the closed forms of thinwall.py except for the final comparison. Eliminating $\\\\Gamma_{out}$ from $\\\\Gamma_{in}-\\\\Gamma_{out}=\\\\sigma R/2$ and $\\\\Gamma_{out}^2-\\\\Gamma_{in}^2=f_{out}-f_{in}$ gives Eq. (7), and the roots of $U^2=0$ in $m$ reproduce Eq. (8): "
           "the maximum deviation over 200 random $(R,\\\\sigma,q)$ is %s, and the two branches meet at $R=R_{dS}$ at $m=(\\\\Delta/3-\\\\sigma^2/4)R_{dS}^3$. The maximum of $m_+(R)$ was recomputed with 120-digit arithmetic (Table S5); at the model's tension $m_c-2=%s$, the series gives %s, and the maximum lies $%s$ below $R_{dS}$." %
           (n(AUD["1_max_deviation_symbolic_vs_closed_form"], 1), AUD["2_physical"]["m_c_minus_limit"][:9], AUD["2_physical"]["series_prediction"][:9], AUD["2_physical"]["one_minus_R_at_max"][:8])))
adds(("table", "Table S5. Critical mass recomputed with 120-digit arithmetic against the small-tension series of Eq. (9) ($R_{dS}=1$).",
      ["$q$", "$\\\\hat\\\\sigma$", "$m_c$ (120 digits)", "Series, Eq. (9)", "$R$ at the maximum"],
      [[n(r["q"], 1), n(r["sigma"], 2), n(r["m_c_hp"], 9), n(r["series"], 9), n(r["R_at_max"], 6)] for r in AUD["2_critical_mass_hp"]], [0.6, 0.8, 1.7, 1.7, 1.5]))
adds(("p", "Horizon roots at $m=m_c=(1+q)R_{dS}$ and mass definitions are in Table S6. For $q>0$ the single exterior horizon is at $R_{dS}$; for $q=0$ it is the Schwarzschild radius; for $q=-0.5$ the roots are the black-hole horizon at $0.618R_{dS}$ and the cosmological horizon at $R_{dS}$, so the equality of horizon areas is with the cosmological horizon there. "
           "The quasi-local Misner–Sharp mass at the wall at $R=R_{dS}$ equals $R_{dS}$ on both sides for every $q$, while the exterior parameter is $(1+q)R_{dS}$. Dimensional restoration used independently chosen constants ($c$, $G$, $M_\\\\odot$, Planck 2018 $H_0$ and $\\\\Omega_\\\\Lambda$): $R_{dS}=%s$ Gpc, and the critical masses agree with thinwall.json to %s." %
           (n(AUD["5_dimensional_restoration"]["R_dS_Gpc"], 4), sci(max(abs(v["rel_diff"]) for v in AUD["5_comparison"].values()) + 1e-16, 0))))
adds(("table", "Table S6. Horizons and masses at $m=m_c$ ($R_{dS}=1$, $8\\\\pi G=1$): exterior parameter, quasi-local mass at the wall at $R=R_{dS}$ and horizon roots; critical mass in solar masses with the quasi-local value.",
      ["$q$", "$m$ (parameter)", "$m_{\\\\rm MS}$ at the wall", "Horizon roots $r/R_{dS}$", "$M_c$ ($M_\\\\odot$)", "$M_{\\\\rm MS}$ ($M_\\\\odot$)"],
      [[n(h_["q"], 2), n(h_["m_c"], 2), n(m_["misner_sharp_exterior_at_wall"], 2), ", ".join(n(x, 3) for x in h_["roots"]), sci(AUD["5_dimensional_restoration"]["q=%g" % h_["q"]]["M_c_Msun"], 2) if ("q=%g" % h_["q"]) in AUD["5_dimensional_restoration"] else "–",
        sci(AUD["5_dimensional_restoration"]["q=%g" % h_["q"]]["M_misner_sharp_at_wall_Msun"], 2) if ("q=%g" % h_["q"]) in AUD["5_dimensional_restoration"] else "–"] for h_, m_ in zip(AUD["3_horizons_at_mc"], AUD["4_mass_definitions"])],
      [0.6, 1.1, 1.3, 1.6, 1.3, 1.3]))
adds(("h1", "S9. Wall-thickness extrapolation and causal diagnostic"))
adds(("p", "Table S7 gives, for each amplitude $A$ (thickness and tension scale together), the wall thickness and tension, the median and maximum absolute deviation of the first trapped radius from the horizon radius of the measured mass, the resolution floor, the shift of the critical mass at the simulated tension and the smallest $m_{\\\\rm ext}/m_c$ among the expanding cores with the limit and with the finite-tension mass (wall_extrapolation.py). "
           "Table S8 gives the causal diagnostic for the two runs of Fig. 5 (causal_diagnostic.py)."))
adds(("table", "Table S7. Wall-thickness dependence of the numerical outputs (normal-core runs with a trapped surface).",
      ["$A$", "Thickness/$R_{dS}$", "$\\\\hat\\\\sigma$", "median |dev|", "max |dev|", "floor $h/r_h$", "$m_c$ shift", "min $m/m_c$ (limit)", "min $m/m_c$ (finite $\\\\hat\\\\sigma$)"],
      [[k_, n(v["thickness_over_RdS"], 4), n(v["sigma_hat_sim"], 3), pc(v["median_abs_dev"], 2), pc(v["max_abs_dev"], 1), pc(v["resolution_floor_median"], 1), pc(v["m_c_shift"], 2), n(v["min_m_over_mc_expanding_core_thin_limit"], 2), n(v["min_m_over_mc_expanding_core_finite_tension"], 2)] for k_, v in WALLX["by_amplitude"].items()],
      [0.5, 1.0, 0.7, 0.9, 0.9, 0.9, 0.8, 1.1, 1.3], 14))
adds(("table", "Table S8. Causal diagnostic on the evolved domain.",
      ["Run", "Rows", "$u_{\\\\rm last}$", "Anti-trapped cells (core / wall / beyond wall)", "Core anti-trapped with exterior in causal past", "Wall drift"],
      [[r_["label"], n(r_["rows_evolved"]), n(r_["u_last"], 2), "%s / %s / %s" % (n(r_["anti_trapped_in_core"]), n(r_["anti_trapped_in_wall"]), n(r_["anti_trapped_beyond_wall"])), n(r_["core_anti_trapped_with_exterior_in_causal_past"]), n(r_["wall_drift_max"], 3)] for r_ in (CAUS["run_normal_core"], CAUS["run_expanding_core"])],
      [1.3, 0.6, 0.8, 2.2, 2.0, 0.8], 14))
'''
wr("content_c.py", s)
print("revision 2 edits applied")
