# -*- coding: utf-8 -*-
"""Section 5: the numerical Einstein-scalar experiment (revision 1)."""
from pcommon import *

TST = VAL["vacuum_tests"]
SC = VAL["scan"]
PER = VAL["per_A"]
CONV = VAL["trapped_convergence"]
TH = VAL["thresholds"]
cls = VAL["classes"]
CON, CONA = VE["constraints"], VE["constraints_AT"]
CO, COA = VE["constraint_orders"], VE["constraint_orders_AT"]
EXP = VE["expansions"]
INC = VE["inconclusive"]
CRE = VE["creation"]
OV = VE["other_vacua"]
n_AT = sum(1 for c in cls if c["T"] and c["AT"])
n_T = sum(1 for c in cls if c["T"] and not c["AT"])
resolved = [r for r in INC if any(v["first_trapped"] is not None for k, v in r.items() if k.startswith("refined"))]
unresolved = [r for r in INC if r not in resolved]
N_TRAP = SC["n_trapped"] + len(resolved)
no_AT_thin = max(t["R_w_no_AT"] for t in TH if t["w"] == 0.05)
AT_thin = min(t["R_w_AT"] for t in TH if t["w"] == 0.05)
no_AT_thick = max(t["R_w_no_AT"] for t in TH if t["w"] > 0.05)
AT_thick = min(t["R_w_AT"] for t in TH if t["w"] > 0.05)
rich = (CONV[0]["r"], CONV[1]["r"], CONV[2]["r"])
rext = rich[2] + (rich[2] - rich[1]) / 3.0
ratio_ok = lambda e: "%s, %s and %s" % (n(e["r"]["ratio"], 1), n(e["r_u"]["ratio"], 1), n(e["r_v"]["ratio"], 1))

add(("h1", "5. Numerical experiment: collapse of a false-vacuum core with a regular centre"))
add(("h2", "5.1 Formulation"))
add(("p", "The thin-wall result assumes an infinitely thin, infinitely light wall. To test its structure with a wall of finite thickness and tension, the coupled Einstein–scalar equations are evolved in "
          "double-null gauge with a regular centre, the standard formulation of spherical scalar-field collapse [@Christodoulou1986; @Choptuik1993], which resolves trapped surfaces and expanding regions and requires no choice of slicing near horizons. "
          "The equations are Eqs. (4) and (5) with the scalar equation and the equation for $f$,"))
add(("eq", "\\phi_{uv}=-\\frac{\\phi_ur_v+\\phi_vr_u}{r}-\\frac{f}{4}U'(\\phi),\\qquad (\\ln f)_{uv}=\\frac{f+4r_ur_v}{2r^2}-\\phi_u\\phi_v", "10"))
add(("p", "with $U(\\phi)=(\\kappa/4)(\\phi^2/A^2-1)^2+s\\phi/A+d$, tuned so that the false minimum has density $\\rho_{in}=3$ (so that $R_{dS}=1$) and the true minimum has density $-q\\rho_{in}$, with $\\kappa=24$ and, in the baseline, $q=1$. "
          "The field amplitude $A$ is the only free scale of the wall: its thickness is $\\ell=A\\sqrt{2/\\kappa}$ and its tension is of the order of $A\\sqrt{\\kappa/2}$, so that $A\\to0$ is "
          "the thin, light wall limit. Only the scalar field is evolved; the dark matter of the full model is not included (Section 2.1). The scheme is the second-order diamond predictor–corrector on a uniform lattice [@Garfinkle1995; @Gundlach2007]: $r$ and $\\phi$ are updated cell by cell, the metric function "
          "$\\sigma=\\ln f$ is integrated from the constraint within six points of the centre, where the evolution equation divides by a small radius, and from its evolution equation beyond; the centre value of "
          "$\\sigma$ follows from elementary flatness, $r_v=\\sqrt f/2$, solved by the secant method because the plain iteration has a gain of modulus greater than one, and the centre value of $\\phi$ from regularity, "
          "$\\phi_t=2\\phi_v$, by an implicit trapezoid. Data are given on the outgoing null cone from the centre in the gauge $r=v/2$: $\\phi(0,v)$ is false vacuum for $v<V_w$ and true vacuum beyond, with a "
          "$\\tanh$ transition of width $w$, and $\\sigma(0,v)$ follows from the constraint. The wall is then an ingoing wall of radius $R_w=V_w/2$, with the false-vacuum core inside it and no matter outside. "
          "The code is dnull.py; the equations were derived symbolically (derive_equations.py)."))
add(("p", "Three limits of the method should be stated at once. The scheme does not follow the evolution beyond the formation of a singularity, where $f$ diverges: a run ends when $\\sigma$ exceeds 12 or a value becomes "
          "non-finite, and the classification below uses only features that cannot be undone by the evolution that follows, a trapped surface and an anti-trapped region. The wall is far thicker and "
          "heavier than the physical wall by thirty and sixty orders of magnitude (Table 1), so the experiments test structure, not the values of the physical model. And the evolution is spherically symmetric "
          "and classical."))
add(("h2", "5.2 Verification"))
add(("p", "The solver was checked against the two exact vacuum solutions with a regular centre. For de Sitter space ($\\rho=3$) and anti-de Sitter space ($\\rho=-3$) the Misner–Sharp mass of Eq. (3) must equal $\\rho r^3/3$ "
          "at every point. Table 3 gives the largest relative error for three grid steps; the errors fall at second order. The radius of the first trapped point of a collapsing run ($R_w=1$, $w=0.1$, $A=0.1$) converges at the same "
          "order, %s, %s and %s at $h=0.02$, $0.01$ and $0.005$, with a Richardson extrapolation of %s." %
          (n(rich[0], 4), n(rich[1], 4), n(rich[2], 4), n(rext, 4))))
add(("table", "Table 3. Verification against exact vacuum solutions: largest relative error of the Misner–Sharp mass with a regular centre ($u\\le0.6$, $v\\le3$).",
     ["Vacuum", "$h=0.04$", "$h=0.02$", "$h=0.01$", "Orders"],
     [["de Sitter"] + [sci(x, 1) for x in TST["dS"]["max_rel_error"]] + [", ".join(n(o, 1) for o in TST["dS"]["orders"])],
      ["anti-de Sitter"] + [sci(x, 1) for x in TST["AdS"]["max_rel_error"]] + [", ".join(n(o, 1) for o in TST["AdS"]["orders"])]],
     [1.6, 1.1, 1.1, 1.1, 1.2]))
add(("p", "Exact vacuum solutions do not test the dynamics of a scalar wall, so the constraints were evaluated on the dynamical solutions. The constraint of Eq. (4) along $v$ and its analogue along $u$ are imposed only "
          "within six points of the centre, and elsewhere the evolved fields are free to violate them. Their root-mean-square residuals over the region $0.05\\le u\\le0.28$, $u+0.4\\le v\\le3$ of the run with a core inside its de Sitter "
          "radius ($R_w=1$) fall with the grid step at orders %s and %s for $C_{vv}$ and %s and %s for $C_{uu}$ between $h=0.02$, $0.01$ and $0.005$, and for the run with a core outside it ($R_w=2$) at orders %s and %s for $C_{vv}$ and %s and %s "
          "for $C_{uu}$ (Fig. S1). The residual of $C_{vv}$ in the first case is larger (RMS %s at $h=0.02$) because the region includes the wall, whose profile is resolved by only a few points at that step. "
          "The null expansions converge as well: the maximum difference of $r$, $r_u$ and $r_v$ at $u=0.2$ between successive resolutions falls with Cauchy ratios %s for the core inside its de Sitter radius (the ratio of $r_v$ is lower, "
          "about order 1.3, because of the sharp wall) and %s for the core outside it, at $u=0.12$ and with the same signs of $r_u$ and $r_v$ at all points." %
          (n(CO["vv"][0], 1), n(CO["vv"][1], 1), n(CO["uu"][0], 1), n(CO["uu"][1], 1), n(COA["vv"][0], 1), n(COA["vv"][1], 1), n(COA["uu"][0], 1), n(COA["uu"][1], 1),
           n(CON[0]["rms_vv"], 2), ratio_ok(EXP["normal_core"]), ratio_ok(EXP["expanding_core"]))))
add(("h2", "5.3 Horizon formation"))
add(("p", "A scan of %d runs varies the field amplitude ($A=0.05$, 0.1 and 0.2), the wall width ($w=0.05$, 0.1 and 0.2) and the core radius $R_w/R_{dS}$ from 0.4 to 3 (Supplementary Table S1). Each run measures the "
          "exterior mass $m_{\\rm ext}$ on the initial cone beyond the wall, from the Misner–Sharp mass (constant across the exterior to within 2%% in 70 of the 72 runs, at most 9.5%%), and the radius of the first trapped point. "
          "At the base resolution a trapped surface forms before the code breaks down in %d runs, and %d runs ended without one. These seven were examined individually (Table S2): the code stopped on non-finite values with $\\sigma$ still "
          "between 1.5 and 2.8, far from the divergence that marks a singularity, so that the stops are numerical. Refinement by factors of 2 and 4 gave a trapped surface in %d of them, so that trapped surfaces form in %d of the %d runs; "
          "the remaining run ($A=0.1$, $w=0.1$, $R_w=0.4R_{dS}$, a mass whose horizon radius, 0.2, is about two wall widths) stays unresolved at the finest step and is reported as such. "
          "The radius of the first trapped point agrees with the horizon radius of the Schwarzschild–anti-de Sitter solution of the measured mass, $r_h^3+r_h=m_{\\rm ext}$ for $q=1$, to a median of %s, a 95th percentile of %s and a maximum of %s "
          "over the %d base-resolution runs (Fig. 2). The agreement is best at $A=0.1$ (median %s) and is not monotonic in $A$, which is expected because the apparent horizon forms in a region where the wall has not yet fully passed; "
          "for the wide high-mass walls resolved only after refinement, the first trapped point lies 19%% to 29%% inside the final horizon radius for the same reason. This agreement checks the identification of the horizon "
          "against the expected exterior solution. It does not validate the daughter-universe interpretation, which rests on the quadrant structure and on the persistence of the expanding region discussed next." %
          (SC["n_runs"], SC["n_trapped"], SC["n_no_trapped"], len(resolved), N_TRAP, SC["n_runs"], pc(abs(SC["dev_median"]), 1), pc(SC["dev_p95"], 1), pc(SC["dev_max_abs"], 1), SC["n_trapped"], pc(PER["0.1"]["median_abs"], 1))))
add(("fig", "Fig2_horizon_check.png", "Fig. 2. (a) Radius of the first trapped point against the horizon radius of the Schwarzschild–anti-de Sitter solution of the mass measured outside the wall, for all runs with a trapped surface at the base resolution; "
                                      "the line is equality. (b) Convergence of the first trapped point with the grid step.", 6.5))
add(("h2", "5.4 Expanding core: persistence or creation"))
add(("p", "The question that decides between outcomes A and C is whether an anti-trapped region survives next to the trapped surface, and the question that decides between persistence and creation is whether it can appear where the initial "
          "data have none. Fig. 3 classifies the runs by core radius and exterior mass. A core that stays inside its own de Sitter radius ($R_w\\le%s\\,R_{dS}$ for the thin walls, $w=0.05$, and $R_w\\le%s\\,R_{dS}$ for the thicker ones, the offset being of the order "
          "of the wall width) is normal everywhere inside the wall in the initial data, and the collapse is outcome A with a trapped region and a shrinking regular core ($%d$ runs with a trapped surface at the base resolution). A core that extends "
          "beyond its own de Sitter radius ($R_w\\ge%s\\,R_{dS}$ for thin walls, $\\ge%s\\,R_{dS}$ for thick walls) contains a region with $r_u>0$ and $r_v>0$ in the initial data, which persists to the end of every such run, next to a trapped region: $%d$ runs. "
          "The thin-wall criterion is a necessary condition in every run: all of the runs with an expanding core have $m_{\\rm ext}>m_c$. Mass alone does not decide the outcome, however: several runs with $m_{\\rm ext}$ up to %s times $m_c$ "
          "and a core inside the de Sitter radius have no expanding region." %
          (n(no_AT_thin, 1), n(no_AT_thick, 1), n_T, n(AT_thin, 1), n(AT_thick, 1), n_AT, n(max(c["m_over_mc"] for c in cls if c["T"] and not c["AT"]), 1))))
add(("p", "The distinction between persistence and creation can be tested directly, and it is the main limitation of the evidence. Of the %d runs, %d start with a core that is normal everywhere in the initial data and %d start with an "
          "expanding core. In none of the %d runs that start normal does an anti-trapped point appear at any later row. The simulations therefore show the persistence of an expanding region that the initial data already contain, "
          "beside a trapped region that forms around it; they do not show the dynamical creation of a daughter universe from initially non-expanding, localized data, and for such data the evidence is the opposite." %
          (CRE["n_runs"], CRE["n_initially_normal_core"], CRE["n_initially_expanding_core"], CRE["n_initially_normal_core"])))
add(("fig", "Fig3_outcome_map.png", "Fig. 3. Outcome of each run against core radius and exterior mass in units of the thin-wall critical mass. Open circles: a trapped surface forms and the core is normal "
                                    "(outcome A); filled squares: a trapped surface forms and an expanding core, present in the initial data, persists; open triangles: a trapped surface forms only after refinement; the cross: unresolved. The dashed line is the core radius equal to $R_{dS}$ and the dotted line is $m=m_c$. "
                                    "Colours distinguish $A=0.05$ (blue), 0.1 (red) and 0.2 (green).", 5.2))
add(("p", "Fig. 4 shows the causal structure of two runs in terms of the signs of the null expansions. For $R_w=0.8R_{dS}$ (panel a) the region inside the ingoing wall stays normal and the region outside it becomes trapped "
          "behind the apparent horizon at $r\\approx0.79$. For $R_w=2R_{dS}$ (panel b) a region of positive expansions in both directions fills the core between the de Sitter horizon and the wall; beyond the wall there is first a "
          "region with $r_u>0>r_v$, the second-exterior quadrant of a bridge, and then a trapped region. No null ray has a minimum of the areal radius, as Proposition 1 requires: "
          "over all %d rows of all %d runs, no outgoing ray has a point at which $r_v$ changes from negative to positive. The runs do not construct a spacelike slice through the bridge, so no throat is claimed "
          "from them; the interpretation of the thin-wall bridge given in Section 4 rests on the exact solution." % (VAL["null_ray_minima"]["rows"], VAL["null_ray_minima"]["runs"])))
add(("fig", "Fig4_sign_maps.png", "Fig. 4. Quadrants of the null expansions in the $(u,v)$ plane for two runs ($A=0.1$, $w=0.1$): (a) a core inside its de Sitter radius, (b) a core outside it. The dashed line is the wall. "
                                  "Grey: normal ($r_u<0<r_v$); red: trapped; blue: anti-trapped (expanding); orange: $r_u>0>r_v$. The runs stop when the code can no longer follow the trapped region.", 6.5))
add(("p", "Two limits of this evidence should be explicit. First, it establishes that the structure predicted by the thin-wall theory is present at finite thickness and tension; it does not establish what lies beyond "
          "the horizon after the code ends, so the statement that an expanding core persists rests on the anti-trapped region, which the dynamics cannot undo, and not on a "
          "global evolution. Second, within the scan the thin-wall limit is approached but not reached; the thickness scale of the wall ranges from 1.4% to 5.8% of $R_{dS}$ (Table 1 gives $10^{-30}$ for the physical wall), so "
          "that the match of the apparent-horizon radius to the analytic value to a fraction of a percent is a test of the code and of the exterior solution, not of the thin-wall limit itself."))
add(("h2", "5.5 Exteriors of other vacuum depth"))
add(("p", "The baseline exterior is anti-de Sitter ($q=1$). To test how much of the picture depends on that choice, the scan was repeated with $A=0.1$, two wall widths and seven core radii for $q=2$ (a deeper anti-de Sitter exterior), "
          "$q=0$ (a Minkowski exterior) and $q=-0.5$ (a de Sitter exterior of energy $0.5V_f$), with the thin-wall threshold $m_c=(1+q)R_{dS}$ in each case (Fig. 5, Supplementary Table S3). For $q=2$ the picture of the baseline is reproduced: "
          "a trapped surface forms in %d of the %d runs, with a radius within a median of %s (maximum %s) of the horizon radius of the exterior, and an expanding core is present in the %d runs whose core extends beyond its de Sitter radius, none of them with a "
          "mass below $m_c$. For $q=0$ trapped surfaces were found in %d of the %d runs (the others ended on non-finite values before one formed), within a median of %s of the Schwarzschild radius $r_h=m_{\\rm ext}$, and again an expanding core "
          "in the %d runs with a core beyond $R_{dS}$ and none below $m_c$. For $q=-0.5$ the exterior has a black-hole horizon only for masses below the Nariai value $m=%s\\,R_{dS}$, and %d of the %d runs have a larger mass; no trapped surface was found "
          "in any of them, and an expanding core is present in the %d runs with a core beyond $R_{dS}$, again with $m_{\\rm ext}>m_c$ in all." %
          (OV["2.0"]["n_trapped"], OV["2.0"]["n"], pc(OV["2.0"]["dev_median"], 1), pc(OV["2.0"]["dev_max"], 1), OV["2.0"]["n_expanding_core"], OV["0.0"]["n_trapped"], OV["0.0"]["n"], pc(OV["0.0"]["dev_median"], 1),
           OV["0.0"]["n_expanding_core"], n(0.5443, 3), OV["-0.5"]["n_no_bh_horizon"], OV["-0.5"]["n"], OV["-0.5"]["n_expanding_core"])))
add(("p", "Two conclusions follow, and a third limits them. The necessary condition $m>m_c$ and the dependence of the expanding core on a core larger than its de Sitter radius do not depend on the sign of the exterior vacuum energy in these runs. "
          "What does depend on it is the horizon structure: for $q\\ge0$ the false-vacuum region is hidden behind a black-hole horizon, whereas for $q<0$ and large mass the exterior has no black-hole horizon at all, so that "
          "the configuration is not a daughter behind a horizon but an expanding region bounded by the exterior's cosmological horizon, and outcome C in the sense of this paper does not apply. The third is the limit: "
          "the $q=0$ and $q=-0.5$ scans are largely unresolved with respect to trapped surfaces, so that they support the statement about the expanding core and the threshold, not a full classification."))
add(("fig", "Fig5_other_vacua.png", "Fig. 5. Outcomes for exteriors of other vacuum depth ($A=0.1$; blue $w=0.05$, red $w=0.1$): core radius against exterior mass in units of $m_c=(1+q)R_{dS}$. The dashed line is $R_w=R_{dS}$ and the dotted line is $m=m_c$.", 6.5))
