# -*- coding: utf-8 -*-
"""Section 5: the numerical Einstein-scalar experiment."""
from pcommon import *

TST = VAL["vacuum_tests"]
SC = VAL["scan"]
PER = VAL["per_A"]
CONV = VAL["trapped_convergence"]
TH = VAL["thresholds"]
cls = VAL["classes"]
n_AT = sum(1 for c in cls if c["T"] and c["AT"])
n_T = sum(1 for c in cls if c["T"] and not c["AT"])
n_inc = sum(1 for c in cls if not c["T"])
# thresholds of the anti-trapped core
no_AT_thin = max(t["R_w_no_AT"] for t in TH if t["w"] == 0.05)
AT_thin = min(t["R_w_AT"] for t in TH if t["w"] == 0.05)
no_AT_thick = max(t["R_w_no_AT"] for t in TH if t["w"] > 0.05)
AT_thick = min(t["R_w_AT"] for t in TH if t["w"] > 0.05)
rich = (CONV[0]["r"], CONV[1]["r"], CONV[2]["r"])
rext = rich[2] + (rich[2] - rich[1]) / 3.0

add(("h1", "5. Numerical experiment: collapse of a false-vacuum core with a regular centre"))
add(("h2", "5.1 Formulation"))
add(("p", "The thin-wall result is exact but assumes an infinitely thin, infinitely light wall. To test its structure with a wall of finite thickness and tension, the coupled Einstein–scalar equations are evolved in "
          "double-null gauge with a regular centre, the standard formulation of spherical scalar-field collapse [@Christodoulou1986; @Choptuik1993], which resolves trapped surfaces and expanding regions and requires no choice of slicing near horizons. The equations are Eqs. (4) and (5) with the scalar equation "
          "and the equation for $f$,"))
add(("eq", "\\phi_{uv}=-\\frac{\\phi_ur_v+\\phi_vr_u}{r}-\\frac{f}{4}U'(\\phi),\\qquad (\\ln f)_{uv}=\\frac{f+4r_ur_v}{2r^2}-\\phi_u\\phi_v", "10"))
add(("p", "with $U(\\phi)=(\\kappa/4)(\\phi^2/A^2-1)^2+s\\phi/A+d$, tuned so that the false minimum has density $\\rho_{in}=3$ (so that $R_{dS}=1$) and the true minimum has density $-q\\rho_{in}$, with $q=1$ and "
          "$\\kappa=24$. The field amplitude $A$ is the only free scale of the wall: its thickness is $\\ell=A\\sqrt{2/\\kappa}$ and its tension is of the order of $A\\sqrt{\\kappa/2}$, so that $A\\to0$ is "
          "the thin, light wall limit. The scheme is the second-order diamond predictor–corrector on a uniform lattice [@Garfinkle1995; @Gundlach2007]: $r$ and $\\phi$ are updated cell by cell, the metric function "
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
          "at every point. Table 2 gives the largest relative error for three grid steps; the errors fall at second order. The radius of the first trapped point of a collapsing run ($R_w=1$, $w=0.1$, $A=0.1$) converges at the same "
          "order, %s, %s and %s at $h=0.02$, $0.01$ and $0.005$, with a Richardson extrapolation of %s." %
          (n(rich[0], 4), n(rich[1], 4), n(rich[2], 4), n(rext, 4))))
add(("table", "Table 2. Verification against exact vacuum solutions: largest relative error of the Misner–Sharp mass with a regular centre ($u\\le0.6$, $v\\le3$).",
     ["Vacuum", "$h=0.04$", "$h=0.02$", "$h=0.01$", "Orders"],
     [["de Sitter"] + [sci(x, 1) for x in TST["dS"]["max_rel_error"]] + [", ".join(n(o, 1) for o in TST["dS"]["orders"])],
      ["anti-de Sitter"] + [sci(x, 1) for x in TST["AdS"]["max_rel_error"]] + [", ".join(n(o, 1) for o in TST["AdS"]["orders"])]],
     [1.6, 1.1, 1.1, 1.1, 1.2]))
add(("h2", "5.3 Horizon formation"))
add(("p", "A scan of %d runs varies the field amplitude ($A=0.05$, 0.1 and 0.2), the wall width ($w=0.05$, 0.1 and 0.2) and the core radius $R_w/R_{dS}$ from 0.4 to 3 (Supplementary Table S1). Each run measures the "
          "exterior mass $m_{\\rm ext}$ on the initial cone beyond the wall, from the Misner–Sharp mass (constant across the exterior to within 2%% in 70 of the 72 runs, at most 9.5%%), and the radius of the first trapped point. In %d runs a trapped surface forms before the code "
          "breaks down, and in %d runs (low-mass wide walls and two high-mass wide walls) the run ended before it could be classified; none of the runs reached the end of its domain without a trapped surface. The radius of the "
          "first trapped point agrees with the horizon radius of the Schwarzschild–anti-de Sitter solution of the measured mass, $r_h^3+r_h=m_{\\rm ext}$ for $q=1$, to a median of %s, a 95th percentile of %s and a maximum of %s (Fig. 2); the agreement is best "
          "at $A=0.1$ (median %s) and is not monotonic in $A$, which is expected because the apparent horizon forms in a region where the wall has not yet fully passed." %
          (SC["n_runs"], SC["n_trapped"], SC["n_no_trapped"], pc(abs(SC["dev_median"]), 1), pc(SC["dev_p95"], 1), pc(SC["dev_max_abs"], 1), pc(PER["0.1"]["median_abs"], 1))))
add(("fig", "Fig2_horizon_check.png", "Fig. 2. (a) Radius of the first trapped point against the horizon radius of the Schwarzschild–anti-de Sitter solution of the mass measured outside the wall, for all runs with a trapped surface; "
                                      "the line is equality. (b) Convergence of the first trapped point with the grid step.", 6.5))
add(("h2", "5.4 Expanding core"))
add(("p", "The question that decides between outcomes A and C is whether an anti-trapped region survives next to the trapped surface. Fig. 3 classifies the runs by core radius and exterior mass. A core that stays "
          "inside its own de Sitter radius ($R_w\\le%s\\,R_{dS}$ for the thin walls, $w=0.05$, and $R_w\\le%s\\,R_{dS}$ for the thicker ones, the offset being of the order of the wall width) is normal everywhere inside the wall ($%d$ runs), and the collapse is outcome A "
          "with a trapped region and a shrinking regular core. A core that extends beyond its own de Sitter radius ($R_w\\ge%s\\,R_{dS}$ for thin walls, $\\ge%s\\,R_{dS}$ for thick walls; the offset is of order the wall width) "
          "contains a region with $r_u>0$ and $r_v>0$ from the start, which persists to the end of every such run, next to a trapped region: $%d$ runs. The thin-wall criterion is a necessary condition in every run: all "
          "of the runs with an expanding core have $m_{\\rm ext}>m_c$ (Fig. 3, filled squares above the dotted line). Mass alone does not decide the outcome, however: several runs with $m_{\\rm ext}$ up to %s times $m_c$ "
          "and a core inside the de Sitter radius have no expanding region." %
          (n(no_AT_thin, 1), n(no_AT_thick, 1), n_T, n(AT_thin, 1), n(AT_thick, 1), n_AT, n(max(c["m_over_mc"] for c in cls if c["T"] and not c["AT"]), 1))))
add(("fig", "Fig3_outcome_map.png", "Fig. 3. Outcome of each run against core radius and exterior mass in units of the thin-wall critical mass. Open circles: a trapped surface forms and the core is normal "
                                    "(outcome A); filled squares: a trapped surface forms and an expanding core persists; crosses: the run ended before classification. The dashed line is the core radius equal to $R_{dS}$ and the dotted line is $m=m_c$. "
                                    "Colours distinguish $A=0.05$ (blue), 0.1 (red) and 0.2 (green).", 5.2))
add(("p", "Fig. 4 shows the causal structure of two runs in terms of the signs of the null expansions. For $R_w=0.8R_{dS}$ (panel a) the region inside the ingoing wall stays normal and the region outside it becomes trapped "
          "behind the apparent horizon at $r\\approx0.79$. For $R_w=2R_{dS}$ (panel b) a region of positive expansions in both directions fills the core between the dS horizon and the wall; beyond the wall there is first a "
          "region with $r_u>0>r_v$, the second-exterior quadrant of a bridge, and then a trapped region. The geometry on the two sides of the wall is that of the neck of Section 4: the areal radius has a minimum "
          "on a slice that crosses the wall, between an expanding interior and an exterior that is trapped. No null ray has a minimum of the areal radius, as Proposition 1 requires: "
          "over all %d rows of all %d runs, no outgoing ray has a point at which $r_v$ changes from negative to positive." % (VAL["null_ray_minima"]["rows"], VAL["null_ray_minima"]["runs"])))
add(("fig", "Fig4_sign_maps.png", "Fig. 4. Quadrants of the null expansions in the $(u,v)$ plane for two runs ($A=0.1$, $w=0.1$): (a) a core inside its de Sitter radius, (b) a core outside it. The dashed line is the wall. "
                                  "Grey: normal ($r_u<0<r_v$); red: trapped; blue: anti-trapped (expanding); orange: $r_u>0>r_v$. The runs stop when the code can no longer follow the trapped region.", 6.5))
add(("p", "Two limits of this evidence should be explicit. First, it establishes that the structure predicted by the thin-wall theory is present at finite thickness and tension; it does not establish what lies beyond "
          "the horizon after the code ends, so the statement that an expanding core is a daughter region rests on the persistence of the anti-trapped region, which the dynamics cannot undo, and not on a "
          "global evolution. Second, within the scan the thin-wall limit is approached but not reached; the thickness of the wall ranges from 1.4% to 5.8% of $R_{dS}$ (Table 1 gives $10^{-30}$ for the physical wall), so "
          "that the match of the apparent-horizon radius to the analytic value to a fraction of a percent is a test of the code and of the exterior solution, not of the thin-wall limit itself."))
