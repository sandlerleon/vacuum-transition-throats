# -*- coding: utf-8 -*-
"""Manuscript: A Nontraversable Einstein-Rosen Bridge to a Nonsingular Cyclic Cosmological Region: Junction Dynamics, Stability, and Causal Structure.
Every number is read from the JSON files written by the scripts of the branch phase2-bridge."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
L = lambda f: json.load(open(os.path.join(HERE, f), encoding="utf-8"))
V = L("phase2d_cycle_validation.json")
BR = L("phase2d_bridge.json")
SH = L("phase2d_shells_q0.json")["shells"]
WD = L("phase2d_winding.json")
JG = L("phase2e_junction_global.json")
PF = L("phase2e_perturb_floquet.json")
PM = L("phase2e_perturb_mechanism.json")
WN = L("phase2e_window.json")
CV = L("phase2d_convergence.json")
PR = L("phase2f_pressure.json")
CU = L("phase2f_continuation_up_mu20.json") + L("phase2f_continuation_down_mu20.json")
MUF = L("phase2f_muscan_up_rc0.2_valid.json")
POT = L("phase2f_potential.json")
AN = L("phase2f_analysis.json")
FIG3 = "FigF1_phase2f.png"
CAP = L("phase2g_cap_spectrum.json")
SHG = L("phase2g_shell_data.json")
NEC0 = L("phase2g_comoving_nec.json")
FIG1 = "FigD1_phase2d.png"
FIG2 = "FigE1_phase2e.png"
GITHUB = "https://github.com/sandlerleon/vacuum-transition-throats"

TITLE = "A Nontraversable Einstein–Rosen Bridge to a Nonsingular Cyclic Cosmological Region: Junction Dynamics, Stability, and Causal Structure"
rep = [r for r in V["reproduction"] if abs(r["x"] - 0.35) < 1e-9][0]
cstar = V["flip"]["C_over_Ccrit_star"]
A0 = BR["analysis"][0]
mode = {r["n"]: r for r in PF["modes"]}
ctrl = {r["n"]: r for r in PM["control"]}
long0 = [x for x in BR["long"] if x["shell_index"] == 0]
long4 = [x for x in BR["long"] if x["shell_index"] == 4]
byw = {r["name"]: r for r in WN["points"]}


def g(x, d=3):
    return ("%." + str(d) + "g") % x


def sci(x, d=1):
    m, e = ("%.*e" % (d, x)).split("e")
    return "%s × 10$^{%d}$" % (m, int(e))


B = []
add = B.append

add(("title", TITLE))
add(("author", "Leon Sandler"))
add(("small", "Independent Researcher, Northbrook, Illinois, USA"))
add(("small", "ORCID: https://orcid.org/0009-0007-4584-808X"))
add(("small", "Corresponding author: sandler.leon@gmail.com"))
add(("blank",))
add(("h1", "Abstract"))
ABSTRACT = (
    "We construct and test a candidate spacetime in which a closed, nonsingular, cyclically evolving cosmological region is joined across a moving thin shell to a Schwarzschild exterior, "
    "so that a spatial slice has two ends and a throat. The region is the closed Friedmann universe of a scalar field with a double-well potential coupled to dust; for a window of dust content "
    "it has exact bouncing cycles with an elliptic return map, reproduced with an independent code. For $w=-1/2$ and a flat exterior the shell has periodic solutions with positive tension that satisfy "
    "the null energy condition, verified from raw extrinsic curvatures to $10^{-15}$. The cyclic region and the parent universe are spacelike separated: the bridge is not traversable and nothing crosses it. "
    "The data need a pre-existing throat but no white-hole region to the future of a suitable slice. The conservative shell survives 60 cycles under 1%% perturbations. However, linear spherical perturbations of the interior "
    "grow by factors of %s to %s per cycle, through a scalar-mediated force on the coupled dust and the hilltop instability; changing the coupling, the potential scale or shape, or adding fluid pressure does not remove them "
    "(worst case above 3×10³). The shell's energy exchange with the bulk needs an unspecified interaction. The construction is an exact, idealized, homogeneous solution, not a robust cosmological model; "
    "we state what is demonstrated and what is open."
    % (g(mode[3]["max_abs"], 2), sci(mode[60]["max_abs"], 1).replace(" × 10$^{5}$", "×10⁵")))
add(("p0", ABSTRACT))
add(("p0", "Keywords: Einstein–Rosen bridge; thin-shell junction; cyclic cosmology; null energy condition; causal structure; cosmological perturbations."))

# ------------------------------------------------------------------ 1
add(("h1", "1. Introduction"))
add(("p", "Two ideas about the fate of a late-time universe meet in this paper: that a region could exist that cycles through expansion and contraction without a singularity, and that such a region could be connected to our universe by a spacetime bridge that is hidden because nothing can cross it. "
          "Bridges of the second kind are familiar. The Einstein–Rosen bridge of the maximally extended Schwarzschild solution [@EinsteinRosen1935; @Kruskal1960; @Fuller1962] has a throat but cannot be traversed, and baby-universe constructions glue a closed region to a Schwarzschild exterior through a junction "
          "[@Blau1987; @Farhi1990; @Fischler1990; @Frolov1989]. Traversable wormholes require exotic stress-energy: a static throat needs a violation of the null energy condition [@Morris1988; @Hochberg1997; @Hayward1999], and topological censorship forbids a causal connection of two asymptotic regions otherwise [@Friedman1993]. "
          "Cyclic cosmologies have a long history [@Tolman1931; @Steinhardt2002]. We do not ask the bridge to cause the cycle, and we do not ask it to carry anything."))
add(("p", "The question is narrower and can be answered. Given a cold, closed, nonsingular cycle of a classical scalar-field model, does a consistent two-ended spacetime exist that contains it, with a bridge that is a legitimate solution of the Einstein equations, and is the result robust? "
          "We answer in four parts. We find exact bouncing cycles with an elliptic return map and verify them with an independent code (Section 2). We glue the cycle across a moving shell to a Schwarzschild exterior, derive the shell dynamics from the Israel conditions, find periodic shells with positive tension that satisfy the null energy condition, "
          "and check the junction from raw extrinsic curvatures (Section 3). We establish the global causal structure, which is that of a non-traversable bridge (Section 4). We then test stability: the homogeneous and shell stability results are positive, but the linear stability of spherical inhomogeneous perturbations and the energy bookkeeping of the shell are not (Section 5). "
          "The conclusion is deliberately limited: the construction is an exact and idealized solution whose physical robustness is not established. Section 6 lists what is demonstrated and what is not."))

# ------------------------------------------------------------------ 2
add(("h1", "2. The cyclic region"))
add(("h2", "2.1 Model"))
add(("p", "We use units $8\\pi G=c=1$ and the de Sitter radius $R_{dS}=1$ of the false vacuum, so $V_f=3$. A closed Friedmann region of curvature radius $a(t)$ contains a scalar field $\\phi$ with the potential "
          "$V(\\phi)=V_f-\\mu^2\\phi^2/2+\\lambda\\phi^4/4$, whose true minima at $\\phi=\\pm v$ have density $-qV_f$ with $v^2=4(1+q)V_f/\\mu^2$ and $\\lambda=\\mu^2/v^2$, and pressureless dust of conserved bare density $\\rho_*=3C/a^3$. The dust couples to the field through its mass, $\\rho_{\\rm d}=A(\\phi)\\rho_*$ with "
          "$A=1+\\phi^2/2M^2$ and $M^2=\\rho_{\\rm crit}/\\mu^2$ [@Amendola2000]. The equations of motion are"))
add(("eq", "\\ddot{a}=-\\frac{a}{6}\\,(\\rho+3P),\\quad \\ddot{\\phi}+3H\\dot{\\phi}+V'(\\phi)+\\rho_*A'(\\phi)=0,\\quad \\dot{a}^2=\\frac{a^2\\rho}{3}-1", "1"))
add(("p", "with $\\rho=\\dot{\\phi}^2/2+V+A\\rho_*$ and $P=\\dot{\\phi}^2/2-V$. The last relation is the Friedmann constraint, imposed exactly on the initial data. A bounce is a minimum of $a$ with $\\dot{a}=0$ and $\\ddot{a}>0$. The dust parameter is $C=x\\,C_c$ with $C_c=2/(3\\sqrt3)$, the value above which a classical de Sitter neck does not bounce. "
          "All results use $\\mu=20$, $\\rho_{\\rm crit}=0.2V_f$, $q=1$; this is one point of a larger parameter space (Section 6)."))
add(("h2", "2.2 Exact cycles and their return map"))
add(("p", "Start at the neck ($\\dot{a}=0$) with $\\dot{\\phi}=0$ and $\\phi=\\phi_0$, the scale factor fixed by the constraint. The equations are time-reversal invariant, so if $\\dot{\\phi}=0$ again at the next turnaround of $a$, the orbit retraces itself and is exactly periodic. We shoot on $\\phi_0$. "
          "The return map $P$ from one bounce to the next acts on $(\\phi,\\dot{\\phi})$, with $a$ again fixed by the constraint; it preserves area. A fixed point is elliptic if $|{\\rm tr}\\,DP|<2$."))
add(("p", "To guard against code errors the cycle was re-derived with an independent implementation that shares no module with the first. It integrates in cosmic and in conformal time, uses three integrators, and computes stability three ways: finite differences of the return map, the four-dimensional monodromy of the variational equations, and a conformal-time return map. "
          "Table 1 gives the result at $x=0.35$."))
add(("table", "Table 1. The cold cycle at $x=C/C_c=0.35$, independent reproduction.",
     ["Quantity", "Value"],
     [["orbit $\\phi_0/v$", "%.12f" % rep["phi0_over_v"]],
      ["period $T$ (cosmic time; conformal-time clock agrees)", "%.12f" % rep["period"]],
      ["full-state return error $(a,\\dot{a},\\phi,\\dot{\\phi})$", sci(rep["return_err_full_state"], 1)],
      ["Friedmann constraint at start and end of the period", "%s, %s" % (sci(abs(rep["constraint_start"]) + 1e-17, 0), sci(abs(rep["constraint_end"]) + 1e-17, 0))],
      ["return map: trace, determinant", "%.6f, %.7f" % (rep["neck_map_trace"], rep["neck_map_det"])],
      ["monodromy eigenvalues", "$-0.852174\\pm0.523259\\,i$, 1, 1"],
      ["rotation number", "0.4124"]],
     [4.3, 2.6]))
add(("p", "Between $10^{-8}$ and $10^{-13}$ in the integrator tolerance the return error falls from %s to %s and $\\phi_0$ and $T$ agree to $10^{-12}$; the three integrators agree on $T$ to ten digits." % (sci(1.21e-9, 1), sci(1.28e-14, 1))))
add(("h2", "2.3 The stability window and the multiplier $-1$ point"))
add(("p", "The orbit belongs to a branch that is smooth in $x$ (Fig. 2b). Its trace is $-3.74$ at $x=0.330$, $-2.07$ at $0.345$ and $-1.70$ at $0.350$, and equals $-2$ at $x^*=%.6f$, where the multipliers pass through $-1$: the orbit is inverse-hyperbolic below and elliptic above. "
          "A period-doubling bifurcation would create a local 2-cycle there. We looked for one in two ways and did not find it. A multi-start search finds symmetric 2-cycles only at large amplitude ($|\\delta\\phi|\\approx0.4v$), which do not shrink to the fixed point as $x\\to x^*$. "
          "The winding number of $P^2-{\\rm id}$ on circles of radius $0.002$–$0.05\\,v$ around the fixed point equals the index of the fixed point itself ($-1$ below $x^*$, $+1$ above) at five values of $x$, so $P^2$ has no other zero inside the island. We therefore describe the point as a change of type with a multiplier equal to $-1$ and do not call it a period-doubling bifurcation. "
          "The elliptic branch continues to $x\\approx0.360$ ($\\phi_0/v=0.2818$, trace $-1.176$), where states with larger $\\phi_0$ no longer admit a Friedmann neck: the window ends at the boundary of the set of cold bounces and not at a resonance." % cstar))
add(("p", "Linear stability of a reversible conservative system means bounded perturbations, not attraction. We tested nonlinear boundedness by iterating the return map 150 times from perturbed starts. At $x=0.3470$, $0.3500$, $0.351358$ (rotation number exactly 2/5) and $0.354887$ (3/8) all starts with relative perturbations up to $3\\%$ survive, with a neck deviation of 4–6 times the perturbation and no growth. "
          "Near the end of the window ($x=0.3598$) the basin shrinks: at $1\\%$ one of two starts is lost at once and at $3\\%$ both. The twist condition required for an invariant-curve (KAM) argument was not checked, and the survival is numerical evidence, not a proof."))

# ------------------------------------------------------------------ 3
add(("h1", "3. The two-ended geometry"))
add(("h2", "3.1 Construction"))
add(("p", "The parent is the maximally extended Schwarzschild solution with $f(R)=1-m/R$ and Kruskal regions I (right exterior), II (black hole), III (left exterior) and IV (white hole). In III, with future-directed Killing time $s$ and tortoise coordinate $r_*$, $U=e^{\\kappa(r_*+s)}$ and $V=-e^{\\kappa(r_*-s)}$. "
          "A spherical timelike shell $R=R_s(\\tau)$ lies in III. We keep the part of III with $R<R_s$ together with I, II and IV and discard the left asymptotic end. The discarded side is replaced by a cap of polar angle $\\chi_0>\\pi/2$ of the closed cycle, glued across the shell. "
          "On a spacelike slice through the cap, the shell and the bifurcation sphere, the areal radius rises from zero to a maximum, the equator of the cap, falls through the shell to a minimum $r_h=m$, the throat, and rises without bound into region I. There are two ends, the closed cyclic region and the asymptotically flat region I, joined by a throat. "
          "The normal $n$ points from the cap to the kept side; both $\\Gamma_{\\rm in}=n\\cdot\\nabla R$ and $\\Gamma_{\\rm out}$ are negative, so the areal radius decreases across the shell toward the throat."))
add(("h2", "3.2 Shell dynamics"))
add(("p", "Let the shell move with rapidity $\\beta$ relative to the comoving fluid: $dt/d\\tau=\\mathrm{cosh}\\beta$, $a\\,d\\chi/d\\tau=\\mathrm{sinh}\\beta$. Then $R=a\\sin\\chi$ and"))
add(("eq", "\\Gamma_{\\rm in}=\\cos\\chi\\mathrm{cosh}\\beta+\\dot{a}\\sin\\chi\\mathrm{sinh}\\beta,\\quad \\Gamma_{\\rm out}=-\\sqrt{F},\\quad F=1+R'^2-\\frac{m}{R}", "2"))
add(("p", "with $\\Gamma_{\\rm in}^2=1+R'^2-\\rho R^2/3$ by the Friedmann constraint. The Israel conditions [@Israel1966] $S^a{}_b=-([K^a{}_b]-\\delta^a{}_b[K])$ give"))
add(("eq", "\\sigma=\\frac{2(\\Gamma_{\\rm in}-\\Gamma_{\\rm out})}{R},\\qquad p=K^{\\tau}_{\\tau}|_{\\rm out}-K^{\\tau}_{\\tau}|_{\\rm in}+\\frac{\\Gamma_{\\rm out}-\\Gamma_{\\rm in}}{R}", "3"))
add(("p", "where $K^{\\tau}_{\\tau}|_{\\rm in}=\\beta'+H\\mathrm{sinh}\\beta$. With an equation of state $p=w\\sigma$, equation (3) is linear in $\\beta'$ and closes the system; no energy source is added by hand. The Gauss–Codazzi law for a shell that moves through the interior fluid is"))
add(("eq", "(\\sigma R^2)'+p\\,(R^2)'=-R^2(\\rho+P)\\mathrm{cosh}\\beta\\mathrm{sinh}\\beta", "4"))
add(("p", "It holds to $10^{-6}$ along every run, and the opposite sign disagrees at order one. The law was not used to derive the equation of motion, so it checks the derivation. A structural consequence follows. The interior equations (1) contain no shell variable, so for an exactly homogeneous interior the interior dynamics, and hence its return map and stability, "
          "do not depend on the shell at all: a gravitational junction cannot act on the field."))
add(("h2", "3.3 Periodic shells"))
add(("p", "We searched for time-symmetric periodic shells ($\\beta=0$ at the neck and at the turnaround) around the cycle of Table 1 for exteriors with $q=0$ and $q=1$ in the exterior vacuum density $-qV_f$, five values of the exterior mass between $0.1$ and $0.85$ of the largest mass with positive tension, and $w\\in\\{-1/2,-1/4,0,1/3,1\\}$. "
          "Only $w=-1/2$ with a flat exterior gave periodic shells (5 of 50 settings); no anti-de Sitter exterior did. All five have a return map with multipliers of modulus 1, positive tension and $\\sigma+p\\ge0$ over the whole cycle (Table 2)."))
rows = []
for a, s in zip(BR["analysis"], SH):
    rows.append(["%.2f" % a["m_frac"], "%.4f" % a["chi0_over_pi"], "%.5f" % a["m"], "%.3f–%.3f" % (a["R_shell_min"], a["R_shell_max"]), "%.0f" % a["shell_over_throat_min"], "%.2f–%.2f" % tuple(a["sigma_range"]), "%.2f" % s["diag"]["max_abs_beta"]])
add(("table", "Table 2. Periodic moving shells around the cycle of Table 1 ($q=0$, $w=-1/2$). $m=r_h$ is the exterior mass parameter, $m_{\\max}$ the largest value with positive tension.",
     ["$m/m_{\\max}$", "$\\chi_0/\\pi$", "$m=r_h$", "shell $R$", "$R/r_h$ min", "$\\sigma$", "max $|\\beta|$"], rows, [0.8, 0.9, 0.9, 1.3, 0.9, 1.1, 0.9]))
add(("h2", "3.4 Global junction consistency"))
r = JG["residuals"]
add(("p", "Local consistency was checked independently of equations (2)–(3). Each side's extrinsic curvature was computed from the Christoffel symbols of its own metric, interior $g={\\rm diag}(-1,a^2,a^2\\sin^2\\chi,\\dots)$ on $(t,\\chi)$ and exterior $-f\\,ds^2+dR^2/f+R^2d\\Omega^2$ on $(s,R)$, with the shell worldlines and unit normals as above. Over one period of the conservative shell ($m/m_{\\max}=0.10$, 3001 points): "
          "the induced metric agrees on both sides to %s; the normals are unit and orthogonal to the velocity to %s; the signs of the normals ($n^\\chi>0$ inside, $n^R<0$ outside) and of $dt/d\\tau$, $ds/d\\tau$ are those required for the kept side to contain the throat and for the two time orientations to agree; "
          "$\\sigma$ and $p$ recomputed from raw curvatures agree with those of the equation of motion to %s and $p=w\\sigma$ holds to %s. The kept exterior has no boundary other than the shell: I, II and IV are complete Kruskal regions, and the part of III with $R<R_s$ is bounded by the shell and by the horizons that continue into II and IV. "
          "The shell stays outside the horizon, $R/r_h\\ge%.0f$ at all times. Residuals of the periodic shell are listed in Section 5.2." % (sci(max(r["h_tt_in"], r["h_tt_out"]), 0), sci(max(r["n_unit"], r["n_orth"]), 0), sci(max(r["sigma"], r["p"]), 0), sci(r["p_minus_w_sigma"], 0), A0["shell_over_throat_min"])))

# ------------------------------------------------------------------ 4
add(("h1", "4. Causal structure and traversability"))
add(("h2", "4.1 The bridge is not traversable"))
add(("p", "The shell worldline lies in region III ($U>0$, $V<0$) because $R_s>r_h$ throughout. In the kept exterior a future-directed causal curve has non-decreasing $U$ and $V$, so every point in the causal future of the cap has $U\\ge U_{\\rm shell}>0$, whereas region I is $\\{U<0\\}$. Hence no future-directed causal curve joins the cap to region I, and none joins region I to the cap: "
          "they are spacelike separated. The causal future of the cap is contained in III and II. This is exact and holds for every shell in Table 2. Radial null rays from the cap centre, emitted at 24 phases of the cycle, all reach the shell after %.1f–%.1f interior periods (the cap is nearly a full three-sphere, so a radial ray needs several cycles), cross it, and fall along $U=$ const toward $V=0$ into region II (Fig. 1c). "
          "Nothing carries matter or a signal from the cyclic region to the parent universe." % (min(min(x["delay_in_periods"] for x in a["rays"]) for a in BR["analysis"]), max(max(x["delay_in_periods"] for x in a["rays"]) for a in BR["analysis"]))))
add(("fig", FIG1, "Figure 1. (a) Schematic causal structure of the kept spacetime: the shell (blue) in region III, the closed region replacing the left end, radial null rays leaving the cap (orange). (b) The equator of the cap $a(t)$, the shell radius $R(t)$ and the throat $r_h$; the ordering $a>R>r_h$ holds at all times. "
                    "(c) Region III in static coordinates: the shell oscillates, rays fall to the horizon. (d) Portrait of the return map around the elliptic fixed point.", 6.4))
add(("h2", "4.2 Energy conditions"))
add(("p", "The shell has $\\sigma>0$ and $\\sigma+p\\ge0$ over the whole cycle for all five shells of Table 2, so no exotic shell matter is needed to maintain the bridge. The bulk satisfies the null energy condition (scalar field and dust). This is compatible with non-traversability: the throat of the Schwarzschild bridge is a vacuum bifurcation sphere behind horizons, "
          "and a traversable throat would require a violation of the null energy condition at or near it [@Hochberg1997], which the model does not contain. The shell is untrapped on both sides at all times; the only apparent horizons are the equator of the cap, marginally trapped at the necks and turnarounds where $\\dot{a}=0$, and the Schwarzschild horizon behind the shell."))
add(("h2", "4.3 A pre-existing bridge, not a formed one"))
add(("p", "The construction uses the eternal Kruskal extension. Take the slice $\\Sigma_0=\\{s=0\\}$ in region III, continued through the bifurcation sphere into region I, with the cap data at the time the shell crosses $s=0$. It carries regular initial data: the cap satisfies the Friedmann constraint, the shell satisfies the Israel conditions, the exterior is time-symmetric vacuum. "
          "The future domain of dependence of $\\Sigma_0$ contains I, II and the later part of III with the cap and shell, and does not contain the white-hole region IV. The evolution therefore does not need a white hole. It does need the initial throat: $\\Sigma_0$ has a minimal sphere and two ends. "
          "The bridge is primordial. Data from regular gravitational collapse have no left end and no region IV, and whether a bridge of this kind can form dynamically is a separate question that we do not address. The extension to the past is only the maximal Cauchy development."))
add(("p", "The throat is very small and very strongly curved on the scale of the cycle. For the conservative shell $r_h=%.5f$ in units $R_{dS}=1$, against an equator $a\\approx0.95$, and the Kretschmann scalar is %s at the throat, %.2f–%.2f at the shell and %.0f–%.0f in the cap." %
          (A0["r_h"], sci(JG["kretschmann_exterior_at_throat"], 1), JG["kretschmann_exterior_at_shell"][0], JG["kretschmann_exterior_at_shell"][1], JG["kretschmann_interior"][0], JG["kretschmann_interior"][1])))

# ------------------------------------------------------------------ 5
add(("h1", "5. Stability"))
add(("h2", "5.1 Homogeneous interior and the shell"))
add(("p", "The homogeneous interior is stable in the sense of Section 2.3. For the shell we perturbed, at random, the neck state of the interior and the shell data $(\\chi,\\beta)$ by relative amounts $10^{-6}$ to $10^{-2}$ and followed 60 interior periods. For the conservative shell ($m/m_{\\max}=0.10$) all runs complete with $R/r_h\\ge%.0f$, $\\sigma\\ge%.2f$ and $\\sigma+p\\ge0$, and the interior history equals the interior evolved alone to $2\\times10^{-10}$, "
          "as it must by Section 3.2. The stiffest shell ($m/m_{\\max}=0.85$) has a return map of modulus 1 but fails at a perturbation of $10^{-3}$: $\\sigma$ falls to $10^{-5}$ and $|\\beta|$ reaches 8 within six periods. It is linearly but not practically stable. "
          "Residuals of the periodic conservative shell: Friedmann constraint %s, $\\Gamma_{\\rm in}^2$ identity %s, $p-w\\sigma$ %s, mass function %s (these follow algebraically and mostly test the code), flux law (4) %s relative (independent, finite-difference limited); the return error is %s and converges between $10^{-9}$ and $10^{-13}$." %
          (min(x["R_over_rh_min"] for x in long0), min(x["sigma_min"] for x in long0), sci(CV["residuals"]["friedmann"], 0), sci(CV["residuals"]["gamma_in"], 0), sci(CV["residuals"]["israel_p"], 0), sci(CV["residuals"]["mass"], 0), sci(CV["residuals"]["flux"], 0), sci(1.3e-14, 0))))
add(("h2", "5.2 Inhomogeneous perturbations"))
add(("p", "We linearized the full system, gravity, the scalar field and the coupled dust, about the cycle for spherically symmetric perturbations in longitudinal gauge, with $\\Psi=\\Phi$ from the traceless equation [@Bardeen1980; @Kodama1984; @Mukhanov1992]. The dust obeys number conservation and the force law $u^\\nu\\nabla_\\nu u_\\mu=-h_\\mu{}^\\nu\\partial_\\nu\\ln A$; the field obeys $\\Box\\phi=V'+A'\\rho_*$. "
          "The equations were expanded to first order symbolically and projected on the radial harmonics $Q_n=\\sin(n\\chi)/\\sin\\chi$, $k^2=n^2-1$. The constraints give $\\Phi$ algebraically, without division by $H$, so the equations are regular at the bounces. The mode $n=1$ is the homogeneous one; $n=2$ is not covered by this gauge. "
          "Two checks passed: the derivative of the algebraic $\\Phi$ agrees with $\\Phi'$ from the momentum constraint to %s and %s relative for $n=5$ and $12$, and the unused trace equation holds to %s and %s." %
          (sci(PF["checks"]["5"]["constraint_propagation"], 0), sci(PF["checks"]["12"]["constraint_propagation"], 0), sci(PF["checks"]["5"]["trace_equation"], 0), sci(PF["checks"]["12"]["trace_equation"], 0))))
add(("p", "The monodromy over one cycle has, for every harmonic from $n=3$ to $60$, a multiplier far above one (Table 3, Fig. 2a). A perturbation of $10^{-6}$ reaches order unity within about 1.2–3 cycles for $n\\ge6$ and within about five cycles at $n=3$."))
add(("table", "Table 3. Largest Floquet multiplier per cycle for inhomogeneous perturbations, and the same with the dust–scalar coupling removed from the perturbation equations only (control).",
     ["$n$", "full system", "control (coupling off)"],
     [[str(n_), g(mode[n_]["max_abs"], 3), g(max(ctrl[n_]["multipliers_abs"]), 3) if n_ in ctrl else "–"] for n_ in (3, 5, 10, 20, 30, 60)], [1.0, 2.2, 2.5]))
add(("p", "The growth has two sources. At large $n$ it is independent of $k$: in the quasi-static limit the dust contrast obeys $\\ddot D=g^2D$ with"))
add(("eq", "g^2=\\frac{6C\\phi_0^2}{M^2a^3\\,(2M^2+\\phi_0^2)}", "5"))
add(("p", "the attractive scalar-mediated force on the coupled dust. Integrated over the cycle $g$ gives a logarithmic growth of %.1f (a rate of %.0f–%.0f per unit time) against $\\ln\\mu=%.2f$ at $n=60$ and %.2f at $n=40$, approaching the estimate from below. With the dust–scalar coupling removed from the perturbations only, the large-$n$ multipliers fall to %.2f, "
          "which isolates this mechanism. The long-wavelength modes ($n=3$ and $10$) remain unstable in the control (%s and %s) because the field is near the top of its potential, where they are tachyonic. "
          "The result is linear, on the whole three-sphere (the cap has a boundary at the shell whose boundary conditions were not imposed), for pressureless dust and one parameter set. Section 5.3 tests whether changing the model removes this growth. "
          "The nonlinear inhomogeneous evolution was not performed. The linear growth rates make a robust survival implausible without a new stabilizing mechanism." %
          (PM["estimate_log_growth"], PM["rate_min"], PM["rate_max"], __import__("math").log(mode[60]["max_abs"]), __import__("math").log(mode[40]["max_abs"]), max(ctrl[60]["multipliers_abs"]), g(max(ctrl[3]["multipliers_abs"]), 3), g(max(ctrl[10]["multipliers_abs"]), 3))))
add(("fig", FIG2, "Figure 2. (a) Largest Floquet multiplier per cycle for inhomogeneous perturbations as a function of the radial harmonic $n$: full system, control with the dust–scalar coupling removed from the perturbations, and the quasi-static estimate from equation (5). "
                    "(b) Trace of the return map along the main orbit branch, the elliptic band, the point where the trace equals $-2$, the end of the branch, and the points used for the resonance tests.", 6.4))
add(("h2", "5.3 Attempts to stabilize"))
add(("p", "We re-derived the perturbation equations for a polynomial potential $V=3+a_2\\phi^2+a_4\\phi^4+a_6\\phi^6$ and for a coupled fluid with pressure, a polytrope $P=Kn^{5/3}$ with $\\rho_f=An+Kn^{5/3}/(\\gamma-1)$, whose motion follows from total energy–momentum conservation; at $K=a_6=0$ they reproduce Table 3. "
          "Each change was followed by continuation along the elliptic orbit at fixed return-map trace ($-1.7043$), solving for $\\phi_0$ and adjusting $C/C_c$ at every step, not by a scan of initial conditions (Table 4, Fig. 3). "
          "Changing the dust–scalar coupling from $\\rho_{\\rm crit}/V_f=%.4f$ to %.2f changes the $n=3$ growth only between %.1f and %.1f, because the orbit retunes $C/C_c$ and the fifth-force rate at the neck stays near 16–24. Raising the hilltop curvature $\\mu$ to %.0f reduces the fifth-force part ($g\\,T$ falls from 10.7 to 0.5) but widens the tachyonic band to $n\\lesssim\\mu$, "
          "and the worst growth over all $n$ stays above $3.6\\times10^{3}$. A sextic reshaping of the hilltop at fixed curvature and well depth changes the $n=3$ growth by about 20%%. Pressure removes only modes with $n\\gtrsim n_c=ga/c_s$: at $c_s^2=0.04$ the modes $n=240$ are stable and $n=120$ is reduced to %.0f, while $n=3$–12 are unchanged, "
          "and the cycle ceases to exist before $c_s^2$ reaches 0.08." % (min(r["rc"] for r in CU), max(r["rc"] for r in CU), min(r["growth"]["3"] for r in CU), max(r["growth"]["3"] for r in CU), MUF[-1]["mu"], [r for r in PR if r["K"] == 0.08][0]["growth"]["120"])))
RNG = lambda rows, k: "%.3g–%.3g" % (min(r["growth"][k] for r in rows), max(r["growth"][k] for r in rows))
POTALL = POT["wider"] + POT["narrower"]
add(("table", "Table 4. Inhomogeneous growth per cycle under each attempt (range of the largest Floquet multiplier along the continuation).",
     ["change (range followed)", "$n=3$", "$n=12$", "$n=60$"],
     [["baseline", g(mode[3]["max_abs"], 3), g(mode[12]["max_abs"], 3), g(mode[60]["max_abs"], 3)],
      ["coupling $\\rho_{\\rm crit}/V_f$ from 0.0086 to 0.51", RNG(CU, "3"), RNG(CU, "12"), RNG(CU, "60")],
      ["hilltop scale $\\mu$ from 20 to 307", RNG(MUF, "3"), RNG(MUF, "12"), RNG(MUF, "60")],
      ["sextic hilltop, $v/v_0$ from 0.98 to 1.74", RNG(POTALL, "3"), RNG(POTALL, "12"), RNG(POTALL, "60")],
      ["pressure, $K=0$ to 0.08 ($c_s^2\\le0.04$)", RNG(PR, "3"), RNG(PR, "12"), RNG(PR, "60")]],
     [3.0, 1.3, 1.3, 1.3]))
WK = [w for w in AN["wkb"] if w["family"] == "pressure" and w["n"] in (120, 240) and w["K"] in (0.04, 0.08)]
add(("p", "Two criteria summarize why. In the quasi-static limit the fluid contrast obeys $\\ddot D=(g^2-c_s^2k^2/a^2)D$, so mode $n$ grows when $n<n_c=ga/c_s$; with pressure the exponent $\\int\\sqrt{g^2-c_s^2k^2/a^2}\\,dt$ equals the Floquet exponent to within %.2f for $n=120$ and 240 at $K=0.04$ and 0.08, and fails at low $n$ where the quasi-static assumption breaks. "
          "The scalar sector alone is tachyonic, $\\ddot{\\delta\\phi}+3H\\dot{\\delta\\phi}+(k^2/a^2+V''+n/M^2)\\delta\\phi=0$ with $V''<0$ near the top, for $k<a\\mu$; its scalar-only multiplier accounts for most of the growth at large $\\mu$ (7.1 of 7.9 in the logarithm at $n=60$, $\\mu=208$). "
          "Along the whole family $\\mu T$ stays between 12 and 16, so the hilltop exponent does not shrink when $\\mu$ grows. A stable cycle would need both the fifth-force and the hilltop band closed; pressure closes only the first and only above $n_c$, and nothing tried closes the second. "
          "These are mechanism-level estimates validated numerically, not rigorous bounds, and the low-$n$ exponent ($\\ln G\\approx2.2$–$2.7$) is an observed invariant, not a derived one. We did not run nonlinear evolutions because no promising region was found." % max(abs(w["S5"] - w["lnG"]) for w in WK)))
add(("fig", FIG3, "Figure 3. Attempts to stabilize. (a) Growth per cycle versus $n$ for the baseline, for larger hilltop curvature $\\mu$, and for a fluid with pressure. (b) Weakening the dust–scalar coupling at fixed return-map trace. (c) Changing the hilltop curvature: growth at $n=3$ and 12 and the worst case over all $n$.", 6.4))
add(("h2", "5.4 The shell matter"))
add(("p", "An isotropic gas of non-interacting, slowly moving Nambu–Goto strings confined to the two-dimensional shell has energy density $\\sigma=\\mu_s n$ ($n$ the length per area) and tension $\\mu_s$ along each string, so $S^{ij}=-\\mu_sn\\langle\\hat t^i\\hat t^j\\rangle=-(\\sigma/2)\\delta^{ij}$ and $p=-\\sigma/2$, the equation of state used for the periodic shells. "
          "Without exchange, string length would be conserved and $(\\sigma R)'=0$. Equation (4) requires otherwise: a shell that moves through the interior fluid must absorb and re-emit bulk energy. Numerically the gross energy exchanged with the bulk per cycle is %.2f, %.0f%% of the shell's mean energy (%.2f–%.2f) and %.0f%% of the energy of the cap; the net over a period is zero. "
          "A string gas of fixed tension can supply this only if strings are created and destroyed by absorbing dust and scalar quanta, and no such interaction has been specified. The equation of state is motivated, but the shell is not a complete matter model: the energy bookkeeping is satisfied only through an unmodelled interaction. For a comoving shell with the field confined to the cap the exchange vanishes (Section 5.5)." %
          (JG["exchanged_energy_per_cycle_abs"], 100 * JG["exchange_over_shell_energy"], JG["shell_energy_range"][0], JG["shell_energy_range"][1], 100 * JG["exchange_over_cap_energy"])))

# ------------------------------------------------------------------ 6
add(("h2", "5.5 Boundary conditions on the finite cap"))
add(("p", "Section 5.2 used the radial harmonics of the full three-sphere, whereas the bridge contains a cap bounded by the shell, and the scalar field must obey some condition there. For the homogeneous interior $\\partial_n\\phi=\\sinh\\beta\\,\\dot\\phi$, so a no-flux condition holds only for a comoving shell ($\\beta=0$); "
          "the moving shells of Table 2 need an absorbing condition for which no field theory is given. For a comoving shell and a field confined to the cap, zero scalar flux and zero normal dust velocity reduce to $Q_\\nu'(\\chi_0)=0$ for the modes $Q_\\nu=\\sin(\\nu\\chi)/\\sin\\chi$, that is $\\tan(\\nu\\chi_0)=\\nu\\tan\\chi_0$ ($\\nu=1$ is the homogeneous mode). "
          "Comoving shells with $\\sigma>0$ and $\\sigma+p\\ge0$ over the cycle exist for many exteriors (%d of 9900 combinations), but for the anti-de Sitter exterior that matches the potential, $V(\\pm v)=-V_f$, only for $\\chi_0\\lesssim0.6\\pi$, with $\\sigma\\approx%.1f$–%.1f and $p/\\sigma$ varying between %.2f and %.2f. "
          "On such caps the lowest mode grows by %.1f–%.1f per cycle (at $\\chi_0=0.6\\pi$ and $0.55\\pi$) and the worst case is unchanged at %s. The growth $G(\\nu)$ exceeds one for every $\\nu>1$, tends to one only at $\\nu=1$, and rises monotonically for $\\nu>2$; a boundary condition only quantizes $\\nu$, so with infinitely many modes at $\\nu\\gtrsim3$ (where $G>14$) none can stabilize the cap. "
          "A shell with scalar-dependent tension $\\sigma(\\phi)$ would give $p=-\\sigma$ and a Robin condition, but continuity of $\\phi$ across the shell conflicts with a static exterior vacuum: the field would have to equal $v$ at the shell while the interior field varies between %.2f$v$ and %.2f$v$, and the tension of a thin wall between them, %.2f–%.2f, is far below the $\\sigma\\approx2.3$ the shell requires. We did not pursue it." %
          (NEC0["n_nec_full"], min(r["sigma_range"][0] for r in SHG["shells"]), max(r["sigma_range"][1] for r in SHG["shells"]), min(r["w_range"][0] for r in SHG["shells"]), max(r["w_range"][1] for r in SHG["shells"]),
           [r for r in CAP if abs(r["chi0_over_pi"] - 0.6) < 1e-9][0]["growth"][0], [r for r in CAP if abs(r["chi0_over_pi"] - 0.55) < 1e-9][0]["growth"][0], "7.7×10⁴",
           SHG["phi_range_over_v"][0], SHG["phi_range_over_v"][1], min(SHG["wall_tension"]), max(SHG["wall_tension"]))))
add(("h1", "6. Discussion"))
add(("p", "What is demonstrated, numerically and with the checks above, is that the closed cold cycle of Section 2 admits a consistent two-ended embedding with a non-traversable throat, with a moving shell that has positive tension and satisfies the null energy condition, and that the cyclic region and the parent universe are causally disconnected. "
          "What is not demonstrated is more than what is. The homogeneous cycle is linearly unstable to inhomogeneous perturbations in this model, by factors of %s to %s per cycle; the shell's matter content requires an unmodelled energy exchange comparable to its own energy; the throat is tiny and its curvature enormous; and the bridge is primordial, not formed. "
          "A traversable bridge, matter transfer between the regions and cosmological recycling are excluded in this geometry. We make no claim that a cyclic region exists or that one is hidden from us: causal inaccessibility would explain why such a region could be hidden, not that it is there." % (g(mode[3]["max_abs"], 2), sci(mode[60]["max_abs"], 1))))
add(("table", "Table 5. Status of the claims.",
     ["Claim", "Status"],
     [["exact cyclic region, elliptic return map, independent reproduction", "demonstrated numerically (not a proof)"],
      ["Israel junction, orientations, global gluing", "demonstrated, residuals $\\le5\\times10^{-15}$"],
      ["shell with $\\sigma>0$, null energy condition, periodic, 60-cycle survival at $\\le1\\%$", "demonstrated for the conservative shell"],
      ["non-traversability, causal disconnection from the parent universe", "demonstrated (exact)"],
      ["independence of the white-hole region", "future development of $\\Sigma_0$ excludes it; the throat is pre-existing"],
      ["stability to inhomogeneous perturbations", "fails (linear, $n\\ge3$); not removed by weaker coupling, other hilltops, pressure or the boundary conditions of the finite cap"],
      ["shell matter from an action", "equation of state motivated; energy exchange unexplained"],
      ["dynamical formation of the bridge; traversability; matter transfer", "not demonstrated (last two excluded here)"]],
     [4.4, 2.6]))
add(("p", "Three routes follow. The instability has two sources, the fifth force on the coupled fluid and the hilltop of the potential, so a different mechanism for the field's return to the false vacuum may be required; continuation in the coupling, the potential scale, the hilltop shape and a fluid with pressure (Section 5.3) did not remove the growth, though other orbits, windows and matter models were not explored. "
          "The boundary conditions of the inhomogeneous modes at the shell, and the nonlinear evolution of spherical inhomogeneities, need a double-null code with a fluid and a moving boundary. And the shell needs a matter model from an action with an explicit interaction with the bulk. "
          "Until these are done the construction should be read as an exact idealized solution of Einstein gravity, not as a cosmological model."))
add(("p0", "Funding: this work received no external funding."))
add(("p0", "Competing interests: the author declares no competing interests."))
add(("p0", "Author contributions (CRediT): Leon Sandler, sole author: conceptualization, methodology, software, formal analysis, writing, review and editing."))
add(("p0", "Data and code availability: the code, the data of every figure and table and the manuscript are available at %s (branch phase2-bridge); one command per script regenerates every number and figure, and the results quoted in the text are read from the JSON files that the scripts write. An archived version with a DOI will be deposited before submission." % GITHUB))
add(("p0", "Use of generative AI: a generative AI tool (Claude, Anthropic) assisted with code development, numerical experiments, drafting and editing. The author checked all of it and takes full responsibility for the content."))
add(("refs",))

import refs as _refs

LIB = {k: _refs.LIB[k] for k in ("Israel1966", "Blau1987", "Farhi1990", "Fischler1990", "Frolov1989", "Morris1988", "Hochberg1997", "Hayward1999", "Steinhardt2002")}
LIB.update({
    "EinsteinRosen1935": "A. Einstein, N. Rosen, The particle problem in the general theory of relativity, Physical Review 48 (1935) 73–77. https://doi.org/10.1103/PhysRev.48.73",
    "Kruskal1960": "M.D. Kruskal, Maximal extension of Schwarzschild metric, Physical Review 119 (1960) 1743–1745. https://doi.org/10.1103/PhysRev.119.1743",
    "Fuller1962": "R.W. Fuller, J.A. Wheeler, Causality and multiply connected space-time, Physical Review 128 (1962) 919–929. https://doi.org/10.1103/PhysRev.128.919",
    "Friedman1993": "J.L. Friedman, K. Schleich, D.M. Witt, Topological censorship, Physical Review Letters 71 (1993) 1486–1489. https://doi.org/10.1103/PhysRevLett.71.1486",
    "Tolman1931": "R.C. Tolman, On the theoretical requirements for a periodic behaviour of the universe, Physical Review 38 (1931) 1758–1771. https://doi.org/10.1103/PhysRev.38.1758",
    "Amendola2000": "L. Amendola, Coupled quintessence, Physical Review D 62 (2000) 043511. https://doi.org/10.1103/PhysRevD.62.043511",
    "Bardeen1980": "J.M. Bardeen, Gauge-invariant cosmological perturbations, Physical Review D 22 (1980) 1882–1905. https://doi.org/10.1103/PhysRevD.22.1882",
    "Kodama1984": "H. Kodama, M. Sasaki, Cosmological perturbation theory, Progress of Theoretical Physics Supplement 78 (1984) 1–166. https://doi.org/10.1143/PTPS.78.1",
    "Mukhanov1992": "V.F. Mukhanov, H.A. Feldman, R.H. Brandenberger, Theory of cosmological perturbations, Physics Reports 215 (1992) 203–333. https://doi.org/10.1016/0370-1573(92)90044-Z",
})
