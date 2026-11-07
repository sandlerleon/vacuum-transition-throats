# -*- coding: utf-8 -*-
"""Front matter, introduction, model, null-energy argument and the thin-wall analysis (revision 1)."""
from pcommon import *

R_dS = TW["R_dS_Gpc"]
M1 = [r for r in TW["mass_table"] if r["q"] == 1.0][0]
SIGMA = TW["sigma_model_hat"]
GRAVP = TW["gravity_parameter_model"]
THICK = TW["thickness_over_RdS"]
PER = VAL["per_A"]
TST = VAL["vacuum_tests"]
SC = VAL["scan"]
FT = TW["finite_tension"]
CRE = VE["creation"]
resolved = sum(1 for r in VE["inconclusive"] if any(v["first_trapped"] is not None for k, v in r.items() if k.startswith("refined")))
N_TRAP = SC["n_trapped"] + resolved
N_UNRES = SC["n_runs"] - N_TRAP

add(("title", TITLE))
add(("author", "Leon Sandler"))
add(("small", "Independent Researcher, Northbrook, Illinois, USA"))
add(("small", "ORCID: https://orcid.org/0009-0007-4584-808X"))
add(("small", "Corresponding author: sandler.leon@gmail.com"))
add(("blank",))
add(("h1", "Abstract"))
ABSTRACT = (
    "Rebirth scenarios invoke a wormhole carrying matter into a new expanding region. We ask whether a density-triggered vacuum transition can form one in a scalar-only gravitational benchmark: a canonical scalar whose false vacuum has the observed "
    "dark-energy density, motivated by a dark-matter coupling that the numerical experiments do not evolve. In this classical model the null energy condition holds, so no ordinary traversable throat can form; nontraversable bridges and quantum or baby-universe "
    "continuations are not excluded. For a thin spherical wall and regular initial data the Israel junction conditions give a finite-tension critical exterior mass parameter, equal to %s up to corrections of order $\hat\sigma^2$, "
    "below which the configurations on the branch analyzed collapse; an independent derivation confirms it. A spherical Einstein–scalar code checked against exact solutions found trapped surfaces in %d of %d runs and in %d more after refinement. "
    "An expanding (anti-trapped) false-vacuum core appears only when the initial data already contain a core larger than its own de Sitter radius; none develops in the %d runs that start smaller: persistence, not creation. "
    "At the observed dark-energy scale the core radius is %s Gpc and the exterior mass parameter %s solar masses (for $q=1$ half of it is quasilocal at the wall), dominated by the vacuum-energy difference. In the computed double-null domain no causal curve carries exterior matter into the core; "
    "the global structure is not computed. The evidence is complete only for anti-de Sitter exteriors and walls much thicker than the physical one."
    % ("$(1+q)R_{dS}$", SC["n_trapped"], SC["n_runs"], resolved, CRE["n_initially_normal_core"], n(R_dS, 1), "%s × 10$^{%d}$" % (n(M1["M_c_Msun"] / 10 ** 23, 1), 23)))
add(("p0", ABSTRACT))
add(("h1", "Keywords"))
add(("p0", "Wormholes; Vacuum decay; Gravitational collapse; Black hole physics; Cosmology: theory; Methods: numerical"))

# ------------------------------------------------------------------------------------------------ 1
add(("h1", "1. Introduction"))
add(("p", "In a universe that dilutes its matter and approaches a cold, dark-energy-dominated state, any mechanism that returns the cosmos to a hot, dense phase has to do two things: change the vacuum, and supply a geometry "
          "through which a new expanding region can be reached. Cyclic models [@Steinhardt2002] and cosmological natural selection [@Smolin1992] place that geometry in a bounce or in a black-hole interior. A more literal "
          "proposal is a wormhole: a spacetime throat that forms near the end of the universe and carries matter back to a new beginning. General relativity does not establish that a wormhole behaves as a gravitational "
          "vacuum cleaner, and ordinary traversable wormholes need stress-energy that violates the null energy condition [@Morris1988; @Hochberg1997; @Hayward1999]. A defensible version of the idea must therefore be "
          "stated as a set of mathematical questions that can fail."))
add(("p", "This paper states and answers such questions for one model. The model is a scalar field $\\phi$ with a double-well potential whose false vacuum is held in place by dark matter through a symmetron-type "
          "coupling [@Hinterbichler2010], so that the observed dark-energy density is carried by the false vacuum while the matter density is above a critical value, and released once the expansion has diluted the matter "
          "(Section 2). The release can start locally, by the nucleation of a bubble of the lower vacuum, and the reverse process, the return of a region to the false vacuum when its matter density rises above the "
          "critical value, can also be local. The question is what the gravitational collapse of such a region can produce. The hypothesis, stated so that it can fail, is that a density-triggered vacuum transition may "
          "create a localized gravitational configuration that, under explicit conditions, develops a spacetime throat connecting a contracting region to a new expanding daughter region."))
add(("p", "Three outcomes are distinguished throughout and are not allowed to blur. In outcome A the collapse is ordinary: a trapped region forms and no regular expanding region does. In outcome B a bridge-like "
          "geometry exists but nothing can cross it. In outcome C an expanding region exists behind a horizon. A collapsing region can produce a trapped surface without producing either a traversable wormhole or a new "
          "universe, and the distinction is the content of the paper. A second distinction is made as sharply: an expanding region that is present in the initial data and persists beside a collapse is not a daughter "
          "universe that the collapse has created. Even a created daughter would establish cosmological reproduction, not a cyclic return of the same universe, which additionally needs a daughter that reconnects "
          "with or replaces its parent; no part of this work provides that. Section 9.1 discusses the alternative that the bridge is primordial rather than formed by the collapse."))
add(("p", "Four terms are used with fixed meanings. An anti-trapped (expanding) region is a region in which both null expansions are positive, $r_u>0$ and $r_v>0$: a local property that the numerical solutions can test. "
          "The expanding core is the false-vacuum interior of the wall when it contains such a region. A bridge is the nontraversable two-ended geometry of the thin-wall solution (Section 4); no bridge is constructed numerically. "
          "A daughter universe is reserved for a global expanding region that the collapse has created and that has a causal continuation; no run produces one. Persistence is the survival of an anti-trapped region that is present in the initial data, "
          "and creation is its appearance from initial data that have none."))
add(("p", "The mechanism is not new in its geometry. Thin false-vacuum bubbles in general relativity were analyzed in detail in the 1980s [@Blau1987; @Berezin1987], the possibility of a baby universe attached to a "
          "black hole by a throat was explored in classical and quantum form [@Frolov1989; @Farhi1990; @Fischler1990], and the obstacle that regular initial data do not lead to one was identified at the time "
          "[@FarhiGuth1987]. Known from that literature are the Israel junction analysis of a false-vacuum bubble, the existence of a critical mass below which regular data collapse, and the obstruction to reaching the baby-universe branch from regular initial data. What this paper adds is quantitative and specific to the model: a null-energy proposition for the scalar-plus-dust model; a finite-tension critical mass with its small-tension limit and a coincidence of horizon areas there, rederived independently (Supplement S8); a regular-center Einstein–scalar experiment that separates persistence from creation, with its dependence on wall thickness quantified (Section 5.6); and dimensioned numbers at the observed dark-energy scale, with the exterior mass parameter distinguished from the quasilocal mass (Section 4). The vacuum energy is the observed dark-energy density, so the critical sizes and masses "
          "are numbers rather than free parameters, and they can be compared with the matter that the universe actually contains. The paper does not add a new bubble solution, a traversable throat or a global spacetime for the numerical solutions. It makes five contributions."))
add(("bullet", "A null-energy argument for the whole model (Section 3): no traversable throat can form, and no null ray has a minimum of the areal radius."))
add(("bullet", "A finite-tension thin-wall criterion (Section 4): a critical exterior mass $m_c(\\hat\\sigma)$ below which a false-vacuum region is confined and collapses, with its small-tension limit $(1+q)R_{dS}$ and a coincidence of "
               "horizon areas at that limit (Section 6)."))
add(("bullet", "A numerical Einstein–scalar experiment with a regular center (Section 5), validated against exact solutions and by the convergence of its constraints, that tests the criterion with thick walls, separates "
               "expanding regions that persist from those that could be created, and compares exteriors of different vacuum depths."))
add(("bullet", "An accounting of the size, mass and energy budget that the criterion requires at the observed dark-energy scale, showing that the critical mass is dominated by the vacuum-energy difference and not by matter (Section 7), and falsification criteria (Section 8)."))
add(("bullet", "Audits of the evidence: an independent derivation of the junction branches, horizons and masses (Supplement S8), a quantitative wall-thickness extrapolation (Section 5.6) and a causal-classification diagnostic that separates the evolved domain from conjecture (Section 5.4)."))

# ------------------------------------------------------------------------------------------------ 2
add(("h1", "2. Model, inputs and scales"))
add(("p", "Units are $8\\pi G = 1$ ($M_{\\rm Pl}$ is the reduced Planck mass) and $c = \\hbar = 1$ unless stated. The action is"))
add(("eq", "S=\\int d^4x\\sqrt{-g}\\left[\\frac{M_{\\rm Pl}^{2}}{2}R-\\frac12(\\nabla\\phi)^2-V_0(\\phi)\\right]+S_{\\rm dm}\\left[A^2(\\phi)g_{\\mu\\nu},\\psi_{\\rm dm}\\right]", "1"))
add(("p", "with $A(\\phi)=1+\\phi^2/2M^2$ the coupling to dark matter alone. For pressureless dark matter of conserved density $\\rho_*\\propto a^{-3}$ the scalar moves in $V_{\\rm eff}=V_0(\\phi)+\\rho_*A(\\phi)$. "
          "While $\\rho_*$ exceeds a critical value $\\rho_{\\rm crit}=\\mu^2M^2$ the symmetric point of $V_0=V_c-\\mu^2\\phi^2/2+\\lambda\\phi^4/4$ is a minimum and holds the field in a false vacuum of density "
          "$V_f$. When the expansion has diluted the matter below $\\rho_{\\rm crit}$ the point becomes a maximum and the field can roll to a lower vacuum. Two numbers fix the gravitational problem, and nothing else of the "
          "particle physics enters it: the false-vacuum density $V_f$, which is set equal to the observed dark-energy density, and the depth of the true vacuum, $V_t=-qV_f$. The baseline is $q=1$; $q>0$ is an anti-de Sitter "
          "true vacuum, $q=0$ a Minkowski one and $-1<q<0$ a lower de Sitter vacuum. The relevant length is the de Sitter radius of the false vacuum,"))
add(("eq", "R_{dS}=\\sqrt{3/V_f}=H_{\\Lambda}^{-1}", "2"))
add(("p", "which for the Planck 2018 values [@Planck2020] $H_0=67.66$ km s$^{-1}$ Mpc$^{-1}$ and $\\Omega_\\Lambda=0.6889$ is %s Gpc ($H_\\Lambda=%s$ km s$^{-1}$ Mpc$^{-1}$). The gravitational effect of a wall between the vacua is set "
          "by its tension $\\sigma$. A scalar whose potential scale is the dark-energy scale, $V_f^{1/4}\\approx 2.2$ meV, has a wall tension of the order of $V_f^{3/4}$, that is of order 10 meV$^3$; a bounce calculation for a "
          "potential of this kind gives 19 meV$^3$ and a thickness of 0.24 mm, and the conclusions depend only on the wall being a probe, not on these values. In units of the de Sitter radius this is "
          "$\\hat\\sigma=\\sigma R_{dS}/M_{\\rm Pl}^2$ = %s and $\\delta/R_{dS}$ = %s, so that the wall's own gravity, measured by $3\\hat\\sigma^2/4$ for $q=1$, is %s. The wall is a probe, and its thickness is invisible to any "
          "gravitational observable of the size of the de Sitter radius. That fact organizes the paper: the thin-wall treatment of Section 4 is controlled by these two numbers, and the numerical experiments of Section 5 can only "
          "reach thick, heavy walls and therefore test the structure of the solution, not the physical values." % (n(R_dS, 2), n(TW["H_Lambda_km_s_Mpc"], 1), sci(SIGMA, 1), sci(THICK, 1), sci(GRAVP, 1))))
add(("h2", "2.1 Inputs and what is derived here"))
add(("p", "The inputs taken from the benchmark dilution-triggered model are the false-vacuum density $V_f$ (the observed dark energy), the depth $q$ of the true vacuum (baseline 1, varied from $-0.5$ to 2 in Section 5.5) "
          "and the order of magnitude of the wall tension. Everything else is derived or computed in this paper: the energy-condition argument (Section 3), the junction analysis (Section 4), the numerical solver and its "
          "experiments (Section 5), and the size, mass and entropy accounting (Sections 6 and 7); none of it relies on a result of the benchmark model's cosmology. One limit of scope is stated at once, because it is the largest gap between the motivating model and the computed one: the numerical experiments are a scalar-only gravitational benchmark. "
          "They evolve the scalar field alone in a prescribed double-well potential, with the false vacuum held by the potential and not by dark matter, so they study the gravitational problem of a false-vacuum region bounded by a wall and not the "
          "matter-triggered process. What depends on the dark matter and what does not is as follows. The null-energy proposition (Section 3) is proved for the scalar together with nonnegative-density dust and does not use the numerical solutions. "
          "The thin-wall criterion (Section 4) depends on the matter only through the exterior mass parameter, to which dust contributes its conserved mass, and through the vacuum energies. The numerical persistence and creation results (Section 5) are results for the scalar-only system: "
          "whether positive-energy dust inside the core or the wall can change the classification of any run was not tested, and a controlled subset with dust, which needs a coupled solver, is the most consequential extension. The claim that the whole process has been simulated is not made."))
add(("h2", "2.2 The exterior vacuum"))
add(("p", "The thin-wall model of Section 4 and the experiments of Section 5 place a false-vacuum region (positive density $V_f$) inside an exterior of lower vacuum energy. That exterior is not the late-time de Sitter state "
          "of the parent universe: the parent's heat-death state is the false vacuum itself. The exterior is the region in which the field has rolled to the true vacuum, which for $q>0$ is anti-de Sitter and recollapsing "
          "and for $q\\le0$ is not. A false-vacuum region inside it arises when a region returns to the false vacuum, as when its matter density rises above $\\rho_{\\rm crit}$ and pins the field, while its surroundings have "
          "already made the transition; the nucleation of a bubble of the true vacuum inside the false vacuum [@Coleman1980] is the opposite configuration, whose interior is anti-de Sitter and ends in a crunch "
          "[@Abbott1985], so that it settles outcome A for that direction without further work. How such a region is reached, and whether a negative exterior is cosmologically viable, are not decided by this "
          "paper. What it can do is test how much the conclusions depend on the sign of the exterior vacuum energy: the thin-wall criterion holds for any $q>-1$ (Section 4), and the numerical experiments are repeated "
          "for a Minkowski exterior and a de Sitter exterior of lower energy in Section 5.5."))
add(("p", "Table 1 lists the scales."))
add(("table", "Table 1. Scales of the benchmark model ($q=1$, Planck 2018 cosmology).",
     ["Quantity", "Symbol", "Value"],
     [["de Sitter radius of the false vacuum", "$R_{dS}=H_\\Lambda^{-1}$", "%s Gpc (%s m)" % (n(R_dS, 2), sci(TW["R_dS_m"], 2))],
      ["Hubble rate of the false vacuum", "$H_\\Lambda$", "%s km s$^{-1}$ Mpc$^{-1}$" % n(TW["H_Lambda_km_s_Mpc"], 1)],
      ["Wall tension in units of $M_{\\rm Pl}^2/R_{dS}$", "$\\hat\\sigma$", sci(SIGMA, 1)],
      ["Wall thickness over $R_{dS}$", "$\\delta/R_{dS}$", sci(THICK, 1)],
      ["Gravitational parameter of the wall", "$3\\hat\\sigma^2/4$", sci(GRAVP, 1)],
      ["de Sitter entropy of the false vacuum", "$\\pi R_{dS}^2/\\ell_P^2$", sci(TW["S_dS"], 1)],
      ["Matter within one Hubble radius today", "$M_H$", "%s $M_\\odot$" % sci(TW["M_Hubble_matter_Msun"], 1)]],
     [3.4, 1.6, 2.2]))

# ------------------------------------------------------------------------------------------------ 3
add(("h1", "3. What the energy conditions allow"))
add(("p", "A spherically symmetric spacetime is written $ds^2=-f(u,v)\\,du\\,dv+r(u,v)^2d\\Omega^2$ with null coordinates $u$ (outgoing rays are $u=$ const) and $v$ (ingoing rays are $v=$ const) and areal radius $r$. The null "
          "expansions of the two families are proportional to $r_u$ and $r_v$. A sphere is normal if $r_v>0>r_u$, trapped if $r_v<0$ and $r_u<0$, and anti-trapped (expanding) if $r_v>0$ and $r_u>0$; the fourth "
          "quadrant, $r_u>0>r_v$, is the second exterior region of a bridge (a Kruskal-type extension). The Misner–Sharp mass [@Misner1964] is"))
add(("eq", "m=r\\left(1+\\frac{4r_ur_v}{f}\\right)", "3"))
add(("p", "in units where it includes the vacuum energy, $m=\\rho r^3/3$ for a vacuum of density $\\rho$, and a Schwarzschild–(anti-)de Sitter exterior of mass parameter $m_{\\rm ext}=2GM$ has $m=m_{\\rm ext}+\\rho_{\\rm out}r^3/3$. "
          "Einstein's equations for the model reduce to (derived symbolically in the repository)"))
add(("eq", "r_{vv}-\\left(\\ln f\\right)_{v}r_v=-\\frac{r}{2}T^{\\rm dust}_{vv}-\\frac{r}{2}\\left(\\phi_v\\right)^2", "4"))
add(("eq", "r_{uv}=-\\frac{r_ur_v}{r}-\\frac{f}{4r}+\\frac{Ufr}{4}", "5"))
add(("p", "where $U$ is the total potential (the dust term is omitted in the numerical experiments, Section 5, where the field alone is evolved) and Eq. (4) is the constraint along $v$ (and the analog along $u$). "
          "The sign of the right-hand side of Eq. (4) is the whole story about throats. For every term of the model the null energy condition holds. For a null vector $k$, "
          "$T_{\\mu\\nu}k^\\mu k^\\nu=(k^\\mu\\partial_\\mu\\phi)^2+\\rho_{\\rm dm}(k\\cdot u)^2\\ge0$: the potential terms $-g_{\\mu\\nu}V$ drop out because $k$ is null, whatever the sign of $V$ (including a negative "
          "true-vacuum density), and the dark matter is dust of nonnegative density in the Einstein frame."))
add(("p", "Proposition 1. In spherical symmetry with $T_{\\mu\\nu}k^\\mu k^\\nu\\ge0$, no null geodesic has a strict local minimum of the areal radius. Proof. A minimum along an outgoing ray requires "
          "$r_v=0$ and $r_{vv}\\ge0$ there. Eq. (4) at $r_v=0$ gives $r_{vv}=-(r/2)T_{vv}\\le0$, so $r_{vv}\\ge0$ forces $T_{vv}=0$ and the minimum is not strict; the same argument with $u$ and $v$ exchanged "
          "applies to ingoing rays. □ This is the local content of the focusing theorem. It means that a throat in the sense of Morris and Thorne, a minimal sphere that a traveler can cross, "
          "or a timelike marginal surface of Hayward's dynamic definition, needs matter that violates the null energy condition [@Hochberg1997; @Hayward1999; @Hayward1994]. "
          "The model has none, so an ordinary traversable throat cannot form in it for any parameter values or any wall profile. The statement is limited to this classical model of a canonical scalar and nonnegative-density dust: it does not exclude nontraversable bridges, quantum or semiclassical continuations that violate the null energy condition locally, or baby-universe scenarios that rest on them. Outcome B survives only as a nontraversable bridge: a minimal sphere on a spacelike slice "
          "(as in the Einstein–Rosen bridge, where the null energy condition holds trivially) that is hidden behind horizons, and outcome C only as an expanding region attached to such a bridge."))
add(("p", "That settles the central obstacle identified at the outset, and it narrows the question to what such a bridge needs. Because the question does not involve exotic matter, it can be put to the ordinary "
          "Einstein–scalar system within a thin-wall treatment and numerically with a thick wall."))

# ------------------------------------------------------------------------------------------------ 4
add(("h1", "4. Thin-wall criterion"))
add(("p", "A thin shell of tension $\\sigma$ at areal radius $R(\\tau)$, with $U=dR/d\\tau$, separates an interior de Sitter region of density $\\rho_{in}=V_f$ from an exterior Schwarzschild–(anti-)de Sitter region of mass "
          "parameter $m$ and density $\\rho_{out}=-q\\rho_{in}$ [@Israel1966; @Blau1987]. With $\\Delta=\\rho_{in}-\\rho_{out}=(1+q)V_f$ the Israel junction condition takes the form"))
add(("eq", "\\Gamma_{in}-\\Gamma_{out}=\\frac{\\sigma R}{2},\\quad \\Gamma_{in}^2=1+U^2-\\frac{R^2}{R_{dS}^2},\\quad \\Gamma_{out}^2=1+U^2-\\frac{m}{R}-\\frac{\\rho_{out}R^2}{3}", "6"))
add(("p", "where $\\Gamma=dR/dn$ is the derivative of the areal radius along the unit normal, which points from the interior to the exterior, so that a negative $\\Gamma_{in}$ means that the interior region extends beyond the equator of its de Sitter 3-sphere. "
          "Eliminating $\\Gamma_{out}$ gives"))
add(("eq", "\\Gamma_{in}=\\frac{\\sigma R}{4}+\\frac{m}{\\sigma R^2}-\\frac{\\Delta R}{3\\sigma}", "7"))
add(("p", "Motion is possible only where $U^2\\ge0$. For $R<R_{dS}$ this reads $m\\ge m_+(R)$ or $m\\le m_-(R)$ with"))
add(("eq", "m_{\\pm}(R)=\\left(\\frac{\\Delta}{3}-\\frac{\\sigma^2}{4}\\right)R^3\\pm\\sigma R^2\\sqrt{1-R^2/R_{dS}^2}", "8"))
add(("p", "and for $R\\ge R_{dS}$ it is always satisfied. The curves $m_\\pm$ bound a band of forbidden radii (Fig. 1). The band reaches up to the finite-tension critical mass"))
add(("eq", "m_c(\\hat\\sigma)=\\max_{R\\le R_{dS}}m_+(R)=(1+q)R_{dS}-\\hat\\sigma^2\\left[\\frac14-\\frac{1}{6(1+q)}\\right]R_{dS}+O(\\hat\\sigma^4)", "9"))
add(("p", "which is exact within the junction equations (the maximum is found numerically; the expansion is that of the maximum of Eq. (8) near $R=R_{dS}$ and is checked against it in Table 2). The simple form "
          "$m_c=(1+q)R_{dS}$ is the limit of vanishing tension, not an additional result: its relative error is of order $\\hat\\sigma^2$ (%s for $q=1$ and the model's tension), and the two are indistinguishable there. "
          "For $m<m_c$ every wall lies either at $R\\le R_1(m)$ or at $R\\ge R_2(m)$, with a forbidden window $R_1<R<R_2$ around $R_M=(3m/\\Delta)^{1/3}$ that no trajectory crosses. For $m>m_c$ there is no window. "
          "A numerical scan of $R$ confirms this at nine combinations of tension and depth: a window exists at $m=m_c(1-10^{-3})$ and none at $m=m_c(1+10^{-3})$ in every case (Table S0 of the supplementary information)."
          % sci(SIGMA ** 2 / 12.0, 1)))
add(("table", "Table 2. Finite-tension critical mass for $q=1$ (units of $R_{dS}$): the numerically maximized $m_c(\\hat\\sigma)$, its small-tension expansion of Eq. (9), the radius of the exterior horizon at $m_c$ and the ratio "
              "of the de Sitter and black-hole horizon entropies, $(R_{dS}/r_h)^2$, which equals 1 only in the limit.",
     ["$\\hat\\sigma$", "$m_c(\\hat\\sigma)$", "Expansion, Eq. (9)", "$r_h(m_c)$", "$S_{dS}/S_{BH}$"],
     [[n(f["sigma"], 2), n(f["m_c_exact"], 5), n(f["m_c_series"], 5), n(f["r_h"], 5), n(f["entropy_ratio_dS_over_BH"], 5)] for f in FT if f["sigma"] > 0] +
     [["%s (model)" % sci(SIGMA, 1), "2 − %s" % sci(SIGMA ** 2 / 6.0, 1), "2 − %s" % sci(SIGMA ** 2 / 6.0, 1), "1 − %s" % sci(SIGMA ** 2 / 24.0, 1), "1 + %s" % sci(SIGMA ** 2 / 12.0, 1)]],
     [1.5, 1.4, 1.4, 1.3, 1.3]))
add(("fig", "Fig1_thin_wall.png", "Fig. 1. Allowed and forbidden wall radii for $q=1$ and two values of the wall tension. The shaded band is forbidden; the dashed line is the finite-tension critical mass $m_c(\\hat\\sigma)$. "
                                   "For $m<m_c$ a horizontal line crosses the band, which a wall cannot cross; for $m>m_c$ it does not.", 6.5))
add(("p", "The consequences follow from which side of the window a wall can be on and what a regular initial slice looks like. A slice that is regular, in the sense that the interior is a de Sitter ball "
          "less than half of its 3-sphere ($\\Gamma_{in}>0$) and the exterior is an asymptotic region ($\\Gamma_{out}>0$), has $R<R_M(1-O(\\hat\\sigma))$. If $m<m_c$ it therefore lies at $R\\le R_1$, "
          "can never reach $R_2$, and collapses to $R=0$, where the singularity theorem applies [@Penrose1965]: outcome A for every configuration of subcritical mass on the branch analyzed here (spherical, thin wall, vacuum exterior of the stated form, regular initial slice), with no assumption on the initial velocity. The branch $R\\ge R_2$ needs $\\Gamma_{in}<0$, "
          "an interior beyond the equator of its own de Sitter space, which is not a regular ball; that is the configuration of the baby-universe solutions of Refs. [@Blau1987; @Farhi1990], and it cannot be "
          "reached from regular data, which is the obstacle of Ref. [@FarhiGuth1987]. For $m>m_c$ the window is absent, so a wall that starts regular can move out through $R_M$, where $\\Gamma_{in}$ changes sign."))
add(("p", "The bridge in these solutions needs to be stated precisely, because the wall is not where the minimal sphere is. On the time-symmetric slice through a turning point with $\\Gamma_{in}<0$ and $\\Gamma_{out}<0$, the "
          "areal radius is zero at the pole of the de Sitter cap, rises to $R_{dS}$ at its equator (a maximal sphere), decreases across the wall, and continues to decrease through the exterior to its minimum at the "
          "bifurcation sphere of the exterior solution, $r=r_h$, beyond which it increases again to the second asymptotic region. The minimal sphere is therefore the bifurcation sphere of the exterior's horizons, "
          "which is a marginal surface and not a throat that can be crossed; the bridge is nontraversable, as Proposition 1 requires. This construction rests on the thin-wall solution and its literature: "
          "the numerical evolution of Section 5 starts on a null cone and ends at the singularity, and does not construct such a slice, so no throat is claimed for it."))
add(("p", "Which mass is meant needs a statement. The $m$ of Eqs. (6)–(9) is the parameter $m=2GM/c^2$ of the exterior metric function $f=1-m/r-\\rho_{out}r^2/3$: for $q>0$ the Schwarzschild–anti-de Sitter mass, for $q=0$ the ADM mass, and for $q<0$ the "
          "Schwarzschild–de Sitter parameter, which is not an asymptotic mass because there is no asymptotically flat region. The quasilocal Misner–Sharp mass at the wall is $m_{\\rm MS}(R)=m+\\rho_{out}R^3/3$, so the two differ by the exterior vacuum energy. "
          "At $R=R_{dS}$ and vanishing tension the interior and exterior quasilocal masses coincide and equal $R_{dS}$ in these units for every $q$, whereas the parameter is $(1+q)R_{dS}$ (Supplementary Table S6). "
          "For $q=1$ the critical mass quoted in this paper, %s $M_\\odot$, is therefore the exterior parameter, twice the quasilocal mass at the wall (%s $M_\\odot$); statements about the mass of the region mean the parameter unless ‘quasilocal’ is written. "
          "An independent derivation (Supplement S8) reproduces Eqs. (7)–(9) from the extrinsic curvatures, the critical mass with 120-digit arithmetic (for the model's tension $m_c-2=$ %s), the horizon roots at $m_c$ and the conversion to solar masses to $10^{-16}$." %
          (sci(AUD["5_dimensional_restoration"]["q=1"]["M_c_Msun"], 2), sci(AUD["5_dimensional_restoration"]["q=1"]["M_misner_sharp_at_wall_Msun"], 2), sci(float(AUD["2_physical"]["m_c_minus_limit"]), 2))))
add(("p", "Two statements follow. First, $m>m_c(\\hat\\sigma)$ is necessary for any regular configuration to reach $\\Gamma_{in}<0$. Second, because $m_c$ is fixed by the vacuum energy through $R_{dS}$, "
          "the necessary mass is a property of the observed dark-energy scale, not of the particle physics. Whether the condition is also sufficient, and what thick walls do, is a question for the numerical evolution of the "
          "next section."))
