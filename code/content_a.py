# -*- coding: utf-8 -*-
"""Front matter, introduction, model, null-energy argument and the exact thin-wall analysis."""
from pcommon import *

R_dS = TW["R_dS_Gpc"]
M1 = [r for r in TW["mass_table"] if r["q"] == 1.0][0]
SIGMA = TW["sigma_model_hat"]
GRAVP = TW["gravity_parameter_model"]
THICK = TW["thickness_over_RdS"]
PER = VAL["per_A"]
TST = VAL["vacuum_tests"]
SC = VAL["scan"]

add(("title", TITLE))
add(("author", "Leon Sandler"))
add(("small", "Independent Researcher, Northbrook, Illinois, USA"))
add(("small", "ORCID: https://orcid.org/0009-0007-4584-808X"))
add(("small", "Corresponding author: leon.sandler@alumni.swinburne.edu"))
add(("blank",))
add(("h1", "Abstract"))
ABSTRACT = (
    "Rebirth scenarios invoke a spacetime throat carrying matter from a collapsing region into a new expanding one. We ask whether this can occur when the only matter is a canonical "
    "scalar field coupled to dark matter, whose false vacuum has the observed dark-energy density and whose true vacuum lies below it. The null energy condition holds, so no traversable throat can form: "
    "no null ray has a minimum of the areal radius. The closest allowed structure is a non-traversable bridge to a black-hole exterior. For a thin wall the Israel "
    "junction conditions give an exact critical exterior mass %s, below which every regular configuration collapses; at that mass the exterior horizon and the interior de Sitter horizon "
    "have equal area, so the threshold is the equality of their entropies. A spherically symmetric Einstein–scalar code, checked against exact solutions, "
    "forms a trapped surface in %d of %d runs, within %s (median) of the analytic horizon radius. An expanding core survives beside the trapped region only when the false-vacuum region extends beyond "
    "its own de Sitter radius. At the observed dark-energy scale that is %s Gpc, and the vacuum-energy deficit of such a region, %s solar masses, is %s times the matter inside a sphere of that radius today. "
    "Localized collapse therefore cannot produce a daughter, and matter is not transferred to it: the daughter is a pre-existing, horizon-sized false-vacuum region that its collapsing surroundings isolate. "
    "The mechanism is not supported as a cosmological reset."
    % ("$m_c=(1+q)R_{dS}$", SC["n_trapped"], SC["n_runs"], pc(abs(SC["dev_median"]), 1), n(R_dS, 1), "%s × 10$^{%d}$" % (n(M1["M_c_Msun"] / 10 ** 23, 1), 23), n(M1["deficit_over_matter_in_RdS_sphere"], 1)))
add(("p0", ABSTRACT))
add(("h1", "Keywords"))
add(("p0", "Spacetime throats; Vacuum decay; False-vacuum bubbles; Baby universes; Null energy condition; Gravitational collapse"))

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
          "universe, and the distinction is the content of the paper. Even outcome C would establish cosmological reproduction, not a cyclic return of the same universe: that additionally needs a daughter that reconnects "
          "with or replaces its parent, which no part of this work provides."))
add(("p", "The mechanism is not new in its geometry. Thin false-vacuum bubbles in general relativity were analysed in detail in the 1980s [@Blau1987; @Berezin1987], the possibility of a baby universe attached to a "
          "black hole by a throat was explored in classical and quantum form [@Frolov1989; @Farhi1990; @Fischler1990], and the obstacle that regular initial data do not lead to one was identified at the time "
          "[@FarhiGuth1987]. What is new here is the application to a specific density-triggered model with a fixed vacuum scale. The vacuum energy is the observed dark-energy density, so the critical sizes and masses "
          "are numbers rather than free parameters, and they can be compared with the matter that the universe actually contains. The paper makes four contributions."))
add(("bullet", "A null-energy argument for the whole model (Section 3): no traversable throat can form, and no null ray has a minimum of the areal radius."))
add(("bullet", "An exact thin-wall criterion (Section 4): a critical exterior mass $m_c$ below which a false-vacuum region is confined and collapses, with a threshold that coincides with the equality of the "
               "horizon entropies of the interior and the exterior (Section 6)."))
add(("bullet", "A numerical Einstein–scalar experiment with a regular centre (Section 5), validated against exact solutions, that tests the criterion with thick walls and maps the outcomes."))
add(("bullet", "An accounting of the size, mass and entropy that the criterion requires at the observed dark-energy scale, with the conclusion that matter cannot supply them (Sections 6 and 7), and falsification criteria (Section 8)."))

# ------------------------------------------------------------------------------------------------ 2
add(("h1", "2. Model and scales"))
add(("p", "Units are $8\\pi G = 1$ ($M_{\\rm Pl}$ is the reduced Planck mass) and $c = \\hbar = 1$ unless stated. The action is"))
add(("eq", "S=\\int d^4x\\sqrt{-g}\\left[\\frac{M_{\\rm Pl}^{2}}{2}R-\\frac12(\\nabla\\phi)^2-V_0(\\phi)\\right]+S_{\\rm dm}\\left[A^2(\\phi)g_{\\mu\\nu},\\psi_{\\rm dm}\\right]", "1"))
add(("p", "with $A(\\phi)=1+\\phi^2/2M^2$ the coupling to dark matter alone. For pressureless dark matter of conserved density $\\rho_*\\propto a^{-3}$ the scalar moves in $V_{\\rm eff}=V_0(\\phi)+\\rho_*A(\\phi)$. "
          "While $\\rho_*$ exceeds a critical value $\\rho_{\\rm crit}=\\mu^2M^2$ the symmetric point of $V_0=V_c-\\mu^2\\phi^2/2+\\lambda\\phi^4/4$ is a minimum and holds the field in a false vacuum of density "
          "$V_f$. When the expansion has diluted the matter below $\\rho_{\\rm crit}$ the point becomes a maximum and the field can roll to a lower vacuum. Two numbers fix the gravitational problem, and nothing else of the "
          "particle physics enters it: the false-vacuum density $V_f$, which is set equal to the observed dark-energy density, and the depth of the true vacuum, $V_t=-qV_f$ with $q>0$ (the baseline has $q=1$). "
          "The relevant length is the de Sitter radius of the false vacuum,"))
add(("eq", "R_{dS}=\\sqrt{3/V_f}=H_{\\Lambda}^{-1}", "2"))
add(("p", "which for the Planck 2018 values [@Planck2020] $H_0=67.66$ km s$^{-1}$ Mpc$^{-1}$ and $\\Omega_\\Lambda=0.6889$ is %s Gpc ($H_\\Lambda=%s$ km s$^{-1}$ Mpc$^{-1}$). The gravitational effect of a wall between the vacua is set "
          "by its tension $\\sigma$. A scalar whose potential scale is the dark-energy scale, $V_f^{1/4}\\approx 2.2$ meV, has a wall tension of the order of $V_f^{3/4}$, that is of order 10 meV$^3$; a bounce calculation for a "
          "potential of this kind gives 19 meV$^3$ and a thickness of 0.24 mm, and the conclusions depend only on the wall being a probe, not on these values. In units of the de Sitter radius this is $\\hat\\sigma=\\sigma R_{dS}/M_{\\rm Pl}^2$ = %s and $\\delta/R_{dS}$ = %s, so that the wall's own gravity, "
          "measured by $3\\hat\\sigma^2/4$ for $q=1$, is %s. The wall is a probe, and its thickness is invisible to any gravitational observable of the size of the de Sitter radius. That fact organizes the "
          "paper: the exact thin-wall treatment of Section 4 is controlled by these two numbers, and the numerical experiments of Section 5 can only reach thick, heavy walls and therefore test the structure of the "
          "solution, not the physical values." % (n(R_dS, 2), n(TW["H_Lambda_km_s_Mpc"], 1), sci(SIGMA, 1), sci(THICK, 1), sci(GRAVP, 1))))
add(("p", "Table 1 lists the scales. Two processes are of interest. The forward process is the nucleation of a bubble of true vacuum (lower density) inside the false vacuum [@Coleman1980]; its interior is anti-de Sitter and "
          "ends in a crunch [@Abbott1985], which settles outcome A for that direction without further work. The reverse process is a region of false vacuum (higher density) inside a true-vacuum exterior, as arises when a "
          "region's matter density returns above $\\rho_{\\rm crit}$ and the field is pinned in its false minimum. The false-vacuum region is then a de Sitter ball in an exterior of lower vacuum energy and "
          "positive mass $m$, which is the configuration of the false-vacuum bubbles of Refs. [@Blau1987; @Farhi1990]. The rest of the paper concerns this direction."))
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
add(("p", "where $U$ is the total potential (the dust term is omitted in the numerical experiments, Section 5, where the field alone is evolved) and Eq. (4) is the constraint along $v$ (and the analogue along $u$). "
          "The sign of the right-hand side of Eq. (4) is the whole story about throats. For every term of the model the null energy condition holds. For a null vector $k$, "
          "$T_{\\mu\\nu}k^\\mu k^\\nu=(k^\\mu\\partial_\\mu\\phi)^2+\\rho_{\\rm dm}(k\\cdot u)^2\\ge0$: the potential terms $-g_{\\mu\\nu}V$ drop out because $k$ is null, whatever the sign of $V$ (including the negative "
          "true-vacuum density), and the dark matter is dust of non-negative density in the Einstein frame."))
add(("p", "Proposition 1. In spherical symmetry with $T_{\\mu\\nu}k^\\mu k^\\nu\\ge0$, no null geodesic has a strict local minimum of the areal radius. Proof. A minimum along an outgoing ray requires "
          "$r_v=0$ and $r_{vv}\\ge0$ there. Eq. (4) at $r_v=0$ gives $r_{vv}=-(r/2)T_{vv}\\le0$, so $r_{vv}\\ge0$ forces $T_{vv}=0$ and the minimum is not strict; the same argument with $u$ and $v$ exchanged "
          "applies to ingoing rays. $\\square$ This is the local content of the focusing theorem. It means that a throat in the sense of Morris and Thorne, a minimal sphere that a traveller can cross, "
          "or a timelike marginal surface of Hayward's dynamic definition, needs matter that violates the null energy condition [@Hochberg1997; @Hayward1999; @Hayward1994]. "
          "The model has none, so a traversable throat cannot form in it for any parameter values or any wall profile. Outcome B survives only as a non-traversable bridge: a minimal sphere on a spacelike slice "
          "(as in the Einstein–Rosen bridge, where the null energy condition holds trivially) that is hidden behind horizons, and outcome C only as an expanding region attached to such a bridge."))
add(("p", "That settles the central obstacle identified at the outset, and it narrows the question to what such a bridge needs. Because the question does not involve exotic matter, it can be put to the ordinary "
          "Einstein–scalar system exactly (thin wall) and numerically (thick wall)."))

# ------------------------------------------------------------------------------------------------ 4
add(("h1", "4. Exact thin-wall criterion"))
add(("p", "A thin shell of tension $\\sigma$ at areal radius $R(\\tau)$, with $U=dR/d\\tau$, separates an interior de Sitter region of density $\\rho_{in}=V_f$ from an exterior Schwarzschild–(anti-)de Sitter region of mass "
          "parameter $m$ and density $\\rho_{out}=-q\\rho_{in}$ [@Israel1966; @Blau1987]. With $\\Delta=\\rho_{in}-\\rho_{out}=(1+q)V_f$ the Israel junction condition takes the form"))
add(("eq", "\\Gamma_{in}-\\Gamma_{out}=\\frac{\\sigma R}{2},\\quad \\Gamma_{in}^2=1+U^2-\\frac{R^2}{R_{dS}^2},\\quad \\Gamma_{out}^2=1+U^2-\\frac{m}{R}-\\frac{\\rho_{out}R^2}{3}", "6"))
add(("p", "where $\\Gamma=dR/dn$ is the derivative of the areal radius along the unit normal, which points from the interior to the exterior, so that a negative $\\Gamma_{in}$ means that the interior region extends beyond the equator of its de Sitter 3-sphere. "
          "Eliminating $\\Gamma_{out}$ gives"))
add(("eq", "\\Gamma_{in}=\\frac{\\sigma R}{4}+\\frac{m}{\\sigma R^2}-\\frac{\\Delta R}{3\\sigma}", "7"))
add(("p", "Motion is possible only where $U^2\\ge0$. For $R<R_{dS}$ this reads $m\\ge m_+(R)$ or $m\\le m_-(R)$ with"))
add(("eq", "m_{\\pm}(R)=\\left(\\frac{\\Delta}{3}-\\frac{\\sigma^2}{4}\\right)R^3\\pm\\sigma R^2\\sqrt{1-R^2/R_{dS}^2}", "8"))
add(("p", "and for $R\\ge R_{dS}$ it is always satisfied. The curves $m_\\pm$ bound a band of forbidden radii (Fig. 1). The band reaches up to the exact critical mass $m_c=\\max_{R\\le R_{dS}}m_+(R)$, "
          "which equals $(1+q)R_{dS}\\,[1+O(\\hat\\sigma^2)]$. For $m<m_c$ every wall lies either at $R\\le R_1(m)$ or at $R\\ge R_2(m)$, with a forbidden window $R_1<R<R_2$ around "
          "$R_M=(3m/\\Delta)^{1/3}$ that no trajectory crosses. For $m>m_c$ there is no window. A numerical scan of $R$ confirms this at nine combinations of tension and depth: a window exists at $m=m_c(1-10^{-3})$ "
          "and none at $m=m_c(1+10^{-3})$ in every case (Table S0 of the supplementary information). For the model's tension the correction is of order $10^{-60}$, so"))
add(("eq", "m_c=(1+q)\\,R_{dS}", "9"))
add(("p", "to all practical purposes: $m_c=%s\\,R_{dS}$ for $q=1$." % n(M1["m_c_over_RdS"], 0)))
add(("fig", "Fig1_thin_wall.png", "Fig. 1. Allowed and forbidden wall radii for $q=1$ and two values of the wall tension. The shaded band is forbidden; the dashed line is the exact critical mass $m_c$. "
                                   "For $m<m_c$ a horizontal line crosses the band, which a wall cannot cross; for $m>m_c$ it does not.", 6.5))
add(("p", "The consequences follow from which side of the window a wall can be on and what a regular initial slice looks like. A slice that is regular, in the sense that the interior is a de Sitter ball "
          "less than half of its 3-sphere ($\\Gamma_{in}>0$) and the exterior is an asymptotic region ($\\Gamma_{out}>0$), has $R<R_M(1-O(\\hat\\sigma))$. If $m<m_c$ it therefore lies at $R\\le R_1$, "
          "can never reach $R_2$, and collapses to $R=0$, where the singularity theorem applies [@Penrose1965]: outcome A for every regular configuration of subcritical mass, with no assumption on the initial velocity. The branch $R\\ge R_2$ needs $\\Gamma_{in}<0$, "
          "an interior beyond the equator of its own de Sitter space, which has a neck on every slice that crosses it. That is the configuration of the baby-universe solutions of Refs. [@Blau1987; @Farhi1990], and it cannot be "
          "reached from regular data; this is the obstacle of Ref. [@FarhiGuth1987]. For $m>m_c$ the window is absent, so a wall that starts regular can move out through $R_M$, where $\\Gamma_{in}$ changes sign. The solution then "
          "develops a neck at the wall: a minimal sphere on a spacelike slice, not forbidden by Proposition 1, whose two sides are the expanding interior and the far side of the exterior bridge."))
add(("p", "Two statements are therefore exact in the thin-wall model. First, $m>m_c$ is necessary for any regular configuration to develop a neck. Second, because $m_c=(1+q)R_{dS}$ is fixed by the vacuum energy, "
          "the necessary mass is a property of the observed dark-energy scale, not of the particle physics. Whether the condition is also sufficient, and what thick walls do, is a question for the numerical evolution of the "
          "next section."))
