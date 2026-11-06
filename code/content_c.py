# -*- coding: utf-8 -*-
"""Sections 6-10: matter and entropy accounting, assembly, phase diagram and falsification, discussion, conclusions, declarations."""
from pcommon import *

M1 = [r for r in TW["mass_table"] if r["q"] == 1.0][0]
R_dS = TW["R_dS_Gpc"]
MT = TW["mass_table"]
SC = VAL["scan"]

add(("h1", "6. Matter transfer and entropy"))
add(("p", "The plan to carry all of a universe's matter through a throat to a new beginning needs an accounting that the geometry can check. In a dynamical spacetime a flux of energy through a "
          "surface is not a conserved mass unless a time-flow vector and a hypersurface are specified. The unambiguous quantity in spherical symmetry is the Misner–Sharp mass of Eq. (3), which is conserved in "
          "vacuum and changes only by the flux of matter, and the exterior mass parameter $m_{\\rm ext}=2GM$ that it defines. The critical mass $m_c=(1+q)R_{dS}$ has a simple meaning: it is "
          "$m_c=\\Delta R_{dS}^3/3$, that is the vacuum-energy difference $\\Delta=(1+q)V_f$ integrated over a de Sitter volume of radius $R_{dS}$, $M_c=(4\\pi/3)\\Delta R_{dS}^3$. A false-vacuum region of radius "
          "$R_w$ has an exterior mass of at least $\\Delta R_w^3/3$ because of its vacuum energy alone, so that $m>m_c$ is automatic for $R_w>R_{dS}$ and no matter is needed to supply it. In the configurations analysed the daughter is vacuum "
          "energy and nothing else: none of the parent's matter enters it, and none has to. Matter is transferred in no sense, and the exterior black hole is what is left of the parent's vacuum-energy deficit."))
add(("p", "The horizon of the exterior at the critical mass has a remarkable property. The exterior has $f(r)=1-m/r+q\\,r^2/R_{dS}^2$ (for $\\rho_{out}=-q\\rho_{in}$ the anti-de Sitter radius satisfies $L^2=R_{dS}^2/q$), and "
          "its horizon obeys $r_h-m+q\\,r_h^3/R_{dS}^2=0$. At $m=(1+q)R_{dS}$ the root is exactly"))
add(("eq", "r_h=R_{dS}\\qquad\\Longrightarrow\\qquad \\frac{S_{dS}}{S_{BH}}=\\left(\\frac{R_{dS}}{r_h}\\right)^{2}=1", "11"))
add(("p", "for every $q$ (the left side of Eq. (11) uses the Gibbons–Hawking entropy of the de Sitter horizon and the Bekenstein–Hawking entropy of the black hole, each $A/4G$ [@Bekenstein1973; @Gibbons1977]; "
          "the thin-wall correction is of order $\\hat\\sigma^2$). At the critical mass the horizon entropies of the interior and of the exterior black hole are equal, and for $m>m_c$ the black hole has the larger entropy. "
          "Below $m_c$ the de Sitter horizon of a daughter would carry more entropy than the black hole that contains it, and that is the regime in which the wall is confined. The geometry therefore "
          "agrees with the generalized second law at the threshold, and the entropy bound [@Bousso2002] gives a reason for the criterion that does not depend on the details of the wall. "
          "The daughter's entropy is $S_{dS}$ = %s in units of $k_B$, the entropy of the observed dark-energy horizon." % sci(TW["S_dS"], 1)))
add(("h1", "7. How large must the region be, and where can it come from?"))
add(("p", "In the benchmark model with $q=1$ the exterior mass must exceed $m_c=%s\\,R_{dS}$, that is a mass $M_c$ = %s $M_\\odot$, for a region of radius $R_{dS}$ = %s Gpc. "
          "This mass is %s times the matter within one Hubble radius today (%s $M_\\odot$) and %s times the matter within a sphere of radius $R_{dS}$ (%s $M_\\odot$), so that today's matter provides only %s of it. "
          "The remainder is the vacuum-energy deficit, and no collection of matter within that sphere can make up the difference. In the past the matter in such a sphere was larger, and equalled $M_c$ at $a=%s$ (Fig. 5b); "
          "but then the interior was dominated by matter and contracting, which has no expanding core, and with dust inside the threshold radius becomes the Hubble radius of the denser interior, a case that "
          "the thin-wall and numerical analyses of this paper did not treat. Table 3 gives the requirement for other depths of the true vacuum; deeper vacua need more." %
          (n(M1["m_c_over_RdS"], 0), sci(M1["M_c_Msun"], 1), n(R_dS, 1), n(M1["M_c_over_matter_in_Hubble_sphere"], 1), sci(TW["M_Hubble_matter_Msun"], 1),
           n(M1["deficit_over_matter_in_RdS_sphere"], 1), sci(TW["M_matter_in_RdS_sphere_Msun"], 1), pc(M1["matter_fraction"], 0), n(TW["a_matter_equals_deficit"], 2))))
add(("table", "Table 3. Critical region for a daughter, as a function of the depth of the true vacuum ($V_t=-qV_f$), with Planck 2018 parameters. $M_H$ is the matter within one Hubble radius today and $M_{R}$ the matter "
              "within a sphere of radius $R_{dS}$.",
     ["$q$", "$m_c/R_{dS}$", "$M_c$ ($M_\\odot$)", "$M_c/M_H$", "$M_c/M_{R}$", "Horizon radius $r_h/R_{dS}$ at $m_c$"],
     [[n(r["q"], 2), n(r["m_c_over_RdS"], 2), sci(r["M_c_Msun"], 2), n(r["M_c_over_matter_in_Hubble_sphere"], 1), n(r["deficit_over_matter_in_RdS_sphere"], 2), n(r["r_h_over_RdS_at_mc"], 3)] for r in MT],
     [0.6, 1.1, 1.4, 1.0, 1.1, 1.9]))
add(("fig", "Fig5_assembly.png", "Fig. 5. (a) Critical mass of the false-vacuum region, in units of the matter within one Hubble radius today (red) and within a sphere of radius $R_{dS}$ (blue), against the depth $q$ of "
                                 "the true vacuum. (b) Matter density of the parent, relative to the vacuum-energy difference $(1+q)V_f$ for $q=1$, against its scale factor.", 6.5))
add(("p", "Three conclusions follow, and each can be wrong in a specific way. First, the size: the region has the radius of the de Sitter horizon of the observed dark-energy density, %s Gpc, %s times the Hubble radius "
          "today; no lump, void or halo is of that size. A false-vacuum region smaller than $R_{dS}$, whatever its mass, is outcome A: this is what the numerical runs with $R_w\\lesssim R_{dS}$ show and, "
          "for $m<m_c$, what the thin-wall analysis proves. Second, the origin: what is left is a region of undecayed false vacuum that is itself larger than $R_{dS}$, which in this model is a patch of the parent universe "
          "that has not made the transition. Such a patch is not created by the collapse; it exists beforehand, and what the collapse of its lower-energy surroundings does is trap it behind a horizon. That does not "
          "produce a reset of the parent: it separates a piece of the parent's own vacuum. Third, causality: while the parent accelerates, its event horizon is of order $R_{dS}$, and a region of "
          "radius $R_{dS}$ cannot be assembled by dynamics inside it, because matter beyond the horizon is out of reach; the requirement can be met only by an acausal initial condition or by a patch that survives. "
          "The special case in which the whole parent returns to the false vacuum in a contracting phase is not a localized collapse: it is the homogeneous restoration that the model assumes elsewhere, and it "
          "needs no throat." % (n(R_dS, 1), n(R_dS / 4.43, 2))))

add(("h1", "8. Phase diagram and falsification criteria"))
add(("p", "Collecting the results, the model has three regimes for a false-vacuum region of core radius $R_w$ and exterior mass $m$, in units of $R_{dS}$ and with $q=1$:"))
add(("bullet", "$m<m_c$ (exact): the wall is confined and collapses; a trapped surface forms (outcome A). This includes every region that is smaller than $R_{dS}$ and has no more than its own vacuum energy."))
add(("bullet", "$m>m_c$ and $R_w\\lesssim R_{dS}$: the numerical runs show a trapped surface and a normal core (outcome A): extra mass by itself does not keep the core from collapsing."))
add(("bullet", "$R_w\\gtrsim R_{dS}$ (which implies $m>m_c$): a trapped surface forms next to a persisting expanding core, with a neck-like minimum of the areal radius on a spacelike slice; this is outcome C, a daughter behind a horizon. "
               "It requires the horizon-sized region of Section 7, and it does not transfer matter."))
add(("p", "Each statement of the paper can be falsified, and the falsification criteria are:"))
add(("num", "A violation of the null energy condition by the matter that the model contains. This would remove Proposition 1 and reopen outcome B. For the canonical scalar and dust of the model it is excluded; for "
            "a model with a non-canonical kinetic term or a ghost it must be checked afresh."))
add(("num", "A numerical solution, with a thick or thin wall, in which a trapped region and an expanding core coexist at $m_{\\rm ext}<m_c$ with regular initial data. The thin-wall result says that none exists; "
            "a counterexample would show that the junction analysis misses a configuration."))
add(("num", "A mechanism that forms a vacuum-dominated false region of radius $\\gtrsim R_{dS}$ by the collapse of smaller regions. None is known, and the numerical runs with $R_w\\lesssim R_{dS}$ give none."))
add(("num", "A daughter whose horizon entropy exceeds that of the exterior black hole. This would conflict with the generalized second law at the threshold of Eq. (11), and either the law or the classical "
            "solution would have to give way."))
add(("num", "A regular, global evolution of the thick-wall problem through the horizon that shows the core is not a separate daughter, for example by closing the anti-trapped region. The anti-trapped region "
            "cannot be undone by the dynamics in the numerical runs, but a global argument has not been given."))

add(("h1", "9. Discussion"))
add(("p", "The result is a conditional negative with one specific positive. The negative: localized collapse in this model does not create a new universe, matter is not carried through a throat, and "
          "the null energy condition removes traversable throats outright. The positive: a regular expanding de Sitter core next to a trapped region does exist when the false-vacuum region "
          "extends beyond its own de Sitter radius, and its existence is sharply characterized by $m>m_c=(1+q)R_{dS}$ together with $R_w\\gtrsim R_{dS}$, a condition that equates the horizon entropies of "
          "the interior and the exterior black hole. What the positive result describes is the isolation of an existing region, not a rebirth."))
add(("p", "Three caveats limit even that. The daughter obtained here is a cold de Sitter region with the vacuum energy $V_f$. Producing a hot big bang inside it needs a reheating "
          "process, and none is included. The daughter also does not reconnect with the parent: to call the process a cycle one would have to show that the daughter replaces or re-enters the parent, which "
          "the geometry does not suggest, because the daughter is behind the horizon of an exterior black hole. And the daughter exists together with an exterior black hole whose mass is at least $M_c$, "
          "so that a parent which isolated such a region would be left with an object of that mass and would differ from the universe observed."))
add(("p", "The work also leaves things undone. The numerical solution is spherically symmetric and classical, the walls are not as thin as the physical wall, the evolution ends at the singularity, and the "
          "dark matter is omitted from the numerical experiments (it is present in the null-energy argument and the thin-wall criterion, which do not depend on it). The tunnelling route to a baby universe "
          "[@Farhi1990; @Fischler1990] is outside the classical treatment, and a quantum or semiclassical stress-energy that violates the null energy condition locally could in principle allow a throat; "
          "the present model has neither. Perturbations of the exterior and of the wall are not studied, and a nonspherical collapse could form different horizons. "
          "A global numerical evolution through the horizon, with excision of the singularity, would settle whether the core is a separate region, and is the most useful extension."))

add(("h1", "10. Conclusions"))
add(("p", "The question was whether a matter-dilution-triggered vacuum transition can produce a non-singular gravitational bridge to an expanding daughter universe, or whether the energy conditions and causal "
          "structure of the model forbid it. For this model the answer has three parts. A traversable throat is forbidden: the null energy condition holds, and no null ray has a minimum of the areal radius. "
          "A non-traversable bridge with an expanding daughter is allowed exactly when the exterior mass exceeds the critical value $m_c=(1+q)R_{dS}$, where the horizon entropies of the interior "
          "and exterior coincide, and in the numerical runs when the false-vacuum core extends beyond its own de Sitter radius. At the observed dark-energy scale that is a region of %s Gpc whose "
          "vacuum-energy deficit, %s $M_\\odot$, is %s times the matter within a sphere of that radius today. A localized collapse cannot supply it, and the matter of the parent does not pass to the daughter. "
          "The mechanism is therefore best described as the isolation of a horizon-sized false-vacuum patch, not as a cyclic reset, and it is not supported as an end-of-universe process of the model." %
          (n(R_dS, 1), sci(M1["M_c_Msun"], 1), n(M1["deficit_over_matter_in_RdS_sphere"], 1))))
add(("h1", "Declarations"))
add(("p0", "Funding: this work received no external funding."))
add(("p0", "Competing interests: the author declares no competing interests."))
add(("p0", "Author contributions (CRediT): Leon Sandler, sole author: conceptualization, methodology, software, formal analysis, writing, review and editing."))
add(("p0", "Data and code availability: the code, the data of every figure and table and the manuscript are archived at https://doi.org/%s (software, all versions) and https://doi.org/%s "
           "(manuscript, all versions), and are available at %s. One command per script regenerates every number and figure; the results quoted in the text are read from the JSON files that the scripts write." % (DOIS["code_concept"], DOIS["paper_concept"], GITHUB)))
add(("p0", "Use of generative AI: a generative AI tool (Claude, Anthropic) assisted with code development, numerical experiments, drafting and editing. The author checked all of it and takes full responsibility for the content."))
add(("refs",))

# ---------------------------------------------------------------------------------------------- supplement
adds(("title", "Supplementary Information"))
adds(("small", "%s" % TITLE))
adds(("small", "Leon Sandler, Independent Researcher, Northbrook, Illinois, USA; ORCID 0009-0007-4584-808X"))
adds(("blank",))
adds(("h1", "S1. Equations and scheme"))
adds(("p", "The metric is $ds^2=-f\\,du\\,dv+r^2d\\Omega^2$ with $f=e^{\\sigma}$, $G_{\\mu\\nu}=T_{\\mu\\nu}$ ($8\\pi G=1$) and $T_{\\mu\\nu}=\\partial_\\mu\\phi\\partial_\\nu\\phi-g_{\\mu\\nu}[(\\nabla\\phi)^2/2+U]$. "
           "A symbolic computation of the Einstein tensor (derive_equations.py) gives"))
adds(("eq", "r_{uu}-\\sigma_ur_u=-\\frac{r}{2}\\phi_u^2,\\qquad r_{uv}=-\\frac{r_ur_v}{r}-\\frac{f}{4r}+\\frac{Ufr}{4}", "S1"))
adds(("eq", "\\phi_{uv}=-\\frac{\\phi_ur_v+\\phi_vr_u}{r}-\\frac{f}{4}U'(\\phi),\\qquad \\sigma_{uv}=\\frac{f+4r_ur_v}{2r^2}-\\phi_u\\phi_v", "S2"))
adds(("p", "The first of Eqs. (S1) is the constraint along $u$ (that along $v$ is Eq. (4) of the paper) and the others are the evolution equations. The scheme advances rows of constant $u$ on a uniform lattice $u_j=jh$, $v_i=ih$. "
           "A cell with corners $W=(j,i-1)$, $A=(j+1,i-1)$, $B=(j,i)$ gives the new point $C=(j+1,i)$ by $X_C=X_A+X_B-X_W+h^2X_{uv}$, with $X_{uv}$ evaluated at the cell centre from the mean values and the central differences, "
           "iterated three times (predictor–corrector). Near the centre the evolution of $\\sigma$ divides by a small radius and its rounding noise is amplified, so within $K=6$ points of the centre $\\sigma$ is integrated from the constraint "
           "$\\sigma_v=(r_{vv}+r\\phi_v^2/2)/r_v$ (trapezoid), and beyond it from the evolution equation. Each row is solved with the centre value of $\\sigma$ as the unknown: elementary flatness requires "
           "$\\sigma_c=\\ln(4r_v^2)$ at the centre, which defines a residual that is solved by the secant method (ten passes at most); an iteration on the plain residual has a gain of modulus greater than one and diverges. "
           "The centre value of $\\phi$ follows from $\\phi_t=2\\phi_v$ at the centre by an implicit trapezoid with the one-sided derivative (the explicit iteration has gain $-1.5$)."))
adds(("p", "The data on the initial cone are $r=v/2$, $\\phi(0,v)=A\\psi(v)$ and $\\sigma(0,v)=\\int_0^v r\\phi_v^2dv$, the integral of the constraint along $v$ with $f=1$ at the centre. The "
           "potential $U=(\\kappa/4)(\\psi^2-1)^2+s\\psi+d$ is tuned by a root finder so that the false minimum has density 3 and the true minimum $-3q$."))
adds(("h1", "S2. Numerical diagnostics"))
adds(("p", "Misner–Sharp mass of exact solutions: Table 2 of the paper. The diagnostics of each run are the signs of $r_u$ (a forward difference in $u$) and $r_v$ (a central difference in $v$) at every "
           "point more than four cells from the centre and eight from the outer boundary; the quadrant map of Fig. 4 is produced by the same classification."))
adds(("h1", "S3. Thin-wall check"))
adds(("p", "Table S0 lists the check of the critical mass. For each of nine combinations of $q$ and $\\hat\\sigma$ the exact value $m_c=\\max_Rm_+(R)$ is computed on a grid of two million points, and the existence of "
           "a forbidden window is tested at $m_c(1\\mp10^{-3})$."))
CHK = TW["critical_mass_check"]
adds(("table", "Table S0. Check of the critical mass (units of $R_{dS}$).",
      ["$q$", "$\\hat\\sigma$", "$m_c$", "Window at $m_c(1-10^{-3})$", "Window at $m_c(1+10^{-3})$"],
      [[n(c["q"], 1), n(c["sigma"], 2), n(c["m_c"], 5), "yes" if c["window_below"] else "no", "yes" if c["window_above"] else "no"] for c in CHK], [0.7, 0.9, 1.0, 1.6, 1.6]))
adds(("h1", "S4. Scan of the collapse runs"))
adds(("p", "Supplementary Table S1 lists all runs: field amplitude $A$, wall width $w$, core radius $R_w=V_w/2$, exterior mass $m_{\\rm ext}$ (standard deviation across the exterior in parentheses), $R_M=(3m_{\\rm ext}/\\Delta)^{1/3}$, "
           "the radius of the first trapped point, the horizon radius $r_h$ of the measured mass, and the number of anti-trapped points in the core on the last row."))
adds(("table", "Table S1. Collapse runs ($q=1$, $\\kappa=24$).",
      ["$A$", "$w$", "$R_w$", "$m_{\\rm ext}$", "$R_M$", "$r$ (first trapped)", "$r_h$", "Anti-trapped core points"],
      [[n(r["A"], 2), n(r["w"], 2), n(r["R_w"], 2), "%s (%s)" % (n(r["m_ext"], 3), n(r["m_ext_sd"], 3)), n(r["R_M"], 3),
        n(r["first_trapped"]["r"], 3) if r["first_trapped"] else "—", n(r["r_h_SAdS"], 3), n(r["at_points_last_row"])] for r in SCAN],
      [0.6, 0.6, 0.6, 1.4, 0.8, 1.2, 0.8, 1.4], 14))
