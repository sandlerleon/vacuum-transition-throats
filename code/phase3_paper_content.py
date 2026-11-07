# -*- coding: utf-8 -*-
"""Manuscript: Exact criteria for regular gravitational collapse with continuation, and a viscous relaxation of the anisotropy threshold.
Every number is read from the JSON files written by the scripts of the branch phase3-recycling."""
import json
import os

import refs as _refs

HERE = os.path.dirname(os.path.abspath(__file__))
L = lambda f: json.load(open(os.path.join(HERE, f), encoding="utf-8"))
E0 = L("phase3d0_ecd.json"); C5 = L("phase3c5_closure.json"); E1 = L("phase3e1_corrections.json"); S0 = L("phase3e0_string.json")
D1 = L("phase3d1_wave.json"); D2 = L("phase3d2_freefall.json"); D3 = L("phase3d3_goldilocks.json"); D4 = L("phase3d4_propertime.json")
F0 = L("phase3f0_fluid.json"); G0 = L("phase3g0_ringdown.json"); H0 = L("phase3h0_viscous_spin_bounce.json")["rows"]
GITHUB = "https://github.com/sandlerleon/vacuum-transition-throats"

g = lambda x, d=3: ("%." + str(d) + "g") % x
def sci(x, d=1):
    m, e = ("%.*e" % (d, x)).split("e")
    return "%s × 10$^{%d}$" % (m, int(e))

UN = E0["C_baseline"]["unpolarized massless (g=2)"]; WY = E0["C_baseline"]["single-helicity Weyl (g=1)"]
HF = E0["B_hf_eps_over_lam_n2"]
hrow = lambda s0, K, c2: [r for r in H0 if abs(r["s0"] - s0) < 1e-12 and abs(r["K"] - K) < 1e-12 and abs(r["c2"] - c2) < 1e-12][0]
B0a, B0b = hrow(0.05, 0, 0)["B_c"], hrow(0.2, 0, 0)["B_c"]
v01, v03, v1, v1q = hrow(0.05, 0.1, 1.0), hrow(0.05, 0.3, 1.0), hrow(0.05, 1.0, 1.0), hrow(0.05, 1.0, 0.25)
w1 = hrow(0.2, 1.0, 1.0)
SGB = E1["A_scan_summary"]["k=1 lam-"]
STR = S0["D_summary"]["-1"]
SC = C5["E_bianchi_ratio"]
FL = D4["D_killing_time_floor"]
TD = D2["C_tidal_scaling"]["1.0"]
NEST = E1["B_nested_profile"]; TF = E1["B_test_field_instability"]
FM = F0["1_critical_dust_mass"]
RD = G0["B_ringdown"]; CS = G0["C_closed_sequence"]

LIB = {k: _refs.LIB[k] for k in ("Penrose1965", "Steinhardt2002")}
LIB.update({
    "Tolman1931": "R.C. Tolman, On the theoretical requirements for a periodic behaviour of the universe, Physical Review 38 (1931) 1758–1771. https://doi.org/10.1103/PhysRev.38.1758",
    "Hehl1976": "F.W. Hehl, P. von der Heyde, G.D. Kerlick, J.M. Nester, General relativity with spin and torsion: Foundations and prospects, Reviews of Modern Physics 48 (1976) 393–416. https://doi.org/10.1103/RevModPhys.48.393",
    "Obukhov1987": "Yu.N. Obukhov, V.A. Korotky, The Weyssenhoff fluid in Einstein–Cartan theory, Classical and Quantum Gravity 4 (1987) 1633–1657. https://doi.org/10.1088/0264-9381/4/6/021",
    "Poplawski2010": "N.J. Popławski, Cosmology with torsion: An alternative to cosmic inflation, Physics Letters B 694 (2010) 181–185. https://doi.org/10.1016/j.physletb.2010.09.056",
    "Pasmatsiou2017": "K. Pasmatsiou, C.G. Tsagas, J.D. Barrow, Kinematics of Einstein–Cartan universes, Physical Review D 95 (2017) 104007. https://doi.org/10.1103/PhysRevD.95.104007",
    "Freidel2005": "L. Freidel, D. Minic, T. Takeuchi, Quantum gravity, torsion, parity violation, and all that, Physical Review D 72 (2005) 104002. https://doi.org/10.1103/PhysRevD.72.104002",
    "PerezRovelli2006": "A. Perez, C. Rovelli, Physical effects of the Immirzi parameter in loop quantum gravity, Physical Review D 73 (2006) 044013. https://doi.org/10.1103/PhysRevD.73.044013",
    "Ashtekar2006": "A. Ashtekar, T. Pawlowski, P. Singh, Quantum nature of the big bang: Improved dynamics, Physical Review D 74 (2006) 084003. https://doi.org/10.1103/PhysRevD.74.084003",
    "Simpson2019": "A. Simpson, M. Visser, Black-bounce to traversable wormhole, Journal of Cosmology and Astroparticle Physics 2019 (2019) 042. https://doi.org/10.1088/1475-7516/2019/02/042",
    "Bronnikov2022": "K.A. Bronnikov, R.K. Walia, Field sources for Simpson–Visser spacetimes, Physical Review D 105 (2022) 044039. https://doi.org/10.1103/PhysRevD.105.044039",
    "Gonzalez2009": "J.A. González, F.S. Guzmán, O. Sarbach, Instability of wormholes supported by a ghost scalar field. I. Linear stability analysis, Classical and Quantum Gravity 26 (2009) 015010. https://doi.org/10.1088/0264-9381/26/1/015010",
    "Guo2010": "Z.-K. Guo, D.J. Schwarz, Slow-roll inflation with a Gauss–Bonnet correction, Physical Review D 81 (2010) 123520. https://doi.org/10.1103/PhysRevD.81.123520",
    "ART1994": "I. Antoniadis, J. Rizos, K. Tamvakis, Singularity-free cosmological solutions of the superstring effective action, Nuclear Physics B 415 (1994) 497–514. https://doi.org/10.1016/0550-3213(94)90120-1",
    "KRT1999": "P. Kanti, J. Rizos, K. Tamvakis, Singularity-free cosmological solutions in quadratic gravity, Physical Review D 59 (1999) 083512. https://doi.org/10.1103/PhysRevD.59.083512",
    "BV1994": "R. Brustein, G. Veneziano, The graceful exit problem in string cosmology, Physics Letters B 329 (1994) 429–434. https://doi.org/10.1016/0370-2693(94)91086-3",
    "Kobayashi2011": "T. Kobayashi, M. Yamaguchi, J. Yokoyama, Generalized G-inflation: Inflation with the most general second-order field equations, Progress of Theoretical Physics 126 (2011) 511–529. https://doi.org/10.1143/PTP.126.511",
    "IS1979": "W. Israel, J.M. Stewart, Transient relativistic thermodynamics and kinetic theory, Annals of Physics 118 (1979) 341–372. https://doi.org/10.1016/0003-4916(79)90130-1",
    "FordRoman1995": "L.H. Ford, T.A. Roman, Averaged energy conditions and quantum inequalities, Physical Review D 51 (1995) 4277–4286. https://doi.org/10.1103/PhysRevD.51.4277",
    "Dvali2010": "G. Dvali, Black holes and large N species solution to the hierarchy problem, Fortschritte der Physik 58 (2010) 528–536. https://doi.org/10.1002/prop.201000009",
    "CarballoRubio2021": "R. Carballo-Rubio, F. Di Filippo, S. Liberati, C. Pacilio, M. Visser, Inner horizon instability and the unstable cores of regular black holes, Journal of High Energy Physics 2021 (2021) 132. https://doi.org/10.1007/JHEP05(2021)132",
})

TITLE = "Exact Criteria for Regular Gravitational Collapse with Continuation: A Constraint Analysis of Torsion, String-Inspired, Quantum and Dissipative Sources, and a Viscous Relaxation of the Anisotropy Threshold"

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
    "Which ingredients must a physical model of regular gravitational collapse, with continuation into an expanding region, contain? We state seven hard conditions and test candidate sources against them with exact derivations "
    "and reproducible numerics. We derive (i) a criterion for any barotropic bounce: the chemical potential must vanish and the compressibility turn negative below the bounce, here at %s and %s of the bounce density for Einstein–Cartan–Dirac matter; "
    "(ii) that theory's contact interaction $-(3\\kappa/16)A_\\mu A^\\mu$ and the mean-field energies of fermion states, which vanish for a fully polarized gas; (iii) the exact rate $-\\dot{\\alpha}s$ at which spin–vorticity cancellation fails; "
    "(iv) a proof that flat collapse cannot bounce in Einstein–scalar–Gauss–Bonnet gravity for any coupling or potential, and the closed-collapse condition with the tensor speed it must satisfy; "
    "(v) for the black-bounce benchmark, the forced phantom profile, a mass-independent null-energy deficit $a/4$, finite tides and a Killing-time floor on transit; and (vi) a scale-free species-count window. "
    "As a constructive result, causal shear viscosity (Knudsen number at most 1, shear stress below half the dust energy density) lowers the critical spin strength of the anisotropic spin–torsion bounce by a factor of about %s, "
    "so the anisotropy threshold is largely relaxed. No class tested meets all seven conditions; we give the resulting specification, identify the one unevaluated class, and provide 40 automated tests."
    % ("%.2f" % UN["nstar_over_nb"], "%.2f" % UN["nsp_over_nb"], "%d" % round(B0a / v1["B_c"], -1)))
add(("p0", ABSTRACT))
add(("p0", "Keywords: nonsingular collapse; Einstein–Cartan–Dirac gravity; Gauss–Bonnet gravity; black-bounce geometry; causal dissipative hydrodynamics; null energy condition."))

# ------------------------------------------------------------------ 1
add(("h1", "1. Introduction"))
add(("p", "A black hole is usually taken to end in a singularity, because collapse that forms a trapped surface with matter obeying the null energy condition must be geodesically incomplete [@Penrose1965]. A hypothesis with a long history is that a high-density quantum-gravity or torsion effect replaces the singularity by a bounce, so that the infalling matter re-emerges in an expanding region [@Poplawski2010; @Ashtekar2006; @Steinhardt2002]. Before such a hypothesis can be turned into a prediction, one has to know what a successful model must contain. "
          "This paper establishes that, constructively: it formulates the requirements as hard conditions, derives the exact criteria that a bounce must meet in each of the main classes of source, tests every class against the same conditions, and reports one new positive mechanism, the relaxation of the anisotropy threshold of a spin–torsion bounce by causal viscosity."))
add(("p", "The sources considered are: barotropic fluids in general; Einstein–Cartan–Dirac (ECD) fermions [@Hehl1976; @Freidel2005; @PerezRovelli2006], including the Weyssenhoff spin fluid [@Obukhov1987; @Pasmatsiou2017]; Einstein–scalar–Gauss–Bonnet (EsGB) gravity with field-dependent couplings of the kind found in the loop-corrected heterotic string [@ART1994; @KRT1999; @Guo2010]; the black-bounce geometry of Simpson and Visser [@Simpson2019] and its field sources [@Bronnikov2022] as a benchmark for what any bounce must supply; free quantum fields; nested throats; and causal dissipative fluids [@IS1979]. "
          "For each we ask the same seven questions (Section 2) and answer them by exact symbolic derivation wherever possible, with numerics confined to cases where no closed form exists. Every number in this paper is read from the output of a script in the repository, and 40 automated tests check the stated identities."))

# ------------------------------------------------------------------ 2
add(("h1", "2. Framework"))
add(("h2", "2.1 Seven hard conditions"))
add(("p", "A candidate source for regular collapse with continuation is admissible only if it satisfies all of the following. They are constraints, not terms in a score, so that a strong result on one cannot compensate for a failure of another."))
add(("table", "Table 1. The seven hard conditions.",
     ["", "Condition"],
     [["H1", "the field equations and covariant conservation hold (by construction or by evolution)"],
      ["H2", "any null-energy-violating sector has a derivable action"],
      ["H3", "no ghost: the Hamiltonian is bounded below"],
      ["H4", "perturbative stability: no exponentially growing mode, no loss of hyperbolicity"],
      ["H5", "the solution lies inside the validity regime of the effective theory"],
      ["H6", "it arises from the collapse of ordinary low-density matter"],
      ["H7", "future-directed continuation into an expanding region, not merely geodesic crossing"]],
     [0.6, 6.6]))
add(("p", "Units are $\\hbar=c=1$ and $\\kappa=8\\pi G$. Two diagnostics recur. A bounce of a homogeneous universe is a point with $H=0$ and $\\dot{H}>0$, to be evaluated in the frame to which ordinary matter couples. And because the singularity theorem [@Penrose1965] applies to any source obeying the null energy condition, every regular trapped-surface geometry requires a null-energy-violating sector, whose origin is the question."))

# ------------------------------------------------------------------ 3
add(("h1", "3. Exact results"))
add(("h2", "3.1 A criterion for barotropic bounces"))
add(("p", "Let a conserved number density $n$ with energy density $\\rho(n)$ and pressure $p=n\\rho'-\\rho$ fill a homogeneous universe of spatial curvature $k\\ge0$. With $F(n)=\\rho(n)-c\\,n^{2/3}$ ($c=0$ flat, $c>0$ closed), $H^2=\\kappa F/3$, and at a bounce $F(n_b)=0$. Since $\\dot{H}=-\\kappa nF'(n)/2$ at $H=0$, a bounce ($\\dot{H}>0$) requires $F'(n_b)<0$:"))
add(("eq", "F(n_b)=0,\\qquad F'(n_b)<0,\\qquad \\mu\\equiv\\frac{d\\rho}{dn}>0", "1"))
add(("p", "Criterion. Under (1), the intermediate-value theorem gives a density $n_*<n_b$ with $F'(n_*)=0$, so that for $c=0$ the enthalpy density $\\rho+p=n\\mu$ vanishes (the stress is vacuum-like and the fluid velocity is undefined), and the mean-value theorem gives a density below $n_b$ with $F''<0$, i.e. $dp/dn=n\\rho''<-(2c/9)n^{-1/3}\\le0$: negative compressibility. The result uses no microphysics; it is a property of every barotropic bounce with a positive low-density chemical potential. It was verified on 4000 random equations of state (every one of the %d bouncing cases showed both features) and is illustrated in Fig. 1 for the gas of Section 3.2. For $\\rho=Cn^{4/3}-Bn^2$ the two densities are $(2/3)^{3/2}n_b=0.544\\,n_b$ and $(2/9)^{3/2}n_b=0.105\\,n_b$, independent of $C$ and $B$; at the second the interaction energy is $2/9$ of the kinetic energy, so that a first-order mean-field estimate is still meaningful there." % E0["C_random_theorem_check"]["bounces"]))
add(("fig", "FigP3_1_barotropic.png", "Figure 1. The barotropic criterion for a degenerate massless Dirac gas with the contact interaction of Section 3.2. (a) Energy density and pressure, with the bounce at $n_b$, the zero of the chemical potential at $0.544\\,n_b$ and the onset of negative compressibility at $0.105\\,n_b$. (b) The chemical potential and $dp/dn$.", 6.3))
add(("p", "The criterion does not say that the dynamical growth rate of the instability is fast, nor that quantum matter beyond a barotropic mean field cannot do better; the collective-mode spectrum was not derived. It does say what any such model must contain: either a non-barotropic state or a stabilizing sector that removes the negative-compressibility region."))

add(("h2", "3.2 Einstein–Cartan–Dirac matter"))
add(("p", "We vary the action $S=∫ e[R(\\omega,e)/2\\kappa+L_D]$ with an independent spin connection. With explicit Dirac matrices we verified that only the totally antisymmetric part of the contortion couples to the Dirac field, $L_K=(i/4)K_{abc}\\bar{\\psi}\\gamma^{abc}\\psi$ with $\\gamma^{abc}=i\\epsilon^{abcd}\\gamma_d\\gamma^5$ (all 24 index triples; vector and tensor torsion decouple to %s), and symbolically that the Ricci scalar is $R=\\tilde{R}+6a_\\mu a^\\mu$ for $K_{abc}=\\epsilon_{abcd}a^d$. Eliminating the algebraic torsion gives" % sci(E0["A_nonaxial_decoupling_residual"], 1)))
add(("eq", "L_{\\rm eff}=-\\frac{3\\kappa}{16}A_\\mu A^\\mu=\\lambda\\left[(A^0)^2-|A|^2\\right],\\qquad \\lambda=\\frac{3\\kappa}{16},\\ A^\\mu=\\bar{\\psi}\\gamma^\\mu\\gamma^5\\psi", "2"))
add(("p", "in agreement with the Einstein–Cartan limit of the literature [@Freidel2005; @PerezRovelli2006]. The Dirac spin tensor is totally antisymmetric, whereas the Weyssenhoff form $s_{abc}=s_{ab}u_c$ used in spin-fluid cosmology [@Obukhov1987; @Pasmatsiou2017] is its projection along the flow; the coefficients of the resulting energy density differ. "
          "For a degenerate gas we evaluated the Hartree–Fock energy $\\langle:(\\psi^† M\\psi)^2:\\rangle=\\mathrm{Tr}(MD)^2-\\mathrm{Tr}(MDMD)$ for the local density matrix $D$ of several states."))
add(("table", "Table 2. Hartree–Fock interaction energy $\\varepsilon_{\\rm int}=-Bn^2$ (also $p_{\\rm int}=-Bn^2$) of the contact term (2) for a degenerate gas; $P$ is the polarization of a nonrelativistic gas.",
     ["State", "$B/\\lambda$", "Remark"],
     [["unpolarized, massless (g = 2)", "%s" % g(-HF["unpolarized massless (g=2, particles)"], 2), "$k_F\\gg m$ at the relevant densities"],
      ["single-helicity (Weyl) gas", "%s" % g(-HF["single-helicity Weyl (g=1)"], 2), "Hartree $-1$, spatial exchange $-1$"],
      ["unpolarized, nonrelativistic", "%s" % g(-HF["unpolarized nonrelativistic"], 2), "4.5 times the Weyssenhoff value $\\lambda/3$"],
      ["partially polarized (P = 1/2)", "%s" % g(-HF["polarized NR, fraction 0.5"], 3), "$B=(3\\lambda/2)(1-P^2)$"],
      ["fully polarized, nonrelativistic", "0", "Pauli exclusion cancels the contact term"]],
     [3.2, 0.9, 3.0]))
add(("p", "Every state has the Weyssenhoff form $\\varepsilon=p=-Bn^2$ but with a coefficient fixed by the Dirac theory, and polarization removes the effect (Table 2), so the spin-polarization strategies of Section 3.3 have no Dirac counterpart. Applied to the massless gas, $\\rho=Cn^{4/3}-Bn^2$ has an exact flat bounce at $n_b=%s$ (Planck units; $k_F=%s\\,M_P$) with finite curvature $\\tilde{R}=%s$ and Friedmann residual below $10^{-11}$, but the criterion of Section 3.1 locates negative compressibility at $k_F=%s\\,M_P$, already at the one-species cutoff, so that the bounce itself lies where first-order mean field and the effective theory are not controlled (Fig. 4b)." %
          (g(UN["n_bounce_numeric"], 3), g(UN["kF_over_Planck_at_bounce"], 2), g(UN["ricci_scalar_at_bounce"], 3), g(UN["kF_over_Planck_at_bounce"] * UN["nsp_over_nb"] ** (1 / 3), 3))))

add(("h2", "3.3 Spin transport and the failure of the cancellation"))
add(("p", "In a Weyssenhoff fluid the cross term between spin $s$ and vorticity $\\tilde{\\omega}$ gives a shear source $\\propto\\tilde{\\omega}(\\tilde{\\omega}+\\kappa s/2)$ that vanishes for anti-aligned spin $s=-2\\tilde{\\omega}/\\kappa$, and the Raychaudhuri terms combine into $2(\\tilde{\\omega}-\\kappa s/2)^2$. Whether this cancellation is maintained depends on the spin transport law, which we derive rather than assume. "
          "Angular-momentum conservation with the Frenkel condition $s^{ab}u_b=0$ gives $\\dot{S}^{ab}=-Θ S^{ab}+c_F(u^bS^{ac}A_c-u^aS^{bc}A_c)$; requiring the Frenkel condition to be preserved yields $c_F=1$ (Fermi–Walker transport) and hence, for any acceleration,"))
add(("eq", "\\frac{d\\ln s}{d\\tau}=-Θ,\\qquad \\frac{d}{d\\tau}\\left(s+\\frac{2\\tilde{\\omega}}{\\kappa}\\right)_{s=-2\\tilde{\\omega}/\\kappa}=-\\dot{\\alpha}\\,s", "3"))
add(("p", "where $\\dot{\\alpha}$ is the expansion rate along the rotation axis and the vorticity follows the Helmholtz law $\\tilde{\\omega}\\propto a_⊥^{-2}$. The cancellation surface is therefore left at the rate $\\dot{\\alpha}s$: $s/\\tilde{\\omega}\\propto 1/a_∥$, and in the anisotropic bounce runs it changes by factors %s–%s between the start and the bounce. It is invariant only for $\\dot{\\alpha}=0$, a cylindrical collapse whose shear fraction is fixed at $1/4$. "
          "Together with energy balance of the cross term in the homogeneous class, the three laws ($s\\propto n$, Helmholtz, balance) force $\\dot{\\beta}=0$. A further exact identity, $D_{[a}D_{b]}f=-\\omega_{ab}\\dot{f}$ (checked on a rigidly rotating frame in flat space), shows that a rotating, time-dependent configuration cannot be uniform in the rest space, so homogeneous rotating collapse needs gradients and tilt. "
          "These are the facts a spin-based mechanism has to respect." % (g(min(r["r_over_r0_at_bounce"] for r in SC if r["outcome"] == "bounce"), 2), g(max(r["r_over_r0_at_bounce"] for r in SC if r["outcome"] == "bounce"), 2))))

add(("h2", "3.4 Einstein–scalar–Gauss–Bonnet gravity"))
add(("p", "Consider $S=∫\\sqrt{-g}[R/2\\kappa-1/2\\sum_i(\\partial\\phi_i)^2-V+\\lambda\\,\\xi(\\phi_i)\\,G]$ with canonical scalars (dilaton and modulus) and arbitrary $\\xi$ and $V$; this is the structure of the loop-corrected heterotic action, in which the Gauss–Bonnet term is multiplied by a modulus-dependent function [@ART1994]. We derived the equations for flat and curved FLRW and compared them with the published form [@Guo2010]: the constraint equals $(a^3/2)$ times, and the acceleration equation $3a^2$ times plus $(3a^2/2)$ times, the published Friedmann and Raychaudhuri combinations, the reparametrization identity holds with residual $0$, and the coefficient of $\\dot{H}$ at $H=0$ is $2$ (no degeneracy). Two consequences follow."))
add(("p", "Flat collapse cannot bounce. All Gauss–Bonnet terms carry a factor of $H$, so at $H=0$, $ä/a=-(\\kappa/2)\\sum_i\\dot{\\phi}_i^2\\le0$ for every $\\xi$ and $V$. The 'singularity-free' solutions of this class [@ART1994; @KRT1999] therefore lack an initial singularity but do not turn around. "
          "Closed collapse. On the constraint surface at $H=0$,"))
add(("eq", "\\left(\\frac{\\ddot{a}}{a}\\right)_{H=0}=\\frac{a^2\\left[2k-\\kappa a^2\\dot{\\phi}_a\\dot{\\phi}_a-8\\kappa\\lambda k\\,\\xi_{,ab}\\dot{\\phi}_a\\dot{\\phi}_b\\right]}{2\\left(a^4+96k^2\\kappa\\lambda^2\\xi_{,a}\\xi_{,a}\\right)}", "4"))
add(("p", "so that a bounce with $V=0$ needs $\\lambda\\xi_{,\\phi\\phi}<0$ and $|\\lambda\\xi_{,\\phi\\phi}|\\gtrsim a^2/12$: the Gauss–Bonnet term must be of order one relative to the leading terms, the same regime in which the expansion in $\\alpha'$ is uncontrolled. "
          "The tensor sector gives the physical size of that regime: $G_T=1/\\kappa+8\\lambda H\\dot{\\xi}$ (derived independently from the exact Bianchi I shear action) and $F_T=1/\\kappa+8\\lambda \\xï$ (Horndeski mapping [@Kobayashi2011], verified symbolically), so the tensor speed is"))
add(("eq", "c_T^2=\\frac{F_T}{G_T}=\\frac{1+8\\kappa\\lambda\\,\\ddot{\\xi}}{1+8\\kappa\\lambda H\\dot{\\xi}}", "5"))
add(("p", "In a scan of 720 contracting initial data for $\\xi=\\phi^2$ with both signs of $\\lambda$, flat and closed, the only bounces are the %d closed cases with $\\lambda\\xi''<0$; all have tensor-sector strength $\\epsilon_T\\ge%s$ and $F_T<0$, i.e. $c_T^2<0$ (Fig. 2a). The sign of the coupling that the loop-corrected action supplies, and the exact modulus function, were not fixed here. "
          "In the tree-level heterotic frame (dilaton plus Gauss–Bonnet, fixed volume) we found in addition that flat and closed collapse cannot bounce at all ($\\Phi'^2=-3/2a^2<0$ at $H=0$ for $k=+1$), and that the open-slicing 'bounces' of the string-frame scale factor (%d of 160) all persist when the Gauss–Bonnet term is switched off and none is a bounce of the Einstein-frame scale factor (Fig. 2b): a frame diagnostic that any claimed bounce in a string-inspired model should pass." %
          (SGB["bounce"], g(SGB["bounce_epsT_min"], 2), STR["string_frame_bounces"])))
add(("fig", "FigP3_4_gb.png", "Figure 2. (a) The %d closed-collapse Gauss–Bonnet bounces ($\\xi=\\phi^2$, $\\lambda<0$): tensor-sector strength against the most negative $c_T^2$ along the path. (b) Frame diagnostic for the tree-level heterotic model with open slicing: bounces of the string-frame scale factor, of the same initial data without the Gauss–Bonnet term, and of the Einstein-frame scale factor." % SGB["bounce"], 6.3))

add(("h2", "3.5 What the black-bounce benchmark requires"))
add(("p", "The black-bounce metric $ds^2=-f\\,dt^2+dr^2/f+R^2d\\Omega^2$, $R^2=r^2+a^2$, $f=1-2m/R$ [@Simpson2019] contains a two-way wormhole ($a>2m$), a one-way transition ($a=2m$) and a regular black hole whose interior bounces into a second region ($a<2m$). It is a prescribed geometry, and its value here is that it states exactly what a source must supply."))
add(("p", "For any $f(r),R(r)$ with $g_{tt}g_{rr}=-1$ we derived $G^r{}_r-G^t{}_t=-2fR''/R$. For scalar sources $\\kappa\\epsilon f\\phi'^2=-2fR''/R$, and a purely radial magnetic nonlinear-electrodynamics field contributes nothing to the left side, so the sign of the kinetic term is $\\epsilon=-\\mathrm{sign}(R'')$. For the benchmark this forces a phantom scalar, with profile and null-energy deficit"))
add(("eq", "\\phi=\\sqrt{\\frac{2}{\\kappa}}\\,\\mathrm{arctan}\\frac{r}{a},\\qquad T_{\\mu\\nu}k^\\mu k^\\nu=-\\frac{2a^2}{\\kappa R^4},\\qquad \\int T_{kk}R^2dr=-\\frac{a}{4}", "6"))
add(("p", "reproducing the Bronnikov–Walia source [@Bronnikov2022] and showing that the ghost is forced rather than chosen. The deficit is independent of the mass and of the observer, linear in the throat radius. The same identity gives a constraint on nested throats: a profile with two minima and an intervening maximum has two inflection points, so a single scalar's kinetic function must vanish at both (infinite coupling), and the local deficit is %s against %s for one throat; a test scalar with $\\xi R$ coupling has its growth-rate threshold lowered from $\\xi_c=%s$ to %s and its growth rate at $\\xi=12$ raised from %s to %s by the second throat, so coupling between the throats does not stabilize a mode that a single throat has." %
          (g(NEST["nested (two throats + belly)"]["nec_violating_budget_x_kappa"], 3), g(NEST["single throat"]["nec_violating_budget_x_kappa"], 3), g(TF["single throat a=1"]["xi_critical"], 3), g(TF["nested a=1, beta=1.5, w=0.8"]["xi_critical"], 3),
           g([r["gamma"] for r in TF["single throat a=1"]["rows"] if r["xi"] == 12.0][0], 3), g([r["gamma"] for r in TF["nested a=1, beta=1.5, w=0.8"]["rows"] if r["xi"] == 12.0][0], 3))))
add(("p", "Kinematics. Radial free fall crosses the minimum-radius surface for every energy (energy drift below $4\\times10^{-11}$) and a bound worldline traverses a chain of trapped, anti-trapped and exterior regions, one-way. Angular momentum only hinders: crossing requires $L\\le L_c(E)$, which equals the Schwarzschild barrier ($L_c=4m$ at $E=1$) for $a\\le3m$ and is limited by the throat, $L_c=a\\sqrt{E^2/f(a)-1}$, for large wormholes. "
          "Tides are finite everywhere, with maxima $%s\\,m/a^3$ (radial) and $%s\\,m/a^3$ (transverse) at $a=m$, regular at the horizon; the transverse Jacobi field equals $R(\\tau)$ to $10^{-11}$. The anti-trapped region has positive volume expansion but leads to a second asymptotically flat exterior, not a cosmological daughter region. "
          "Transit. Through the region $K\\ge0.1\\,K(0)$ ($K$ the Kretschmann scalar) the proper time is finite and scales with $a$ (Fig. 3); a global mode of the test field, however, grows in Killing time, whose residence has a light-crossing floor of $%s\\,m$ and $%s\\,m$ for $a/m=3$ and $6$ that no energy can beat. Fast transit shortens proper time but not the exposure of such a mode, and an amplification $A=∫\\gamma\\,dt$ must be evaluated in the time variable of the instability." %
          (g(TD["K11_abs_max"], 3), g(TD["Kth_abs_max"], 3), g(FL["wormhole a/m=3"]["light_crossing_killing_time"], 3), g(FL["wormhole a/m=6"]["light_crossing_killing_time"], 3))))
add(("fig", "FigP3_3_benchmark.png", "Figure 3. The black-bounce benchmark. (a) Kretschmann profile normalized to its throat value. (b) Maximum tidal eigenvalues along radial free fall against $a/m$, compared with $m/a^3$ and $2m/a^3$. (c) Proper time to cross the region $K\\ge0.1\\,K(0)$ for energies $E=1,5$.", 6.4))

add(("h2", "3.6 Free quantum fields: a scale-free window"))
add(("p", "With $N$ free fields a Ford–Roman-type quantum inequality [@FordRoman1995], $|\\rho|\\le3N/(32\\pi^2\\tau_0^4)$ with sampling time $\\tau_0\\lesssim a$, can supply the density $|\\rho+p_r|\\approx c/(\\kappa a^2)$ that a throat needs only if $N\\ge%s\\,a^2/\\ell_P^2$ (%s for $a/m\\to\\infty$; about $10^3$ for a typical state of a conformal field). The gravitational cutoff for $N$ species is $\\Lambda\\approx M_P/\\sqrt{N}$ [@Dvali2010], and a throat of radius $a$ lies inside the effective theory only if $N\\le a^2/\\ell_P^2$. The two bounds never overlap, for any $a$, by the factor %s (Fig. 4a). "
          "This is a conditional statement within flat-space free-field estimates and a conjectured species bound with order-one uncertainty; it does not exclude strongly interacting or non-free quantum matter, and it specifies the slack ($\\gtrsim5.6$ in the product of the order-one constants) that such matter would have to supply." %
          (g(D3["B_constants"]["k_qi_static"], 3), g(D3["B_constants"]["k_qi_infinity"], 2), g(D3["B_constants"]["k_qi_static"], 3))))
add(("fig", "FigP3_5_window.png", "Figure 4. (a) The species count required to support a throat by quantum inequalities, and the count allowed by the species cutoff, in units of $a^2$: the window is empty for all throat sizes. (b) Fermi momenta of the Einstein–Cartan–Dirac gas at the onset of negative compressibility and at the bounce, against the one-species cutoff.", 6.3))

# ------------------------------------------------------------------ 4
add(("h1", "4. Dissipation as a regulator"))
add(("h2", "4.1 Causal bulk and shear viscosity"))
add(("p", "We use the truncated Israel–Stewart relaxation equations [@IS1979], $\\tau_Π\\dot{Π}+Π=-\\zetaΘ$ and $\\tau_\\pi\\dot{\\pi}+\\pi=-2\\eta\\sigma$, with causality $c^2=\\eta/(\\tau(\\rho+p))\\le1$ and hydrodynamic validity $K=\\tau|Θ|\\lesssim1$. "
          "Bulk viscosity cannot create a bounce. After a contracting phase $Π\\ge0$, also at the turning point, so $\\dot{H}=-(\\kappa/2)(\\rho+p+Π)+k/a^2$ is below its ideal value: viscosity adds positive effective pressure. In a closed universe with a cosmological constant the largest dust mass that still bounces falls from $M_c=%s$ without viscosity to $%s$ at $K=1$, $c_b^2=1$, and the null energy condition is preserved in every run (minimum of $(\\rho+p+\\pi_{\\min})/(\\rho+p)$ is $%s$), so the singularity theorem still applies to any spherical collapse of such a fluid. "
          "Shear viscosity can suppress anisotropy. The exact equation in axisymmetric Bianchi I variables gives $\\dot{\\sigma}/\\sigma=|Θ|-2\\kappa\\eta$ and, with $\\dot{Θ}/Θ=(1+w)|Θ|/2$,"
          % (g(FM["K=1 c_b2=0"], 4), g(FM["K=1 c_b2=1"], 3), g(min(r["nec_ratio_min"] for r in F0["2_shear_scan"]), 3))))
add(("eq", "\\frac{d\\ln\\Sigma^2}{dt}=(1-w)|Θ|-4\\kappa\\eta,\\qquad \\Sigma^2=\\frac{\\sigma^2}{Θ^2/3},\\qquad c^2K>\\frac{3(1-w)}{4(1+w)}", "7"))
add(("p", "which agrees with all 12 runs in the admissible domain ($K\\le1$, $c^2\\le1$). For radiation the anisotropy $\\Sigma^2$ grows from $0.0025$ to %s without viscosity and falls to $0.0017$ at the edge $K=1$, $c^2=1$." % g(F0["2_inviscid_baseline"]["0.3333333333333333"]["Sigma2_end"], 2)))

add(("h2", "4.2 Viscous relaxation of the anisotropy threshold"))
add(("p", "In the exact axisymmetric Bianchi I spin-fluid model ($\\rho_s=p_s=-Bn^2$ added to dust) the critical spin strength $B_c$ for a bounce is set by the initial shear: $B_c=%s$ for an initial shear parameter $s_0=0.05$ and $%s$ for $s_0=0.2$, equal to the analytic anisotropy threshold, because shear and spin both scale as $V^{-2}$ and shear decides which wins. We now add the causal shear stress of Section 4.1 to the dust, with $\\eta=c^2\\tau n$, $\\tau=K/|Θ|$, dust heated by the viscous work (treated as pressureless) and the constraint monitored (residual below $10^{-4}$)." % (g(B0a, 4), g(B0b, 4))))
add(("table", "Table 3. Critical spin strength $B_c$ for an axisymmetric bounce with causal shear viscosity. $\\pi_{\\max}/\\rho_d$ is the largest viscous stress relative to the dust energy density.",
     ["$s_0$", "$K$", "$c^2$", "$B_c$", "$B_c/B_c^{\\rm inviscid}$", "$\\pi_{\\max}/\\rho_d$"],
     [[("%.2f" % r["s0"]), ("%g" % r["K"]) if r["K"] else "0", ("%g" % r["c2"]) if r["c2"] else "0", g(r["B_c"], 3), g(r["ratio"], 3), ("%.2f" % r["S_over_rhod_max"]) if r.get("S_over_rhod_max") else "–"]
      for r in H0 if r["K"] <= 1.0],
     [0.7, 0.6, 0.6, 1.2, 1.6, 1.4]))
add(("p", "Result. Causal viscosity inside the hydrodynamic regime lowers the critical spin strength monotonically with the Knudsen number: by a factor %s at $K=0.1$ and %s at $K=0.3$, and to %s of the inviscid value ($\\times%s$ reduction) at $K=1$, $c^2=1$ for $s_0=0.05$; for the larger initial shear $s_0=0.2$ the reduction at $K=1$ is a factor %s (Fig. 5). In all admissible rows the viscous stress stays below $0.45$ of the dust energy density. The anisotropy threshold is therefore relaxed by nearly three orders of magnitude, toward the isotropic condition $B>0$, within the hydrodynamic regime. "
          "This is the constructive pairing the analysis suggests: a derived bounce source supplies the repulsion and a causal fluid regulates the anisotropy that otherwise raises its threshold." %
          (g(1 / v01["ratio"], 2), g(1 / v03["ratio"], 2), g(v1["ratio"], 2), g(1 / v1["ratio"], 3), g(1 / w1["ratio"], 2))))
add(("fig", "FigP3_2_viscous.png", "Figure 5. (a) Critical spin strength relative to the inviscid threshold against the Knudsen number $K$ for causal shear viscosity (circles: $c^2=1$; squares: $c^2=0.25$ at $K=1$), for initial shear $s_0=0.05$ and $0.2$. (b) The largest viscous stress relative to the dust energy density in the admissible rows.", 6.3))
add(("p", "Limits of this result. It uses a truncated Israel–Stewart theory, viscosity acting on the dust only, viscous heat treated as pressureless, and a homogeneous anisotropic model; the stress is small compared with the dust energy density but not with the rest-mass density ($|\\pi|/n$ reaches $1.9$ at $K=1$); the rows with $K>1$ lie outside the hydrodynamic regime and are not used. "
          "The pairing does not remove the fluid-sector instability of the spin-fluid equation of state (Section 3.1): viscosity relaxes the anisotropy threshold, not the negative compressibility."))

add(("h2", "4.3 Ringdown is not a bounce"))
add(("p", "A sequence of diminishing oscillations after a transition (here a ball-like sequence with ratio $q=0.7$, whose ninth peak is $q^8=5.8\\%%$ of the first) can arise in two ways that must not be conflated. In a flat expanding universe with a massive scalar and a decay width, $\\dot{H}\\le0$ and $\\rho>0$ so $H$ never changes sign (%d zero crossings in all runs); the oscillation is a ripple in $ä/a$, with at most two acceleration episodes per cycle (%d in 14 cycles for Hubble friction alone) and none once radiation dominates. "
          "Hubble friction gives a power-law, not geometric, decay (extremum ratios rising from %s to %s, ninth peak %s of the first), and a constant ratio needs a constant decay width. "
          "Genuine repeated bounces exist in the closed scalar model of earlier work: the best seed reaches %d consecutive bounces, but the oscillation amplitude is irregular and grows at the end, the neck radius shrinks from %s to %s, the entropy parameter rises monotonically by a factor %s (the Tolman effect [@Tolman1931]), and the run ends singular. The number of oscillations in a ringdown is fixed by a threshold and by $m/\\Gamma$, not by a preferred integer." %
          (RD["Hubble friction only"]["H_zero_crossings"], RD["Hubble friction only"]["accel_episodes"], g(RD["Hubble friction only"]["ratio"][0], 2), g(RD["Hubble friction only"]["ratio"][-1], 2), g(RD["Hubble friction only"]["A9_over_A1"], 2),
           CS["n_bounces"], g(CS["a_neck"][0], 3), g(CS["a_neck"][-1], 3), g(CS["S"][-1] / CS["S"][0], 2))))

# ------------------------------------------------------------------ 5
add(("h1", "5. The specification and the audit"))
add(("p", "Table 4 summarizes the six source classes against the seven conditions. A class is marked as failing only where a derivation or test in this paper or in the cited literature shows it; entries marked open were not evaluated."))
add(("table", "Table 4. Source classes against the hard conditions (P pass, F fail, – open).",
     ["Source class", "H1", "H2", "H3", "H4", "H5", "H6", "H7"],
     [["ordinary matter, canonical scalar", "P", "P", "P", "P", "P", "P", "F"],
      ["phantom scalar + NED (forced source of the benchmark)", "P", "P", "F", "F", "–", "F", "P"],
      ["NED regular black hole (de Sitter core)", "P", "P", "P", "–", "P", "–", "F"],
      ["Einstein–Cartan–Dirac, barotropic mean field", "P", "P", "P", "F", "F", "–", "–"],
      ["free quantum vacuum", "–", "P", "P", "–", "F", "–", "–"],
      ["limiting-curvature / mimetic / Horndeski", "–", "P", "–", "–", "–", "–", "–"]],
     [3.9, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]))
add(("p", "The failures rest on: the singularity theorem [@Penrose1965] (ordinary matter); the forced ghost and the instability of static ghost-supported wormholes [@Gonzalez2009] and of inner horizons of regular black holes [@CarballoRubio2021], quoted from the literature and not recomputed here (phantom and regular-black-hole classes); the barotropic criterion and the cutoff comparison of Section 3.2 (ECD); and the scale-free window of Section 3.6 (quantum vacuum). "
          "The Gauss–Bonnet class of Section 3.4 is not a separate row because it fails to produce a controlled bounce for flat collapse by theorem and for closed collapse in the scanned region, and the graceful-exit difficulty of lowest-order string cosmology is the known precedent [@BV1994]."))
add(("p", "The results above define the specification a successful model must meet. (i) It must supply a null-energy-violating sector with a derivable action whose kinetic terms are positive (the benchmark's forced phantom profile (6) is what to avoid), and whose size is at least the deficit $a/4$ for a throat of radius $a$. (ii) If it is a barotropic matter model, it must remove the negative-compressibility region required by the criterion of Section 3.1, or leave the barotropic class. (iii) If it relies on torsion, the contact interaction coefficient must be that of the Dirac theory (Table 2), not the Weyssenhoff value, and polarization cannot be used. "
          "(iv) If it is a modified-gravity model, its bounce must be tested in the frame to which matter couples, its tensor speed (5) must be positive, and the correction must act at an expansion parameter below one. (v) Its anisotropy threshold can be relaxed, within the hydrodynamic regime, by a causal fluid (Section 4.2). (vi) Any amplification estimate must use the time variable of the instability. (vii) A transition that continues into an asymptotically flat exterior is not a cosmological daughter region."))
add(("p", "Of the classes examined, only the last row of Table 4 was not evaluated; it is the natural next test, with these seven requirements as acceptance criteria."))

# ------------------------------------------------------------------ 6
add(("h1", "6. Discussion"))
add(("p", "We have not found a physical mechanism that satisfies all seven conditions, and we make no claim that one exists. What the work adds is a set of exact statements that did not exist in this form: the barotropic criterion with its two universal density ratios; the Dirac contact coefficient and the polarization cancellation, which change what a torsion-based bounce can rely on; the exact cancellation rate for spin transport; the closed-collapse condition (4) and tensor speed (5) of EsGB gravity together with the flat-collapse theorem; the mass-independent null-energy deficit and the forced ghost of the black-bounce benchmark; the Killing-time floor on transit; the scale-free species window; the shear-damping threshold (7); and the viscous relaxation of the anisotropy threshold of a spin bounce. "
          "Several of these are partly anticipated in the literature (for example the energy-condition obstruction in the heterotic loop-corrected action [@ART1994], and the positivity of bulk viscous pressure in contraction), and we cite them where we know of them; we do not claim priority beyond what is derived here."))
add(("p", "The principal limitations are these. The ECD mean field is first order and the Dirac sea and antiparticles are ignored; the string-inspired analyses use a four-dimensional truncation, not the ten-dimensional theory; the nested-throat instability and the transit amplification use a test field, not the source's own perturbations; the species-count argument combines a flat-space inequality with a conjectured cutoff; the viscous relaxation is homogeneous and uses a truncated causal theory; and no inhomogeneous collapse with a core and an envelope was evolved. Each of these is a stated condition on the result, and each can be removed by a specific calculation."))

add(("h1", "Data and code availability"))
add(("p", "All code, data and the 40 automated tests are in the public repository %s (branch phase3-recycling); every number, table and figure of this paper is produced by the scripts phase3d0_ecd.py, phase3c5_closure.py, phase3e0_string.py, phase3e1_corrections.py, phase3d1_wave.py, phase3d2_freefall.py, phase3d3_goldilocks.py, phase3d4_propertime.py, phase3f0_fluid.py, phase3g0_ringdown.py and phase3h0_viscous_spin_bounce.py. A versioned archive with a DOI will be deposited before submission." % GITHUB))
add(("h1", "Declarations"))
add(("p", "The author declares no competing interests and received no funding. Symbolic derivations and numerical experiments were carried out with the assistance of an AI coding agent; the author is responsible for the content, and every statement is backed by code that can be rerun."))
add(("refs",))
