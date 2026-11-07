# -*- coding: utf-8 -*-
"""Phase 3D.3: source-first Goldilocks feasibility.  G = c = hbar = 1 (Planck units), kappa = 8 pi.

Step 1  Hard physical constraints H1-H7 (below).  Step 2  source census: what field content does the Simpson-Visser benchmark require, and which candidate sources can provide it?
Part A  (symbolic) the NEC-violating part of ANY source of the benchmark is a phantom (ghost) scalar with profile phi = sqrt(2/kappa) arctan(r/a): derived from G^r_r - G^t_t.
Part B  Goldilocks map for the only source class with a computable two-sided competition (free quantum fields): vacuum-energy SUPPLY (quantum inequality) versus EFT VALIDITY (species cutoff); scale-free closure margin.
Part C  census table of candidate sources against the hard constraints, with the evidence used (computed here / earlier phases / literature, flagged).
Part D  decision: no admissible source model established -> stop before a parameter scan.
"""
import json
import os

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
kappa = 8 * np.pi

# ------------------------------------------------------------------------------------------------ Part A: field content forced by the benchmark
t, r, th, ph = sp.symbols("t r theta phi", real=True)
a, m = sp.symbols("a m", positive=True)
Rs = sp.sqrt(r ** 2 + a ** 2); fs = 1 - 2 * m / Rs
X = [t, r, th, ph]
g = sp.diag(-fs, 1 / fs, Rs ** 2, Rs ** 2 * sp.sin(th) ** 2); gi = g.inv()
Gam = [[[sum(gi[i, l] * (sp.diff(g[l, j], X[k]) + sp.diff(g[l, k], X[j]) - sp.diff(g[j, k], X[l])) for l in range(4)) / 2 for k in range(4)] for j in range(4)] for i in range(4)]
def Rup(i, j, k, l):
    e = sp.diff(Gam[i][j][l], X[k]) - sp.diff(Gam[i][j][k], X[l])
    return e + sum(Gam[i][k][s_] * Gam[s_][j][l] - Gam[i][l][s_] * Gam[s_][j][k] for s_ in range(4))
Ric = sp.Matrix(4, 4, lambda j, l: sum(Rup(i, j, i, l) for i in range(4)))
# G^r_r - G^t_t = R^r_r - R^t_t (the Ricci-scalar terms cancel)
diff_rt = sp.simplify(sum(gi[1, 1] * Ric[1, 1] for _ in [0]) - sum(gi[0, 0] * Ric[0, 0] for _ in [0]))
print("A: (G^r_r - G^t_t) =", sp.simplify(diff_rt), " ; expected -2 f a^2 / R^4 :", sp.simplify(diff_rt + 2 * fs * a ** 2 / Rs ** 4))
out["A_Grr_minus_Gtt_residual"] = str(sp.simplify(diff_rt + 2 * fs * a ** 2 / Rs ** 4))
# a minimally coupled scalar  L = -(eps/2) (d phi)^2 - V :  T^r_r - T^t_t = eps f phi'^2 ; a purely radial magnetic NED has T^t_t = T^r_r.  Hence:  kappa * eps * f * phi'^2 = G^r_r - G^t_t = -2 f a^2/R^4
phip = sp.symbols("phip")
eps_sol = sp.solve(sp.Eq(kappa * (-1) * fs * phip ** 2, diff_rt), phip)
print("A: with eps = -1 (phantom) the profile is phi' =", [sp.simplify(x) for x in eps_sol], "  ; eps = +1 (canonical) would need phi'^2 < 0")
out["A_phi_prime"] = [str(sp.simplify(x)) for x in eps_sol]
phi_expr = sp.sqrt(2 / sp.Float(kappa)) * sp.atan(r / a)
chk = sp.simplify(sp.diff(phi_expr, r) ** 2 * kappa - 2 * a ** 2 / Rs ** 4)
print("A: phi = sqrt(2/kappa) arctan(r/a): kappa phi'^2 - 2 a^2/R^4 =", chk)
out["A_phi_profile_residual"] = str(chk)
out["A_ghost_statement"] = ("Any source of the benchmark must contain a radial-NEC-violating sector; a canonical scalar (eps=+1) cannot give T^r_r - T^t_t < 0, NED with only a radial magnetic field has T^t_t = T^r_r. "
                            "The simplest realization is a PHANTOM scalar (eps=-1) with phi = sqrt(2/kappa) arctan(r/a) (Bronnikov-Walia 2022, verified here by the profile above): its Hamiltonian is not bounded below.")

# ------------------------------------------------------------------------------------------------ Part B: Goldilocks map for free quantum fields
W1 = json.load(open(os.path.join(HERE, "phase3d1_wave.json")))
k_qi_static = W1["C_QI_species_coefficient"]                    # N >= k a^2  (a/m = 6 value, uses rho + p_r of the static observer)
k_qi_infty = (32 * np.pi ** 2 / 3) * 2.0 / kappa                  # a/m -> infinity limit (f -> 1): kappa a^2 |rho + p_r| = 2
k_typical = k_qi_static * (32 * 1920 / 3.0 / 32 * 1.0)            # typical-state (Page/DeWitt-Schwinger magnitude N/(1920 pi^2 a^4)) instead of the extremal QI bound: (3/(32 pi^2)) / (1/(1920 pi^2)) = 180
k_typical = k_qi_static * 180.0
print("B: supply: N >= %.2f a^2 (QI, a/m=6), %.2f a^2 (a/m -> inf), %.0f a^2 (typical state, 1/(1920 pi^2 a^4) per species)" % (k_qi_static, k_qi_infty, k_typical))
# validity: species cutoff  Lambda = M_P / sqrt(N)  -> require curvature scale 1/a <= Lambda  <=>  N <= a^2 (G = 1 so M_P = 1);  reduced Planck mass version: Lambda = M_P,red/sqrt(N), M_P,red = 1/sqrt(8 pi)
conv = {"M_P": 1.0, "M_P_reduced": 1.0 / np.sqrt(8 * np.pi)}
la = np.linspace(0.0, 40.0, 401); lN = np.linspace(0.0, 90.0, 901)
A_, N_ = np.meshgrid(10 ** la, 10 ** lN)
results = {}
for supply_name, kk in (("QI a/m=6", k_qi_static), ("QI a/m->inf", k_qi_infty), ("typical state", k_typical)):
    for cn, Mp in conv.items():
        supply = N_ >= kk * A_ ** 2
        valid = N_ <= (Mp * A_) ** 2                                   # a >= sqrt(N)/Mp
        adm = supply & valid
        margin = kk / Mp ** 2
        results["%s | cutoff %s" % (supply_name, cn)] = dict(admissible_grid_points=int(adm.sum()), grid_points=int(adm.size), closure_margin_k_over_Mp2=float(margin))
        print("B: supply %-14s validity %-12s: admissible grid points = %d / %d ; the window would need the combined O(1) constants to be %.3g times more favourable (scale-free, independent of a)" % (supply_name, cn, adm.sum(), adm.size, margin))
out["B_goldilocks_map"] = results
out["B_constants"] = dict(k_qi_static=float(k_qi_static), k_qi_infinity=float(k_qi_infty), k_typical=float(k_typical))
# robustness against the O(1) constants: log-uniform slack factors in [1/3, 3] on the supply coefficient and on the cutoff relation (prior is arbitrary: this only measures how much slack is needed)
rng = np.random.default_rng(1)
u = np.exp(rng.uniform(-np.log(3), np.log(3), size=(2, 400000)))
frac = float(np.mean(k_qi_static * u[0] <= 1.0 * u[1]))
frac_inf = float(np.mean(k_qi_infty * u[0] <= 1.0 * u[1]))
print("B: Monte Carlo over O(1) constant uncertainty (log-uniform [1/3,3] each): fraction with a non-empty window: %.4f (a/m=6 coefficient), %.4f (a/m->inf)" % (frac, frac_inf))
out["B_slack_fraction_nonempty"] = dict(qi_a_over_m_6=frac, qi_a_over_m_inf=frac_inf)
# ECD scale check from Phase 3D.0
E0 = json.load(open(os.path.join(HERE, "phase3d0_ecd.json")))["C_baseline"]
ecd = {}
for k_, v in E0.items():
    kF_b = v["kF_over_Planck_at_bounce"]; kF_sp = kF_b * v["nsp_over_nb"] ** (1 / 3)
    ecd[k_] = dict(kF_bounce=kF_b, kF_spinodal=float(kF_sp), over_cutoff_MP=float(kF_sp / 1.0), over_cutoff_MPred=float(kF_sp / conv["M_P_reduced"]))
    print("B: ECD %s: k_F at spinodal = %.2f M_P (%.1f M_P,red), at bounce %.2f M_P: already at/above the one-species cutoff" % (k_, kF_sp, kF_sp / conv["M_P_reduced"], kF_b))
out["B_ecd_vs_cutoff"] = ecd

# ------------------------------------------------------------------------------------------------ Part C: source census (hard constraints H1-H7)
H = {
    "H1": "satisfies the field equations (by construction or by evolution)",
    "H2": "the NEC-violating sector, if needed, has a derivable action",
    "H3": "no ghost: Hamiltonian bounded below",
    "H4": "perturbatively stable (no exponentially growing mode, no loss of hyperbolicity)",
    "H5": "inside the effective-field-theory validity regime",
    "H6": "arises dynamically from collapse of ordinary low-density matter",
    "H7": "future-directed continuation into an expanding region (not merely geodesic crossing)",
}
P, F, O = "pass", "FAIL", "open"
census = {
    "S1 ordinary matter / canonical scalar (GR)": dict(H1=P, H2=P, H3=P, H4=P, H5=P, H6=P, H7=F,
        why="NEC holds => Penrose theorem: trapped surface => singularity (3A: 65/72 trapped, 0/42 expanding regions, LTB 96 shells all trapped); no regular transition", evidence="computed 3A + theorem"),
    "S2 phantom scalar + magnetic NED (Bronnikov-Walia source of the benchmark)": dict(H1=P, H2=P, H3=F, H4=F, H5=O, H6=F, H7=P,
        why="the profile phi = sqrt(2/kappa) arctan(r/a) with eps = -1 is FORCED (Part A); ghost; static Ellis-Bronnikov-type wormholes carry an unstable radial mode; not produced from ordinary matter", evidence="Part A computed; instability from literature (search-verified secondary sources), not recomputed"),
    "S3 NED regular black hole (de Sitter or Minkowski core)": dict(H1=P, H2=P, H3=P, H4=O, H5=P, H6=O, H7=F,
        why="satisfies NEC-compatible regular cores but no bounce into an expanding second region; inner-horizon (mass-inflation) instability generically reported; 'recycling' channel absent", evidence="literature (Carballo-Rubio et al.), not recomputed"),
    "S4 Einstein-Cartan-Dirac, barotropic mean field": dict(H1=P, H2=P, H3=P, H4=F, H5=F, H6=O, H7=O,
        why="3D.0: bounce at k_F = 2-3 M_P (above cutoff) and negative compressibility from 0.105 n_b (theorem); spinodal itself at k_F = 1.45 M_P", evidence="computed 3D.0"),
    "S5 free-field quantum vacuum (semiclassical)": dict(H1=O, H2=P, H3=P, H4=O, H5=F, H6=O, H7=O,
        why="supply N >= 5.6 a^2 versus validity N <= a^2: empty for every a (Part B), margin 5.6; no specified quantum state solved", evidence="Part B computed (order of magnitude, heuristic constants)"),
    "S6 modified gravity (limiting-curvature / mimetic / Horndeski)": dict(H1=O, H2=P, H3=O, H4=O, H5=O, H6=O, H7=O,
        why="derivable actions exist; for homogeneous dust the limiting-curvature bounce reproduces the 3B effective fluid (rho = n - B n^2) already shown to have the c_s^2 / mu = 0 pathology; stability and inhomogeneous collapse not evaluated here", evidence="3B/3C (homogeneous), otherwise not evaluated"),
}
out["C_hard_constraints"] = H
out["C_census"] = census
for k_, v in census.items():
    fails = [h for h in H if v[h] == F]
    opens = [h for h in H if v[h] == O]
    print("C: %-70s FAIL: %s | open: %s" % (k_, ",".join(fails) or "-", ",".join(opens) or "-"))
passes_all = [k_ for k_, v in census.items() if all(v[h] == P for h in H)]
out["C_sources_passing_all_hard_constraints"] = passes_all
out["D_decision"] = ("NO-GO. No source in the census satisfies H1-H7; the only quantitatively closed competition (free quantum vacuum, Part B) has an empty window for every throat size, with a scale-free closure margin of 5.6 or more; "
                     "the only derivable source of the benchmark is a ghost. A parameter scan would only map the benchmark metric, so none was run.")
print("D:", out["D_decision"])
json.dump(out, open(os.path.join(HERE, "phase3d3_goldilocks.json"), "w"), indent=1, default=str)
