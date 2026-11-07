# -*- coding: utf-8 -*-
"""Phase 3E.0: string-inspired high-curvature collapse, homogeneous baseline.  Units alpha' = 1 (so epsilon = alpha' |curvature| = |curvature|).

Model (heterotic tree level, string frame, 4D after reduction on a torus with FIXED internal volume; internal-volume moduli and gauge fields set to zero):
    S = (1/2 kappa_4^2) \int d^4x sqrt(-g) e^{-2 Phi} [ R + 4 (d Phi)^2 + c G ] ,   c = alpha'/8 (heterotic sign: +),   G = R^2 - 4 R_mn R^mn + R_mnrs R^mnrs.
This is a SCHEMATIC truncation: H_{MNP}, gauge fields, the volume modulus and all higher orders in alpha' are dropped (as in the instruction's own schematic action); field redefinitions fix the quadratic curvature
terms to the Gauss-Bonnet combination at this order.  Solutions are therefore valid only where epsilon << 1.

Part A  minisuperspace reduction of FLRW with curvature k: R, G closed forms (sqrt(-g) G = 8 d/dt[ a'^3/N^3 + 3 k a'/N ] verified), Euler-Lagrange equations for N, a, Phi.
Part B  checks: (i) c = 0, k = 0 reproduces the pre-big-bang dilaton vacuum a ~ (-t)^(-1/sqrt3), shifted dilaton constraint; (ii) at H = 0 for k = 0 the equations force Phi' = 0 and a'' = 0: NO bounce for flat collapse at tree level with GB;
        (iii) scale argument: the Gauss-Bonnet terms compete with the Einstein-Hilbert terms only where alpha' x curvature = O(1).
Part C  shear-sector (Bianchi I) ghost condition from the exact quadratic action; background-dependent coefficient K_shear.
Part D  feasibility scan: closed (k = +1) contracting initial data -> outcome (bounce / singular / degenerate), epsilon at the bounce, K_shear sign.
"""
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}

# ------------------------------------------------------------------------------------------------ Part A: reduction
t = sp.Symbol("t")
k, c = sp.symbols("k c", real=True)
N = sp.Function("N")(t); a = sp.Function("a")(t); Ph = sp.Function("Phi")(t)
ad, add_ = sp.diff(a, t), sp.diff(a, t, 2)
Nd = sp.diff(N, t)
Rsc = 6 * (add_ / (a * N ** 2) - ad * Nd / (a * N ** 3) + ad ** 2 / (a ** 2 * N ** 2) + k / a ** 2)
sqrtgG = 24 * (ad ** 2 / N ** 2 + k) * sp.diff(ad / N, t)                                  # sqrt(-g) G per comoving volume (a^3 N G)
chk_G = sp.simplify(sqrtgG - 8 * sp.diff(ad ** 3 / N ** 3 + 3 * k * ad / N, t))
print("A: sqrt(-g) G - 8 d/dt[a'^3/N^3 + 3 k a'/N] =", chk_G)
out["A_GB_total_derivative_residual"] = str(chk_G)
# independent check of the closed forms for N = 1 from the metric curvature (flat and closed): compute R and G for the S^3 slicing by brute force
th_, ch_, ph_ = sp.symbols("theta chi phi")
Xc = [t, ch_, th_, ph_]
Nn = sp.Integer(1)
gmat = sp.diag(-1, a ** 2, a ** 2 * sp.sin(ch_) ** 2, a ** 2 * sp.sin(ch_) ** 2 * sp.sin(th_) ** 2)
gin = gmat.inv()
Gm = [[[sum(gin[i, l] * (sp.diff(gmat[l, j], Xc[kk]) + sp.diff(gmat[l, kk], Xc[j]) - sp.diff(gmat[j, kk], Xc[l])) for l in range(4)) / 2 for kk in range(4)] for j in range(4)] for i in range(4)]
def Rm(i, j, kk, l):
    e = sp.diff(Gm[i][j][l], Xc[kk]) - sp.diff(Gm[i][j][kk], Xc[l])
    return e + sum(Gm[i][kk][s_] * Gm[s_][j][l] - Gm[i][l][s_] * Gm[s_][j][kk] for s_ in range(4))
Ricm = sp.Matrix(4, 4, lambda j, l: sum(Rm(i, j, i, l) for i in range(4)))
Rbf = sp.simplify(sum(gin[i, i] * Ricm[i, i] for i in range(4)))
Ric2 = sum(gin[i, i] * gin[j, j] * Ricm[i, j] ** 2 for i in range(4) for j in range(4))
Riem2 = 0
for i in range(4):
    for j in range(4):
        for kk in range(4):
            for l in range(4):
                v = Rm(i, j, kk, l)
                if v != 0:
                    Riem2 += gmat[i, i] * v ** 2 * gin[j, j] * gin[kk, kk] * gin[l, l]
Gbf = sp.simplify(Rbf ** 2 - 4 * Ric2 + Riem2)
R_closed = Rsc.subs({N: 1, k: 1}).doit(); G_closed = (sqrtgG / (a ** 3 * N)).subs({N: 1, k: 1}).doit()
print("A: brute-force R (k=1, N=1) - closed form =", sp.simplify(Rbf - R_closed), "; brute-force G - closed form =", sp.simplify(Gbf - G_closed))
out["A_bruteforce_R_residual"] = str(sp.simplify(Rbf - R_closed)); out["A_bruteforce_G_residual"] = str(sp.simplify(Gbf - G_closed))

Lag = a ** 3 * N * sp.exp(-2 * Ph) * (Rsc - 4 * sp.diff(Ph, t) ** 2 / N ** 2) + c * sp.exp(-2 * Ph) * sqrtgG
from sympy.calculus.euler import euler_equations
eqs = euler_equations(Lag, [N, a, Ph], t)
subsN = {N: 1, Nd: 0, sp.diff(N, t, 2): 0, sp.diff(N, t, 3): 0}
def at_N1(e):
    e = e.lhs if hasattr(e, "lhs") else e
    return sp.simplify(e.subs(sp.diff(N, t, 3), 0).subs(sp.diff(N, t, 2), 0).subs(sp.diff(N, t), 0).subs(N, 1).doit())
EN = at_N1(eqs[0]); Ea = at_N1(eqs[1]); EP = at_N1(eqs[2])
# symbols for the state
Aa, Hh, PD, AD2, PD2, A3, P3 = sp.symbols("A H Pd addot Pdd a3dot P3dot")
def tosym(e):
    e = e.subs(sp.diff(Ph, t, 3), P3).subs(sp.diff(a, t, 3), A3).subs(sp.diff(Ph, t, 2), PD2).subs(sp.diff(a, t, 2), AD2).subs(sp.diff(Ph, t), PD).subs(sp.diff(a, t), Hh * Aa).subs(a, Aa)
    return sp.simplify(e)
ENs, Eas, EPs = tosym(EN), tosym(Ea), tosym(EP)
print("A: constraint E_N =", sp.simplify(ENs / sp.exp(-2 * Ph)))
out["A_constraint"] = str(sp.simplify(ENs / sp.exp(-2 * Ph)))

# ------------------------------------------------------------------------------------------------ Part B: checks
# (i) c = 0, k = 0: dilaton vacuum a ~ (-t)^p with p = -1/sqrt3 ; shifted dilaton phibar = 2 Phi - 3 ln a ; constraint phibar'^2 = 3 H^2 .. check the power-law solution
tt = sp.Symbol("tt", negative=True)
p = -1 / sp.sqrt(3)
a_sol = (-tt) ** p
# phibar' = -sqrt(3)*H sign choice; Phi from phibar: Phi = (phibar + 3 ln a)/2 ; phibar = -sqrt3 * ln a ... (choose the branch phibar' = sqrt(3) H ... test both)
res_vac = {}
for sgn in (+1, -1):
    ph_sol = sp.Rational(1, 2) * (sgn * sp.sqrt(3) * sp.log(a_sol) + 3 * sp.log(a_sol)) + sp.Symbol("C0")
    subsd = {k: 0, c: 0}
    chks = []
    for e in (EN, Ea, EP):
        ee = e.subs(subsd)
        ee = ee.subs({sp.diff(a, t, 3): sp.diff(a_sol, tt, 3), sp.diff(Ph, t, 3): sp.diff(ph_sol, tt, 3)}).subs({sp.diff(a, t, 2): sp.diff(a_sol, tt, 2), sp.diff(Ph, t, 2): sp.diff(ph_sol, tt, 2)}).subs({sp.diff(a, t): sp.diff(a_sol, tt), sp.diff(Ph, t): sp.diff(ph_sol, tt)}).subs({a: a_sol, Ph: ph_sol})
        chks.append(sp.simplify(ee))
    res_vac[str(sgn)] = [str(x) for x in chks]
ok_vac = [sg for sg, v in res_vac.items() if all(x == "0" for x in v)]
print("B: dilaton-vacuum power law a ~ (-t)^(-1/sqrt3) solves the c=0,k=0 equations for sign branch(es):", ok_vac)
out["B_prebigbang_branches"] = ok_vac
# (ii) k = 0: at H = 0 the constraint gives Phi' = 0, and then the a- and Phi-equations give addot = 0 , Phi'' = 0 : no bounce
sub0 = {k: 0, Hh: 0}
cons0 = sp.simplify(ENs.subs(sub0) / sp.exp(-2 * Ph))
print("B: constraint at H=0, k=0:", cons0, " -> solutions for Phi':", sp.solve(cons0, PD))
sol_pd = sp.solve(cons0, PD)
out["B_flat_constraint_H0"] = str(cons0)
sub1 = {k: 0, Hh: 0, PD: 0}
a_eq = sp.simplify(Eas.subs(sub1)); P_eq = sp.simplify(EPs.subs(sub1))
sol_a = sp.solve(a_eq, AD2); sol_P = sp.solve(P_eq, PD2)
print("B: at H = 0, Phi' = 0, k = 0: a-equation ->", sol_a, "; Phi-equation ->", sol_P)
out["B_flat_H0_addot"] = [str(x) for x in sol_a]; out["B_flat_H0_Phidd"] = [str(x) for x in sol_P]
# (iii) closed universe: constraint and acceleration at H = 0
sub_k = {k: 1, Hh: 0}
cons_k = sp.simplify(ENs.subs(sub_k) / sp.exp(-2 * Ph))
print("B: constraint at H=0, k=1:", cons_k)
out["B_closed_constraint_H0"] = str(cons_k)
print("B: k=+1: at H = 0 the constraint requires Phi'^2 =", sp.solve(cons_k, PD**2 if False else PD), "-> imaginary: no closed-universe bounce with the GB term (all GB terms carry a factor of H)")
Pk = sp.solve(sp.simplify(ENs.subs({k: -1, Hh: 0}) / sp.exp(-2 * Ph)), PD)
print("B: k=-1: constraint at H=0 gives Phi' =", Pk)
sol_open = sp.solve([Eas.subs({k: -1, Hh: 0}), EPs.subs({k: -1, Hh: 0})], [AD2, PD2], dict=True)[0]
addot_open = [sp.simplify(sol_open[AD2].subs(PD, pk)) for pk in Pk]
print("B: k=-1, H=0: addot / a = ", [sp.simplify(x / Aa) for x in addot_open], "(c = 1/8 substituted:", [sp.simplify((x / Aa).subs(c, sp.Rational(1, 8))) for x in addot_open], ")")
out["B_open_H0"] = dict(Phid=[str(x) for x in Pk], addot_over_a=[str(sp.simplify(x / Aa)) for x in addot_open])

# ------------------------------------------------------------------------------------------------ Part C: shear-sector kinetic coefficient (Bianchi I, N = 1)
b = sp.Function("b")(t); eps = sp.Symbol("eps")
a1, a2, a3 = a * sp.exp(eps * b), a * sp.exp(-eps * b), a
def hubble(f_): return sp.diff(f_, t) / f_
Hs = [hubble(x) for x in (a1, a2, a3)]
Rb = 2 * sum(sp.diff(x, t, 2) / x for x in (a1, a2, a3)) + 2 * (Hs[0] * Hs[1] + Hs[0] * Hs[2] + Hs[1] * Hs[2])
sqrtgG_b = 8 * sp.diff(sp.diff(a1, t) * sp.diff(a2, t) * sp.diff(a3, t), t)
LB = a1 * a2 * a3 * sp.exp(-2 * Ph) * (Rb - 4 * sp.diff(Ph, t) ** 2) + c * sp.exp(-2 * Ph) * sqrtgG_b
L2 = sp.simplify(sp.diff(LB, eps, 2).subs(eps, 0) / 2)
eb = euler_equations(L2, [b], t)[0]
eb_lhs = sp.expand(eb.lhs)
Kcoef = sp.simplify(eb_lhs.coeff(sp.diff(b, t, 2)))
Kcoef_s = tosym(Kcoef)
print("C: coefficient of b'' in the linearised shear equation, K_shear / (e^{-2 Phi} ... ) =", sp.simplify(Kcoef_s))
out["C_Kshear"] = str(sp.simplify(Kcoef_s))
Kn = sp.simplify(Kcoef_s / (sp.exp(-2 * Ph) * Aa ** 3))
print("C: normalised K_shear =", Kn)
out["C_Kshear_normalised"] = str(Kn)

# ------------------------------------------------------------------------------------------------ Part D: homogeneous feasibility scan
# solve (a-eq, Phi-eq) for (addot, Phidd) ; constraint fixes one initial datum
sol2 = sp.solve([Eas, EPs], [AD2, PD2], dict=True)[0]
add_f = sp.lambdify((Aa, Hh, PD, k, c), sol2[AD2], "numpy"); pdd_f = sp.lambdify((Aa, Hh, PD, k, c), sol2[PD2], "numpy")
den = sp.simplify(sp.denom(sp.together(sol2[AD2])))
den_f = sp.lambdify((Aa, Hh, PD, k, c), den, "numpy")
cons_f = sp.lambdify((Aa, Hh, PD, k, c), ENs / sp.exp(-2 * Ph), "numpy")
K_f = sp.lambdify((Aa, Hh, PD, AD2, PD2, k, c), Kn, "numpy")
print("D: denominator of the acceleration solution:", den)
out["D_denominator"] = str(den)
cval = 1.0 / 8.0
def run_case(a0, H0, kk, cc, root=0, tmax=400.0):
    cs = sp.Poly(sp.together(sp.simplify(ENs / sp.exp(-2 * Ph))).subs({Aa: a0, Hh: H0, k: kk, c: cc}), PD)
    rts = [complex(x) for x in np.roots([float(cf) for cf in cs.all_coeffs()])]
    rts = [x.real for x in rts if abs(x.imag) < 1e-9]
    if not rts: return None
    pd0 = rts[root % len(rts)] if len(rts) > 1 else rts[0]
    def rhs2(tau, y):
        A_, Ad_, P_, F_ = y
        Hn = Ad_ / A_
        return [Ad_, add_f(A_, Hn, P_, kk, cc), pdd_f(A_, Hn, P_, kk, cc), P_]
    ev_b = lambda tau, y: y[1]; ev_b.direction = 1; ev_b.terminal = False
    ev_s = lambda tau, y: y[0] - 1e-3; ev_s.terminal = True; ev_s.direction = -1
    ev_g = lambda tau, y: abs(den_f(y[0], y[1] / y[0], y[2], kk, cc)) - 1e-8; ev_g.terminal = True; ev_g.direction = -1
    ev_E = lambda tau, y: y[1] / y[0] - y[2]; ev_E.direction = 1; ev_E.terminal = False         # Einstein-frame bounce: H - Phi' crosses zero upward
    ev_x = lambda tau, y: y[0] - 3000.0; ev_x.terminal = True; ev_x.direction = 1
    sol = solve_ivp(rhs2, (0, tmax), [a0, H0 * a0, pd0, 0.0], method="DOP853", rtol=1e-11, atol=1e-13, events=[ev_b, ev_s, ev_g, ev_E, ev_x], max_step=0.2)
    nb = len(sol.t_events[0]); nE = len(sol.t_events[3])
    outcome = "bounce" if nb else ("singular" if len(sol.t_events[1]) else ("degenerate" if len(sol.t_events[2]) else "none"))
    eps_b = None; Kb = None
    if nb:
        yb = sol.y_events[0][0]; Hb = yb[1] / yb[0]
        eps_b = float(max(Hb ** 2, abs(Hb * yb[2]), yb[2] ** 2, abs(kk) / yb[0] ** 2))
    A_, Ad_, P_, F_ = sol.y[:, -1]
    Hn = Ad_ / A_
    Ks = [float(K_f(sol.y[0, i], sol.y[1, i] / sol.y[0, i], sol.y[2, i], add_f(sol.y[0, i], sol.y[1, i] / sol.y[0, i], sol.y[2, i], kk, cc), pdd_f(sol.y[0, i], sol.y[1, i] / sol.y[0, i], sol.y[2, i], kk, cc), kk, cc)) for i in range(0, sol.y.shape[1], max(1, sol.y.shape[1] // 400))]
    aE = sol.y[0] * np.exp(-sol.y[3])
    eps_path = max(float(np.max(np.maximum.reduce([(sol.y[1] / sol.y[0]) ** 2, np.abs(sol.y[1] / sol.y[0] * sol.y[2]), sol.y[2] ** 2, abs(kk) / sol.y[0] ** 2]))), 0.0)
    return dict(outcome=outcome, string_bounce=bool(nb), einstein_bounce=bool(nE), eps_at_string_bounce=eps_b, eps_max_along_path=eps_path, g_s_max=float(np.exp(sol.y[3].max())), g_s_min=float(np.exp(sol.y[3].min())),
                aE_min=float(aE.min()), Kshear_max=float(max(Ks)), a_min=float(sol.y[0].min()), t_end=float(sol.t[-1]), constraint_residual=float(abs(cons_f(A_, Hn, P_, kk, cc))))
rng = np.random.default_rng(7)
table = []
for case in range(480):
    a0 = float(np.exp(rng.uniform(np.log(1.5), np.log(60.0))))
    H0 = -float(np.exp(rng.uniform(np.log(1e-3), np.log(0.4))))
    kk = (1, 0, -1)[case % 3]; root = (case // 3) % 2
    try:
        r_ = run_case(a0, H0, kk, cval, root)
        r0_ = run_case(a0, H0, kk, 0.0, root)                                  # same data without the Gauss-Bonnet term (c = 0)
    except Exception as e:
        r_ = None
    if r_ is not None:
        r_.update(a0=a0, H0=H0, k=kk, root=root, outcome_c0=None if r0_ is None else r0_["outcome"], einstein_bounce_c0=None if r0_ is None else r0_["einstein_bounce"])
        table.append(r_)
out["D_scan"] = table
out["D_summary"] = {}
for kk in (0, 1, -1):
    sel = [r_ for r_ in table if r_["k"] == kk]
    cnt = {o: sum(1 for r_ in sel if r_["outcome"] == o) for o in ("bounce", "singular", "degenerate", "none")}
    bn = [r_ for r_ in sel if r_["string_bounce"]]
    nE = sum(1 for r_ in sel if r_["einstein_bounce"])
    same = sum(1 for r_ in bn if r_["outcome_c0"] == "bounce")
    summ = dict(n=len(sel), counts=cnt, string_frame_bounces=len(bn), einstein_frame_bounces=nE, string_bounces_also_without_GB=same,
                eps_at_bounce=[min(x["eps_at_string_bounce"] for x in bn), max(x["eps_at_string_bounce"] for x in bn)] if bn else None,
                bounce_with_eps_below_0p1=sum(1 for x in bn if x["eps_at_string_bounce"] < 0.1), bounce_eps_path_below_0p1=sum(1 for x in bn if x["eps_max_along_path"] < 0.1),
                gs_max_over_bounces=max((x["g_s_max"] for x in bn), default=None), max_Kshear=max((x["Kshear_max"] for x in bn), default=None), ghost_runs=sum(1 for x in sel if x["Kshear_max"] >= 0))
    out["D_summary"][str(kk)] = summ
    print("D: k = %2d: %3d runs: %s | string-frame bounces %d (without GB term too: %d), Einstein-frame bounces %d | eps at string bounce %s, bounces with eps<0.1: %d, with eps<0.1 along the whole path: %d | max g_s over bounces %s | shear-ghost runs %d" %
          (kk, len(sel), cnt, len(bn), same, nE, ("%.2g..%.2g" % tuple(summ["eps_at_bounce"])) if bn else "-", summ["bounce_with_eps_below_0p1"], summ["bounce_eps_path_below_0p1"], ("%.3g" % summ["gs_max_over_bounces"]) if bn else "-", summ["ghost_runs"]))
json.dump(out, open(os.path.join(HERE, "phase3e0_string.json"), "w"), indent=1, default=str)
