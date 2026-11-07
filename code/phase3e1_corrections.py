# -*- coding: utf-8 -*-
"""Phase 3E.1: one-loop-type (field-dependent) Gauss-Bonnet couplings, nested throats, and effective propagation speeds.  Einstein-frame EsGB:
    S = \int d^4x sqrt(-g) [ R/(2 kappa) - (1/2) sum_i (d phi_i)^2 - V(phi) + lambda xi(phi) G ]            (lambda xi = coefficient of the Gauss-Bonnet scalar)
This is the structure of the loop-corrected heterotic action of Antoniadis-Rizos-Tamvakis (modulus-dependent Gauss-Bonnet coupling) with arbitrary xi(phi_i), V and any number of canonical scalars; the exact form of xi is NOT needed for Part A.
Frame statement: for canonical scalars minimally coupled to this metric the Einstein frame is the physical one; matter coupling to a different (string-frame) metric is not considered.

Part A  (symbolic) homogeneous FLRW(k): constraint, acceleration at H = 0 for arbitrary xi(phi,chi), V(phi,chi): flat collapse cannot bounce; closed collapse bounce needs the Gauss-Bonnet term at O(1) relative strength;
        numerical closed-collapse scan for xi = phi^2 (even power, Kanti-Rizos-Tamvakis type): bounce, relative GB strength eps_T at the bounce, tensor ghost / gradient conditions along the path.
Part B  nested throats: identity G^r_r - G^t_t = -2 f R''/R for generic f(r), R(r); two-throat profile; required sign changes of the scalar kinetic function; test-field instability of single vs nested throat.
Part C  effective propagation speeds: tensor modes G_T = 1/kappa + 8 lambda H xi', F_T = 1/kappa + 8 lambda xi_ddot (G_T derived independently from Bianchi I; F_T from the Horndeski mapping, checked symbolically); c_T^2 = F_T/G_T.
"""
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.linalg import eigh_tridiagonal
from sympy.calculus.euler import euler_equations

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}

# ------------------------------------------------------------------------------------------------ Part A: symbolic reduction, two scalars, arbitrary xi and V
t = sp.Symbol("t")
kap, lam, k = sp.symbols("kappa lambda k", real=True)
N = sp.Function("N")(t); a = sp.Function("a")(t); ph = sp.Function("phi")(t); ch = sp.Function("chi")(t)
xi = sp.Function("xi"); V = sp.Function("V")
ad = sp.diff(a, t); Nd = sp.diff(N, t)
Rsc = 6 * (sp.diff(a, t, 2) / (a * N ** 2) - ad * Nd / (a * N ** 3) + ad ** 2 / (a ** 2 * N ** 2) + k / a ** 2)
sqrtgG = 24 * (ad ** 2 / N ** 2 + k) * sp.diff(ad / N, t)
Lag = a ** 3 * N * (Rsc / (2 * kap) + (sp.diff(ph, t) ** 2 + sp.diff(ch, t) ** 2) / (2 * N ** 2) - V(ph, ch)) + lam * xi(ph, ch) * sqrtgG
eqs = euler_equations(Lag, [N, a, ph, ch], t)
def at1(e):
    e = e.lhs
    for n_ in (3, 2, 1):
        e = e.subs(sp.diff(N, t, n_), 0)
    return e.subs(N, 1).doit()
EN, Ea, Ep, Ec = [sp.simplify(at1(e)) for e in eqs]
A, H, P, C_, App, Pdd, Cdd = sp.symbols("A H Pd Cd addot Pdd Cdd")
Ph_, Ch_ = sp.symbols("Phi Chi")
def sym(e):
    for fn, (d1, d2) in ((ph, (P, Pdd)), (ch, (C_, Cdd))):
        e = e.subs(sp.diff(fn, t, 2), d2).subs(sp.diff(fn, t), d1)
    e = e.subs(sp.diff(a, t, 2), App).subs(sp.diff(a, t), H * A).subs(a, A)
    return e.subs(ph, Ph_).subs(ch, Ch_).doit()
ENs, Eas, Eps, Ecs = [sym(e) for e in (EN, Ea, Ep, Ec)]
cons0 = sp.simplify(ENs.subs(H, 0))
print("A: constraint at H = 0 (arbitrary xi(phi,chi), V(phi,chi)):", cons0)
out["A_constraint_H0"] = str(cons0)
# solve the a, phi, chi equations for (addot, Pdd, Cdd) at H = 0 and eliminate V via the constraint
sol = sp.solve([Eas.subs(H, 0), Eps.subs(H, 0), Ecs.subs(H, 0)], [App, Pdd, Cdd], dict=True)[0]
add0 = sp.simplify(sol[App] / A)
print("A: addot/a at H = 0 =", add0)
out["A_addot_over_a_H0"] = str(add0)
# the potential enters Raychaudhuri only through derivatives that multiply the field accelerations; use the constraint to eliminate V
Vsym = sp.Symbol("V0")
cons_V = sp.solve(cons0.subs(V(Ph_, Ch_), Vsym), Vsym)[0]
add0_k0 = sp.simplify(add0.subs(k, 0))
print("A: k = 0: addot/a at H=0 =", add0_k0)
out["A_addot_over_a_H0_flat"] = str(add0_k0)
# explicit check: for k = 0 and the constraint (phi'^2 + chi'^2)/2 + V = 0: Hdot = addot/a (H=0) = -(kappa/2)(phi'^2 + chi'^2) ; compare
target = -kap / 2 * (P ** 2 + C_ ** 2)
diff0 = sp.simplify(add0_k0.subs(V(Ph_, Ch_), -(P ** 2 + C_ ** 2) / 2) - target)
print("A: k = 0: addot/a - [-(kappa/2)(Pd^2 + Cd^2)] on the constraint surface =", diff0)
out["A_flat_residual"] = str(diff0)
# closed case: Hdot at H = 0 on the constraint surface
add0_k = sp.simplify(add0.subs(V(Ph_, Ch_), sp.solve(cons0, V(Ph_, Ch_))[0]))
print("A: general k: addot/a at H=0 on the constraint surface =", sp.simplify(add0_k))
out["A_closed_addot_over_a_H0"] = str(sp.simplify(add0_k))

# ------- single-field numerical scan: xi = phi^2 (even power), V = 0, kappa = 1 ; closed (k = +1) and flat contraction
phi_s = sp.Symbol("phi")
single = {}
xi_expr = phi_s ** 2
L1 = Lag.subs({ch: 0, sp.diff(ch, t): 0}).doit()
def build(xi_f, V_f):
    Ls = a ** 3 * N * (Rsc / (2 * kap) + sp.diff(ph, t) ** 2 / (2 * N ** 2) - V_f(ph)) + lam * xi_f(ph) * sqrtgG
    es = euler_equations(Ls, [N, a, ph], t)
    es = [sp.simplify(at1(e)) for e in es]
    def sy(e):
        e = e.subs(sp.diff(ph, t, 2), Pdd).subs(sp.diff(ph, t), P).subs(sp.diff(a, t, 2), App).subs(sp.diff(a, t), H * A).subs(a, A).subs(ph, Ph_).doit()
        return e
    return [sy(e) for e in es]
E1N, E1a, E1p = build(lambda x: x ** 2, lambda x: 0 * x)
sol1 = sp.solve([E1a, E1p], [App, Pdd], dict=True)[0]
add_f = sp.lambdify((A, H, P, Ph_, k, lam, kap), sol1[App], "numpy"); pdd_f = sp.lambdify((A, H, P, Ph_, k, lam, kap), sol1[Pdd], "numpy")
cons_f = sp.lambdify((A, H, P, Ph_, k, lam, kap), E1N / (2 * A ** 2) if False else E1N, "numpy")
cons_poly = sp.Poly(sp.expand(E1N), P)
cons_cf = [sp.lambdify((A, H, Ph_, k, lam, kap), cf, "numpy") for cf in cons_poly.all_coeffs()]
den1 = sp.simplify(sp.denom(sp.together(sol1[App])))
den_f = sp.lambdify((A, H, P, Ph_, k, lam, kap), den1, "numpy")
print("A: single-field xi = phi^2: constraint is a polynomial of degree %d in phi'; acceleration denominator:" % cons_poly.degree(), den1)
GT = lambda Hn, xid, kp, lm: 1 / kp + 8 * lm * Hn * xid
FT = lambda xidd, kp, lm: 1 / kp + 8 * lm * xidd
def run1(a0, H0, Phi0, kk, lm, root, kp=1.0, tmax=300.0):
    cfs = [f_(a0, H0, Phi0, kk, lm, kp) for f_ in cons_cf]
    rts = [x.real for x in np.roots(cfs) if abs(x.imag) < 1e-9]
    if not rts: return None
    pd0 = rts[root % len(rts)]
    def rhs(tau, y):
        A_, Ad_, P_, F_ = y
        Hn = Ad_ / A_
        return [Ad_, add_f(A_, Hn, P_, F_, kk, lm, kp), pdd_f(A_, Hn, P_, F_, kk, lm, kp), P_]
    ev_b = lambda tau, y: y[1]; ev_b.direction = 1; ev_b.terminal = True
    ev_s = lambda tau, y: y[0] - 1e-3; ev_s.terminal = True; ev_s.direction = -1
    ev_g = lambda tau, y: abs(den_f(y[0], y[1] / y[0], y[2], y[3], kk, lm, kp)) - 1e-9; ev_g.terminal = True; ev_g.direction = -1
    sl = solve_ivp(rhs, (0, tmax), [a0, H0 * a0, pd0, Phi0], method="DOP853", rtol=1e-11, atol=1e-13, events=[ev_b, ev_s, ev_g], max_step=0.1)
    outcome = "bounce" if len(sl.t_events[0]) else ("singular" if len(sl.t_events[1]) else ("degenerate" if len(sl.t_events[2]) else "none"))
    # tensor coefficients along the path: G_T = 1/kap + 8 lam H xi', F_T = 1/kap + 8 lam xi'' with xi' = d xi/dt = 2 phi phidot, xi'' = 2 phidot^2 + 2 phi phiddot
    yy = sl.y; Hn = yy[1] / yy[0]
    xid = 2 * yy[3] * yy[2]
    phidd = np.array([pdd_f(yy[0][i], Hn[i], yy[2][i], yy[3][i], kk, lm, kp) for i in range(yy.shape[1])])
    xidd = 2 * yy[2] ** 2 + 2 * yy[3] * phidd
    G_ = GT(Hn, xid, kp, lm); F_ = FT(xidd, kp, lm)
    dev = np.maximum(np.abs(8 * kp * lm * Hn * xid), np.abs(8 * kp * lm * xidd))                     # relative GB strength in the tensor sector (control parameter)
    return dict(outcome=outcome, a_min=float(yy[0].min()), eps_T_max=float(dev.max()), eps_T_end=float(dev[-1]), GT_min=float(G_.min()), FT_min=float(F_.min()), cT2_min=float((F_ / G_).min()), cT2_max=float((F_ / G_).max()),
                t_end=float(sl.t[-1]), phi_end=float(yy[3][-1]), H_end=float(Hn[-1]))
rng = np.random.default_rng(11)
scan = []
for case in range(720):
    kk = (1, 0)[case % 2]
    a0 = float(np.exp(rng.uniform(np.log(2.0), np.log(40.0)))); H0 = -float(np.exp(rng.uniform(np.log(1e-3), np.log(0.3))))
    Phi0 = float(rng.uniform(-3, 3)); lm = float(np.exp(rng.uniform(np.log(1e-3), np.log(20.0)))) * (1 if (case // 6) % 2 == 0 else -1); root = (case // 2) % 3
    try:
        r_ = run1(a0, H0, Phi0, kk, lm, root)
    except Exception:
        r_ = None
    if r_ is not None:
        r_.update(k=kk, a0=a0, H0=H0, Phi0=Phi0, lam=lm, root=root); scan.append(r_)
out["A_scan"] = scan
summ = {}
for kk, sg in ((0, 1), (0, -1), (1, 1), (1, -1)):
    sel = [r_ for r_ in scan if r_["k"] == kk and np.sign(r_["lam"]) == sg]
    bn = [r_ for r_ in sel if r_["outcome"] == "bounce"]
    summ["k=%d lam%s" % (kk, "+" if sg > 0 else "-")] = dict(n=len(sel), bounce=len(bn), singular=sum(r_["outcome"] == "singular" for r_ in sel), degenerate=sum(r_["outcome"] == "degenerate" for r_ in sel),
                         bounce_epsT_min=min((r_["eps_T_max"] for r_ in bn), default=None), bounce_epsT_max=max((r_["eps_T_max"] for r_ in bn), default=None),
                         bounce_with_epsT_below_0p3=sum(1 for r_ in bn if r_["eps_T_max"] < 0.3), bounce_with_ghost_or_gradient_instability=sum(1 for r_ in bn if r_["GT_min"] <= 0 or r_["FT_min"] <= 0))
    sm = summ["k=%d lam%s" % (kk, "+" if sg > 0 else "-")]
    print("A: scan k = %d, sign(lambda xi'') %s: %d runs: bounces %d, singular %d, degenerate %d ; eps_T at bounce %s ; bounces with eps_T < 0.3: %d ; bounces with G_T<=0 or F_T<=0 along the path: %d" %
          (kk, "+" if sg > 0 else "-", len(sel), sm["bounce"], sm["singular"], sm["degenerate"], ("%.3g..%.3g" % (sm["bounce_epsT_min"], sm["bounce_epsT_max"])) if bn else "-", sm["bounce_with_epsT_below_0p3"], sm["bounce_with_ghost_or_gradient_instability"]))
out["A_scan_summary"] = summ

# ------------------------------------------------------------------------------------------------ Part B: nested throats
r = sp.Symbol("r", real=True)
fF = sp.Function("f")(r); RF = sp.Function("R")(r); th_, vp_ = sp.symbols("theta varphi")
Xc = [sp.Symbol("t"), r, th_, vp_]
g = sp.diag(-fF, 1 / fF, RF ** 2, RF ** 2 * sp.sin(th_) ** 2); gi = g.inv()
Gam = [[[sum(gi[i, l] * (sp.diff(g[l, j], Xc[kk_]) + sp.diff(g[l, kk_], Xc[j]) - sp.diff(g[j, kk_], Xc[l])) for l in range(4)) / 2 for kk_ in range(4)] for j in range(4)] for i in range(4)]
def Rup(i, j, kk_, l):
    e = sp.diff(Gam[i][j][l], Xc[kk_]) - sp.diff(Gam[i][j][kk_], Xc[l])
    return e + sum(Gam[i][kk_][s_] * Gam[s_][j][l] - Gam[i][l][s_] * Gam[s_][j][kk_] for s_ in range(4))
Ric = sp.Matrix(4, 4, lambda j, l: sum(Rup(i, j, i, l) for i in range(4)))
diff_rt = sp.simplify(gi[1, 1] * Ric[1, 1] - gi[0, 0] * Ric[0, 0])
print("B: G^r_r - G^t_t for generic f(r), R(r) =", diff_rt)
assert sp.simplify(diff_rt + 2 * fF * sp.diff(RF, r, 2) / RF) == 0
Rsc_gen = sp.simplify(sum(gi[i, i] * Ric[i, i] for i in range(4)))
print("B: Ricci scalar for generic f, R =", Rsc_gen)
out["B_Grr_minus_Gtt"] = str(diff_rt)
out["B_ricci_scalar_general"] = str(Rsc_gen)
# kinetic function of a single scalar source: T^r_r - T^t_t = eps f phi'^2 -> kappa eps f phi'^2 = -2 f R''/R  => sign(eps) = -sign(R'')
def profiles(a_, beta, w):
    Rr = lambda x: np.sqrt(x ** 2 + a_ ** 2) + beta * np.exp(-(x / w) ** 2)
    return Rr
def derivs(Rr, x, h=1e-4):
    R0 = Rr(x); R1 = (Rr(x + h) - Rr(x - h)) / (2 * h); R2 = (Rr(x + h) - 2 * R0 + Rr(x - h)) / h ** 2
    return R0, R1, R2
nest = {}
for name, (a_, beta, w) in {"single throat": (1.0, 0.0, 1.0), "nested (two throats + belly)": (1.0, 1.5, 0.8)}.items():
    Rr = profiles(a_, beta, w); xs = np.linspace(-6, 6, 120001); R0, R1, R2 = derivs(Rr, xs)
    mins = [float(xs[i]) for i in range(1, len(xs) - 1) if R1[i - 1] < 0 <= R1[i]]
    maxs = [float(xs[i]) for i in range(1, len(xs) - 1) if R1[i - 1] > 0 >= R1[i]]
    infl = [float(xs[i]) for i in range(1, len(xs) - 1) if R2[i - 1] * R2[i] < 0]
    # NEC budget ~ integral of (-T_kk) R^2 dr = (2/kappa) int R'' R dr over regions with R'' > 0 ; and the part with R'' < 0 (NEC-satisfying belly)
    budget_viol = float(np.trapezoid(np.where(R2 > 0, R2 * R0, 0), xs) * 2); budget_ok = float(np.trapezoid(np.where(R2 < 0, -R2 * R0, 0), xs) * 2)
    nest[name] = dict(local_minima_of_R=mins, local_maxima_of_R=maxs, inflection_points=infl, nec_violating_budget_x_kappa=budget_viol, nec_satisfying_budget_x_kappa=budget_ok, Rmin=float(R0[np.argmin(np.abs(xs))]))
    print("B: %-30s minima of R at %s, maxima at %s, R'' sign changes at %s ; int(-kappa T_kk R^2 dr)_(violating) = %.3f, satisfying part %.3f" % (name, [round(x, 3) for x in mins], [round(x, 3) for x in maxs], [round(x, 3) for x in infl], budget_viol, budget_ok))
out["B_nested_profile"] = nest
out["B_kinetic_sign_statement"] = ("For a static metric with g_tt g_rr = -1 sourced by canonical-or-phantom scalars, kappa eps f phi'^2 = -2 f R''/R, so the kinetic sign eps = -sign(R'') (f>0). Each throat (R''>0) needs eps = -1 (phantom), each belly (R''<0) eps = +1: "
                                    "a single scalar with kinetic function h(phi) must change sign at every inflection point (h = 0 there): the kinetic term of perturbations vanishes (infinite coupling); two throats need at least two such zero crossings.")
# test-field instability single vs nested: xi R coupled scalar, l = 0, f = 1 (m = 0) ; V = f[xi Ricci] + (1/R)(f R')' f  with f = 1
def ground_eig(Rr, xi_c, Lb=60.0, h=0.01):
    xs = np.arange(-Lb, Lb + h / 2, h)
    R0, R1, R2 = derivs(Rr, xs, 1e-4)
    Ric_num = -4 * 0 + (-2 * R2 / R0 - 2 * R1 ** 2 / R0 ** 2 + 2 / R0 ** 2)               # f = 1: Ricci = -4 R''/R + ... computed below from the general formula
    Ric_num = (-4 * R2 / R0 - 2 * R1 ** 2 / R0 ** 2 + 2 / R0 ** 2)
    V = xi_c * Ric_num + R2 / R0
    return float(eigh_tridiagonal(2 / h ** 2 + V, -np.ones(len(xs) - 1) / h ** 2, select="i", select_range=(0, 0), eigvals_only=True)[0])
# check the f = 1 Ricci scalar against the symbolic general expression
chk = sp.simplify(Rsc_gen.subs(fF, 1).doit() - (-4 * sp.diff(RF, r, 2) / RF - 2 * sp.diff(RF, r) ** 2 / RF ** 2 + 2 / RF ** 2))
print("B: Ricci scalar (f = 1) check vs closed form:", chk)
out["B_ricci_f1_check"] = str(chk)
cmp = {}
for name, (a_, beta, w) in {"single throat a=1": (1.0, 0.0, 1.0), "nested a=1, beta=1.5, w=0.8": (1.0, 1.5, 0.8), "nested a=1, beta=3, w=0.6": (1.0, 3.0, 0.6)}.items():
    Rr = profiles(a_, beta, w)
    rows = []
    for xi_c in (0.0, 0.5, 1.0, 3.0, 12.0):
        lam0 = ground_eig(Rr, xi_c)
        rows.append(dict(xi=xi_c, lowest_eigenvalue=lam0, gamma=float(np.sqrt(-lam0)) if lam0 < -1e-6 else 0.0))
    lo, hi = 0.0, 40.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if ground_eig(Rr, mid) < -1e-5: hi = mid
        else: lo = mid
    cmp[name] = dict(rows=rows, xi_critical=hi)
    print("B: %-30s xi_crit(l=0) = %.3f ; gamma(xi) = %s" % (name, hi, [(r_["xi"], round(r_["gamma"], 3)) for r_ in rows]))
out["B_test_field_instability"] = cmp

# ------------------------------------------------------------------------------------------------ Part C: propagation characteristics
# Horndeski mapping of lambda xi(phi) G : G5 = -4 xi' ln X, G4 = 4 xi'' X (2 - ln X) + 1/(2 kappa), G3 = 4 xi''' X (7 - 3 ln X), G2 = 8 xi'''' X^2 (3 - ln X) ; G_T = 2[G4 - 2X G4X - X (H phidot G5X - G5phi)], F_T = 2[G4 - X (phiddot G5X + G5phi)]
X, lnX = sp.symbols("X lnX", positive=True)
x1, x2, Hs, pd, pdd = sp.symbols("xi1 xi2 H phid phidd")
Xs = pd ** 2 / 2
G5 = -4 * x1 * sp.log(X); G5X = sp.diff(G5, X); G5p = -4 * x2 * sp.log(X)
G4 = 4 * x2 * X * (2 - sp.log(X)) + 1 / (2 * kap); G4X = sp.diff(G4, X)
GT_h = sp.simplify((2 * (G4 - 2 * X * G4X - X * (Hs * pd * G5X - G5p))).subs(X, Xs))
FT_h = sp.simplify((2 * (G4 - X * (pdd * G5X + G5p))).subs(X, Xs))
print("C: Horndeski-mapped G_T =", sp.simplify(GT_h), " ; F_T =", sp.simplify(FT_h))
xid = x1 * pd; xidd = x2 * pd ** 2 + x1 * pdd
rG = sp.simplify(GT_h - (1 / kap + 8 * Hs * xid)); rF = sp.simplify(FT_h - (1 / kap + 8 * xidd))
print("C: G_T - (1/kappa + 8 H xi_dot) =", rG, " ; F_T - (1/kappa + 8 xi_ddot) =", rF, "  (lambda absorbed in xi)")
out["C_GT_residual"] = str(rG); out["C_FT_residual"] = str(rF)
out["C_cT2"] = "c_T^2 = F_T/G_T = (1 + 8 kappa lam xi_ddot)/(1 + 8 kappa lam H xi_dot)"
# independent derivation of G_T from the Bianchi I quadratic action (Einstein frame, single scalar)
b = sp.Function("b")(t); e_ = sp.Symbol("e")
aI = [a * sp.exp(e_ * b), a * sp.exp(-e_ * b), a]
Hs_ = [sp.diff(x, t) / x for x in aI]
RI = 2 * sum(sp.diff(x, t, 2) / x for x in aI) + 2 * (Hs_[0] * Hs_[1] + Hs_[0] * Hs_[2] + Hs_[1] * Hs_[2])
LI = aI[0] * aI[1] * aI[2] * (RI / (2 * kap) + sp.diff(ph, t) ** 2 / 2) + lam * xi(ph) * 8 * sp.diff(sp.diff(aI[0], t) * sp.diff(aI[1], t) * sp.diff(aI[2], t), t)
L2 = sp.simplify(sp.diff(LI, e_, 2).subs(e_, 0) / 2)
eb = sp.expand(euler_equations(L2, [b], t)[0].lhs)
Kc = sp.simplify(eb.coeff(sp.diff(b, t, 2)))
Kc_n = sp.simplify(Kc / (-2 * a ** 3))
GTb = sp.simplify(Kc_n)
expected = 1 / kap + 8 * lam * (sp.diff(a, t) / a) * sp.diff(xi(ph), t)
print("C: independent Bianchi-I shear-sector coefficient / (-2 a^3) =", GTb, "  vs expected 1/kappa + 8 lam H d xi/dt ; difference:", sp.simplify((GTb - expected).doit()))
out["C_GT_bianchi"] = str(GTb); out["C_GT_bianchi_residual"] = str(sp.simplify((GTb - expected).doit()))
json.dump(out, open(os.path.join(HERE, "phase3e1_corrections.json"), "w"), indent=1, default=str)
