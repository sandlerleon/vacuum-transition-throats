# -*- coding: utf-8 -*-
"""Phase 2E Test B (linear part): inhomogeneous (spherically symmetric, radial harmonic n = 3, 4, ...) perturbations of the closed cold cycle with scalar + coupled dust, longitudinal gauge.
Equations derived symbolically in phase2e_perturb_derive.py / phase2e_perturb_reduce.py.  Evolved variables y = (f = d phi, f', D = number-density contrast, v = velocity potential) per mode,
Phi = Psi from the constraints (00) and (0 chi), no division by H.  Floquet multipliers = eigenvalues of the 4x4 monodromy over one period of the background cycle.
Correctness checks: constraint propagation  d/dt Phi_alg - Phi'_(0chi) = 0  and the unused trace equation (chi chi) along a solution.

    python phase2e_perturb_floquet.py          (writes phase2e_perturb_floquet.json)
"""
import json
import os
import pickle
import sys

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

from phase2d_cycle_validation import Cycle

HERE = os.path.dirname(os.path.abspath(__file__))
res = pickle.load(open(os.path.join(HERE, "phase2e_perturb_reduced.pkl"), "rb"))
t = sp.Symbol("t")
C, M2, mu, lam = sp.symbols("C M2 mu lam", positive=True)
k2 = sp.Symbol("k2", positive=True)
a_ = sp.Function("a")(t); p0_ = sp.Function("p0")(t); ph_ = sp.Function("ph")(t); f_ = sp.Function("f")(t); D_ = sp.Function("Dd")(t); v_ = sp.Function("v")(t)
A_, AD, ADD, P0, P0D, P0DD, PH, PHD, PHDD, F, FD, FDD, DD, DDD, V, VD = sp.symbols("A AD ADD P0 P0D P0DD PH PHD PHDD F FD FDD DD DDD V VD")
Qs, Qps = sp.symbols("Q Qp")


def sym(e):
    e = e.subs({sp.Derivative(ph_, (t, 2)): PHDD, sp.Derivative(f_, (t, 2)): FDD, sp.Derivative(a_, (t, 2)): ADD, sp.Derivative(p0_, (t, 2)): P0DD})
    e = e.subs({sp.Derivative(ph_, t): PHD, sp.Derivative(f_, t): FD, sp.Derivative(a_, t): AD, sp.Derivative(p0_, t): P0D, sp.Derivative(D_, t): DDD, sp.Derivative(v_, t): VD})
    return e.subs({ph_: PH, f_: F, D_: DD, v_: V, a_: A_, p0_: P0})


E00 = sym(sp.expand(res[(0, 0)]).coeff(Qs))
E01 = sym(sp.expand(res[(0, 1)]).coeff(Qps))
E11 = sym(sp.expand(res[(1, 1)]).coeff(Qs))
KG = sym(sp.expand(res["KG"]).coeff(Qs))
NC = sym(sp.expand(res["NC"]).coeff(Qs) / (3 * C * sp.sin(sp.Symbol("chi")) ** 2 * sp.sin(sp.Symbol("th"))))
FF = sym(sp.expand(res["F"]).coeff(Qps))
rhod = 3 * C * (1 + P0 ** 2 / (2 * M2)) / A_ ** 3
PHD_expr = -AD / A_ * PH + (F * P0D - rhod * V) / 2                         # (0 chi) constraint
PH_sol = sp.solve(E00.subs(PHD, PHD_expr), PH)[0]                           # (00) constraint: algebraic Phi
VD_sol = sp.solve(FF, VD)[0]
DDD_sol = sp.solve(NC, DDD)[0].subs(PHD, PHD_expr)
FDD_sol = sp.solve(KG, FDD)[0].subs(PHD, PHD_expr)
args = (A_, AD, ADD, P0, P0D, P0DD, F, FD, DD, V, mu, lam, C, M2, k2)
fPH = sp.lambdify(args, PH_sol, "numpy")
fPHD = sp.lambdify(args + (PH,), PHD_expr, "numpy")
fVD = sp.lambdify(args + (PH,), VD_sol, "numpy")
fDDD = sp.lambdify(args + (PH,), DDD_sol, "numpy")
fFDD = sp.lambdify(args + (PH,), FDD_sol, "numpy")
fE11 = sp.lambdify(args + (PH, PHD, PHDD), E11, "numpy")


class Pert:
    def __init__(self, cyc, n):
        self.c, self.n = cyc, n
        self.k2 = n * n - 1.0
        self.par = (cyc.mu, cyc.lam, cyc.C, cyc.M2)

    def base(self, bg, y):
        a, ad, p0, p0d = bg
        d = self.c.f(0, bg)
        return (a, ad, d[1], p0, p0d, d[3], y[0], y[1], y[2], y[3]) + self.par + (self.k2,)

    def L(self, bg, y):
        base = self.base(bg, y)
        ph = fPH(*base)
        return np.array([y[1], fFDD(*base, ph), fDDD(*base, ph), fVD(*base, ph)]), ph

    def rhs_aug(self, tt, Y):
        bg = Y[:4]
        M = Y[4:].reshape(4, 4)
        Lm = np.array([self.L(bg, np.eye(4)[j])[0] for j in range(4)]).T
        return np.concatenate([self.c.f(tt, bg), (Lm @ M).ravel()])

    def monodromy(self, bg0, T, rtol=1e-11):
        sol = solve_ivp(self.rhs_aug, (0.0, T), np.concatenate([bg0, np.eye(4).ravel()]), method="DOP853", rtol=rtol, atol=1e-13)
        return sol.y[4:, -1].reshape(4, 4), sol


def checks(c, bg0, T, n):
    P = Pert(c, n)
    rng = np.random.default_rng(n)
    y0 = rng.normal(size=4) * 1e-3

    def rhs(tt, Y):
        return np.concatenate([c.f(tt, Y[:4]), P.L(Y[:4], Y[4:])[0]])
    sol = solve_ivp(rhs, (0.0, T), np.concatenate([bg0, y0]), method="DOP853", rtol=1e-13, atol=1e-16, dense_output=True)
    h = 2e-5
    r1, r2 = [], []
    for x in np.linspace(0.05 * T, 0.95 * T, 160):
        def st(z):
            Y = sol.sol(z)
            return P.L(Y[:4], Y[4:])[1], Y
        pc, Y = st(x)
        pp, _ = st(x + h)
        pm, _ = st(x - h)
        base = P.base(Y[:4], Y[4:])
        phd_alg = (pp - pm) / (2 * h)
        phdd_alg = (pp - 2 * pc + pm) / h ** 2
        phd01 = fPHD(*base, pc)
        r1.append(abs(phd_alg - phd01) / max(abs(phd01), 1e-3 * abs(pc) + 1e-18))
        r2.append(abs(fE11(*base, pc, phd01, phdd_alg)) / (abs(pc) * 10 + 1e-18))
    return float(np.max(r1)), float(np.max(r2))


if __name__ == "__main__":
    c = Cycle(0.35)
    p = c.orbit(0.2229, span=0.05)
    bg0 = np.array([c.neck_a(p, 0.0), 0.0, p, 0.0])
    T = c.flow_to_bounce(p, 0.0)[0]
    out = dict(T=T, checks={}, modes=[])
    for n in (5, 12):
        r1, r2 = checks(c, bg0, T, n)
        out["checks"][n] = dict(constraint_propagation=r1, trace_equation=r2)
        print("n=%d: constraint-propagation residual %.2e (relative to Phi'), trace-equation residual %.2e (relative)" % (n, r1, r2), flush=True)
    for n in list(range(3, 13)) + [15, 20, 30, 40, 60]:
        M, sol = Pert(c, n).monodromy(bg0, T)
        ev = np.linalg.eigvals(M)
        rec = dict(n=n, k2=n * n - 1, multipliers_abs=sorted([float(abs(e)) for e in ev], reverse=True), det=float(np.linalg.det(M)), max_abs=float(max(abs(ev))))
        out["modes"].append(rec)
        print("n=%2d: |multipliers| = %s  det %.6f" % (n, ", ".join("%.5g" % z for z in rec["multipliers_abs"]), rec["det"]), flush=True)
    json.dump(out, open(os.path.join(HERE, "phase2e_perturb_floquet.json"), "w"), indent=1)
