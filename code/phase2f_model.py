# -*- coding: utf-8 -*-
"""Phase 2F: generalised closed-cycle model (polynomial potential, coupled polytropic fluid with pressure) and its inhomogeneous perturbation Floquet analysis.

V(phi) = 3 + a2 phi^2 + a4 phi^4 + a6 phi^6;  coupled fluid: n = 3C/a^3 (m0 = 1), rho_f = A n + K n^gam/(gam-1), P_f = K n^gam, A = 1 + phi^2/(2 M2).  K = 0, a6 = 0 is the Phase 2D/2E model.
Background:  a'' = -(a/6)(2 w^2 - 2V + A n + e + 3 P_f),  phi'' = -3 H phi' - V' - n phi/M2,  a'^2 = a^2 rho/3 - 1,  rho = w^2/2 + V + A n + e,  e = P_f/(gam-1).
"""
import os
import pickle

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from phase2d_cycle_validation import Cycle, VF, CC

HERE = os.path.dirname(os.path.abspath(__file__))


class Model2(Cycle):
    def __init__(self, x, mu=20.0, rc=0.2, q=1.0, a6=0.0, v=None, K=0.0, gam=5.0 / 3.0):
        super().__init__(x, q=q, rc=rc, mu=mu)
        self.K, self.gam = K, gam
        self.a2 = -0.5 * mu ** 2
        if a6 == 0.0 and v is None:
            self.a4 = 0.25 * mu ** 2 / self.v2
            self.a6 = 0.0
        else:
            # sextic family with the same curvature at the top (-mu^2) and the same true-vacuum depth -q V_f at |phi| = v (free): solve a4, a6
            vv = self.v if v is None else v
            self.v2, self.v = vv ** 2, vv
            A_ = np.array([[vv ** 4, vv ** 6], [4 * vv ** 4, 6 * vv ** 6]])
            rhs_ = np.array([-(1 + q) * VF - self.a2 * vv ** 2, -2 * self.a2 * vv ** 2])
            self.a4, self.a6 = np.linalg.solve(A_, rhs_)
            if a6 != 0.0:
                self.a6 = a6
        self.lam = 4 * self.a4

    def V(self, p):
        return VF + self.a2 * p ** 2 + self.a4 * p ** 4 + self.a6 * p ** 6

    def Vp(self, p):
        return 2 * self.a2 * p + 4 * self.a4 * p ** 3 + 6 * self.a6 * p ** 5

    def fluid(self, a):
        n = 3.0 * self.C / a ** 3
        P = self.K * n ** self.gam
        return n, P, (P / (self.gam - 1.0) if self.K else 0.0)

    def rho(self, a, p, w):
        n, P, e = self.fluid(a)
        return 0.5 * w * w + self.V(p) + (1.0 + p * p / (2.0 * self.M2)) * n + e

    def f(self, t, y):
        a, b, p, w = y
        n, P, e = self.fluid(a)
        A = 1.0 + p * p / (2.0 * self.M2)
        return [b, -(a / 6.0) * (2.0 * w * w - 2.0 * self.V(p) + A * n + e + 3.0 * P), w, -3.0 * b / a * w - self.Vp(p) - n * p / self.M2]

    def neck_a(self, p, w):
        g = lambda a: a * a * self.rho(a, p, w) / 3.0 - 1.0
        lo, hi = (self.C / 2.0) ** (1.0 / 3.0) * 1.0001, 10.0
        return brentq(g, lo, hi, xtol=1e-15) if g(lo) * g(hi) < 0 else None


# ----------------------------------------------------------------------------------------- perturbations
_S = None


def _build():
    global _S
    if _S is not None:
        return _S
    res = pickle.load(open(os.path.join(HERE, "phase2f_perturb.pkl"), "rb"))
    t = sp.Symbol("t")
    C, M2, K, gam = sp.symbols("C M2 K gam", positive=True)
    a2, a4, a6 = sp.symbols("a2 a4 a6", real=True)
    k2 = sp.Symbol("k2", positive=True)
    a_ = sp.Function("a")(t); p0_ = sp.Function("p0")(t); ph_ = sp.Function("ph")(t); f_ = sp.Function("f")(t); D_ = sp.Function("Dd")(t); v_ = sp.Function("v")(t)
    A_, AD, ADD, P0, P0D, P0DD, PH, PHD, PHDD, F, FD, FDD, DD, DDD, V, VD = sp.symbols("A AD ADD P0 P0D P0DD PH PHD PHDD F FD FDD DD DDD V VD")
    Qs, Qps = sp.symbols("Q Qp")

    def sym(e):
        e = e.subs({sp.Derivative(ph_, (t, 2)): PHDD, sp.Derivative(f_, (t, 2)): FDD, sp.Derivative(a_, (t, 2)): ADD, sp.Derivative(p0_, (t, 2)): P0DD})
        e = e.subs({sp.Derivative(ph_, t): PHD, sp.Derivative(f_, t): FD, sp.Derivative(a_, t): AD, sp.Derivative(p0_, t): P0D, sp.Derivative(D_, t): DDD, sp.Derivative(v_, t): VD})
        return e.subs({ph_: PH, f_: F, D_: DD, v_: V, a_: A_, p0_: P0})
    E00 = sym(sp.expand(res[(0, 0)]).coeff(Qs)); E01 = sym(sp.expand(res[(0, 1)]).coeff(Qps)); E11 = sym(sp.expand(res[(1, 1)]).coeff(Qs))
    KG = sym(sp.expand(res["KG"]).coeff(Qs)); FF = sym(sp.expand(res["F"]).coeff(Qps))
    NC = sym(sp.expand(res["NC"]).coeff(Qs) / (3 * C * sp.sin(sp.Symbol("chi")) ** 2 * sp.sin(sp.Symbol("th"))))
    PHD_expr = sp.solve(E01, PHD)[0]
    PH_sol = sp.solve(E00.subs(PHD, PHD_expr), PH)[0]
    VD_sol = sp.solve(FF, VD)[0].subs(PHD, PHD_expr)
    DDD_sol = sp.solve(NC, DDD)[0].subs(PHD, PHD_expr)
    FDD_sol = sp.solve(KG, FDD)[0].subs(PHD, PHD_expr)
    args = (A_, AD, ADD, P0, P0D, P0DD, F, FD, DD, V, a2, a4, a6, C, M2, K, gam, k2)
    L = lambda e, extra=(): sp.lambdify(args + extra, e, "numpy")
    _S = dict(fPH=L(PH_sol), fPHD=L(PHD_expr, (PH,)), fVD=L(VD_sol, (PH,)), fDDD=L(DDD_sol, (PH,)), fFDD=L(FDD_sol, (PH,)), fE11=L(E11, (PH, PHD, PHDD)), args=args)
    return _S


class Pert2:
    def __init__(self, m, n):
        self.m, self.n = m, n
        self.k2 = n * n - 1.0
        self.S = _build()
        self.par = (m.a2, m.a4, m.a6, m.C, m.M2, m.K, m.gam)

    def base(self, bg, y):
        d = self.m.f(0, bg)
        return (bg[0], bg[1], d[1], bg[2], bg[3], d[3], y[0], y[1], y[2], y[3]) + self.par + (self.k2,)

    def L(self, bg, y):
        base = self.base(bg, y)
        ph = self.S["fPH"](*base)
        return np.array([y[1], self.S["fFDD"](*base, ph), self.S["fDDD"](*base, ph), self.S["fVD"](*base, ph)]), ph

    def rhs_aug(self, tt, Y):
        bg = Y[:4]
        M = Y[4:].reshape(4, 4)
        Lm = np.array([self.L(bg, np.eye(4)[j])[0] for j in range(4)]).T
        return np.concatenate([self.m.f(tt, bg), (Lm @ M).ravel()])

    def monodromy(self, bg0, T, rtol=1e-11):
        sol = solve_ivp(self.rhs_aug, (0.0, T), np.concatenate([bg0, np.eye(4).ravel()]), method="DOP853", rtol=rtol, atol=1e-13)
        return sol.y[4:, -1].reshape(4, 4)


def cycle_data(m, guess_over_v, span=0.05):
    """Orbit, period, background start, homogeneous neck-map trace/det."""
    p = m.orbit(guess_over_v, span=span)
    if p is None:
        return None
    J = m.jac_neck(p, 0.0)
    if J is None:
        return None
    T = m.flow_to_bounce(p, 0.0)[0]
    return dict(phi0=p, phi0_over_v=p / m.v, T=T, bg0=np.array([m.neck_a(p, 0.0), 0.0, p, 0.0]), trace=float(np.trace(J)), det=float(np.linalg.det(J)))


def growth(m, cd, ns=(3, 6, 12, 30, 60)):
    out = {}
    for n in ns:
        M = Pert2(m, n).monodromy(cd["bg0"], cd["T"])
        out[n] = float(max(abs(np.linalg.eigvals(M))))
    return out
