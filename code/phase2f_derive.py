# -*- coding: utf-8 -*-
"""Phase 2F: GENERALISED derivation of the spherically symmetric linear perturbations of the closed cycle, for (a) a polynomial potential V = 3 + a2 phi^2 + a4 phi^4 + a6 phi^6 and (b) a coupled fluid
WITH PRESSURE (polytrope): particle number density n (m0 = 1), energy density rho_f = A(phi) n + K n^gam/(gam-1), pressure P = K n^gam, A = 1 + phi^2/(2 M2).  K = 0 is the dust of Phase 2E.
Field equation  Box phi = V' + A' n ;  fluid:  number conservation  nabla(n u) = 0  and  (rho + P) u^nu nabla_nu u_mu = - h_mu^nu (d_nu P + A' n d_nu phi)  (from total energy-momentum conservation).
Longitudinal gauge ds^2 = -(1+2 eps Psi) dt^2 + a^2 (1-2 eps Phi) (dchi^2 + sin^2 chi dOmega^2).   Output: phase2f_perturb.pkl (reduced mode equations, harmonic basis Q_n, k2 = n^2 - 1).
"""
import pickle
import sympy as sp

t, chi, th, ph, eps = sp.symbols("t chi th ph epsilon")
C, M2, K, gam = sp.symbols("C M2 K gam", positive=True)
a2, a4, a6 = sp.symbols("a2 a4 a6", real=True)
a = sp.Function("a")(t); p0 = sp.Function("p0")(t)
Psi = sp.Function("Psi")(t, chi); Phi = sp.Function("Phi")(t, chi); df = sp.Function("df")(t, chi); dl = sp.Function("dl")(t, chi); U = sp.Function("U")(t, chi)
X = [t, chi, th, ph]
g = sp.diag(-(1 + 2 * eps * Psi), a ** 2 * (1 - 2 * eps * Phi), a ** 2 * (1 - 2 * eps * Phi) * sp.sin(chi) ** 2, a ** 2 * (1 - 2 * eps * Phi) * sp.sin(chi) ** 2 * sp.sin(th) ** 2)
gi = g.inv()
n_ = 4
Gam = [[[sum(gi[i, l] * (sp.diff(g[l, j], X[k]) + sp.diff(g[l, k], X[j]) - sp.diff(g[j, k], X[l])) for l in range(n_)) / 2 for k in range(n_)] for j in range(n_)] for i in range(n_)]
def Ric(j, k):
    r = 0
    for i in range(n_):
        r += sp.diff(Gam[i][j][k], X[i]) - sp.diff(Gam[i][j][i], X[k])
        for l in range(n_):
            r += Gam[i][i][l] * Gam[l][j][k] - Gam[i][k][l] * Gam[l][j][i]
    return r
R_ = {(j, k): Ric(j, k) for j in range(n_) for k in range(j, n_)}
Rs = sum(gi[j, j] * R_[(j, j)] for j in range(n_))
Gt = {(j, k): R_[(j, k)] - g[j, k] * Rs / 2 for (j, k) in R_}

xs = sp.Symbol("x")
Vx = 3 + a2 * xs ** 2 + a4 * xs ** 4 + a6 * xs ** 6
Ax = 1 + xs ** 2 / (2 * M2)
Vf = lambda x: Vx.subs(xs, x)
Af = lambda x: Ax.subs(xs, x)
dV = lambda x: sp.diff(Vx, xs).subs(xs, x)
dA = lambda x: sp.diff(Ax, xs).subs(xs, x)
phi = p0 + eps * df
n0 = 3 * C / a ** 3
nn = n0 * (1 + eps * dl)
rho_f = Af(phi) * nn + K * nn ** gam / (gam - 1)
P_f = K * nn ** gam
dphi = [sp.diff(phi, x) for x in X]
kin = sum(gi[i, i] * dphi[i] ** 2 for i in range(n_)) / 2
u_dn = [-(1 + eps * Psi), eps * U, 0, 0]
T = {}
for (j, k) in Gt:
    T[(j, k)] = dphi[j] * dphi[k] - g[j, k] * (kin + Vf(phi)) + (rho_f + P_f) * u_dn[j] * u_dn[k] + P_f * g[j, k]
lin = lambda e: sp.simplify(sp.diff(e, eps).subs(eps, 0))
E1 = {key: lin(Gt[key] - T[key]) for key in [(0, 0), (0, 1), (1, 1), (2, 2)]}
E0 = {key: sp.simplify((Gt[key] - T[key]).subs(eps, 0)) for key in [(0, 0), (1, 1)]}
sq = sp.sqrt(-g.det())
box = sum(sp.diff(sq * gi[i, i] * dphi[i], X[i]) for i in range(n_)) / sq
KG = box - dV(phi) - dA(phi) * nn
KG1 = lin(KG); KG0 = sp.simplify(KG.subs(eps, 0))
u_up = [gi[0, 0] * u_dn[0], gi[1, 1] * u_dn[1]]
NC = sp.diff(sq * nn * u_up[0], t) + sp.diff(sq * nn * u_up[1], chi)
NC1 = lin(NC)
def cov(nu, mu_):
    return sp.diff(u_dn[mu_], X[nu]) - sum(Gam[l][nu][mu_] * u_dn[l] for l in range(2))
u_up_full = [gi[0, 0] * u_dn[0], gi[1, 1] * u_dn[1], 0, 0]
acc_chi = sum(u_up_full[nu] * cov(nu, 1) for nu in range(2))
h_chi = lambda nu: (1 if nu == 1 else 0) + u_dn[1] * u_up_full[nu]
force = (rho_f + P_f) * acc_chi + sum(h_chi(nu) * (sp.diff(P_f, X[nu]) + dA(phi) * nn * dphi[nu]) for nu in range(2))
F1 = lin(force)
print("derived; reducing")
# --------------------------------------------------------------------------- reduction to the harmonic basis
k2 = sp.Symbol("k2", positive=True)
ph_ = sp.Function("ph")(t); f_ = sp.Function("f")(t); Dd = sp.Function("Dd")(t); v_ = sp.Function("v")(t)
Qf = sp.Function("Q")(chi)
Q, Qp = sp.symbols("Q Qp")
cot = sp.Symbol("cot")
rep = {Phi: ph_ * Qf, Psi: ph_ * Qf, df: f_ * Qf, dl: Dd * Qf, U: v_ * sp.diff(Qf, chi)}
def red(e):
    e = e.subs(rep).doit()
    e = e.subs(sp.tan(chi), 1 / cot)
    e = e.subs(sp.Derivative(Qf, (chi, 3)), sp.diff(-2 * sp.cos(chi) / sp.sin(chi) * sp.diff(Qf, chi) - k2 * Qf, chi))
    e = e.subs(sp.Derivative(Qf, (chi, 2)), -2 * cot * sp.diff(Qf, chi) - k2 * Qf)
    e = e.subs(sp.Derivative(Qf, chi), Qp).subs(Qf, Q)
    return sp.expand(e)
out = {}
for key in [(0, 0), (0, 1), (1, 1)]:
    out[key] = red(E1[key])
out[(2, 2)] = red(sp.simplify(E1[(2, 2)] / sp.sin(chi) ** 2))
out["KG"] = red(KG1); out["F"] = red(F1)
out["NC"] = red(sp.simplify(NC1.subs(sp.sqrt(a ** 6 * sp.sin(chi) ** 4 * sp.sin(th) ** 2), a ** 3 * sp.sin(chi) ** 2 * sp.sin(th))))
out["aniso"] = sp.simplify(out[(2, 2)] - out[(1, 1)])
out["E0"] = E0; out["KG0"] = KG0
print("traceless check (should be 0):", out["aniso"])
pickle.dump(out, open("phase2f_perturb.pkl", "wb"))
print("saved")
