# -*- coding: utf-8 -*-
"""Derive (sympy) the first-order equations for spherically symmetric perturbations of the closed FRW cold cycle with scalar + dust coupled by the mass function A(phi), in longitudinal gauge.
Units 8 pi G = 1, G_{mu nu} = T_{mu nu}.  Output: pickled linearised equations in the harmonic basis Q_n(chi) = sin(n chi)/sin(chi), k^2 = n^2 - 1."""
import pickle, sympy as sp

t, chi, th, ph, eps = sp.symbols("t chi th ph epsilon")
C, M2, mu, lam = sp.symbols("C M2 mu lam", positive=True)
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

phi = p0 + eps * df
Q_pot = lambda x: 1 + 0 * x
Vf = lambda x: 3 - mu ** 2 * x ** 2 / 2 + lam * x ** 4 / 4
Af = lambda x: 1 + x ** 2 / (2 * M2)
rs0 = 3 * C / a ** 3
rstar = rs0 * (1 + eps * dl)
dphi = [sp.diff(phi, x) for x in X]
kin = sum(gi[i, i] * dphi[i] ** 2 for i in range(n_)) / 2
u_dn = [-(1 + eps * Psi), eps * U, 0, 0]            # u_t, u_chi (lower), normalised to O(eps)
T = {}
for (j, k) in Gt:
    T[(j, k)] = dphi[j] * dphi[k] - g[j, k] * (kin + Vf(phi)) + Af(phi) * rstar * u_dn[j] * u_dn[k]
E = {key: sp.simplify(sp.diff(Gt[key] - T[key], eps).subs(eps, 0)) for key in [(0, 0), (0, 1), (1, 1), (2, 2)]}
E0 = {key: sp.simplify((Gt[key] - T[key]).subs(eps, 0)) for key in [(0, 0), (0, 1), (1, 1), (2, 2)]}
print("background eqs:", E0)
# scalar field: box phi - V' - A' rstar = 0
sq = sp.sqrt(-g.det())
box = sum(sp.diff(sq * gi[i, i] * dphi[i], X[i]) for i in range(n_)) / sq
KG = box - sp.diff(Vf(sp.Symbol("x")), sp.Symbol("x")).subs(sp.Symbol("x"), phi) - sp.diff(Af(sp.Symbol("x")), sp.Symbol("x")).subs(sp.Symbol("x"), phi) * rstar
KG1 = sp.simplify(sp.diff(KG, eps).subs(eps, 0)); KG0 = sp.simplify(KG.subs(eps, 0))
# dust: number conservation d_mu (sqrt(-g) rstar u^mu) = 0, u^mu = g^{mu nu} u_nu
u_up = [gi[0, 0] * u_dn[0], gi[1, 1] * u_dn[1]]
NC = sp.diff(sq * rstar * u_up[0], t) + sp.diff(sq * rstar * u_up[1], chi)
NC1 = sp.simplify(sp.diff(NC, eps).subs(eps, 0)); NC0 = sp.simplify(NC.subs(eps, 0))
# dust force (chi component): u^nu nabla_nu u_chi = - h_chi^nu d_nu ln A ;  nabla_nu u_chi = d_nu u_chi - Gamma^l_{nu chi} u_l
def cov(nu, mu_):
    return sp.diff(u_dn[mu_], X[nu]) - sum(Gam[l][nu][mu_] * u_dn[l] for l in range(2))
u_up_full = [gi[0, 0] * u_dn[0], gi[1, 1] * u_dn[1], 0, 0]
acc_chi = sum(u_up_full[nu] * cov(nu, 1) for nu in range(2))
lnA = sp.log(Af(phi))
h_chi = lambda nu: sum(0 for _ in [0]) + (1 if nu == 1 else 0) + u_dn[1] * u_up_full[nu]
force = acc_chi + sum(h_chi(nu) * sp.diff(lnA, X[nu]) for nu in range(2))
F1 = sp.simplify(sp.diff(force, eps).subs(eps, 0))
pickle.dump(dict(E=E, E0=E0, KG1=KG1, KG0=KG0, NC1=NC1, NC0=NC0, F1=F1), open("phase2e_perturb_eqs.pkl", "wb"))
for k, v in E.items(): print("E1", k, ":", v)
print("KG1:", KG1); print("NC1:", NC1); print("F1:", F1)
print("KG0:", KG0); print("NC0:", NC0)
