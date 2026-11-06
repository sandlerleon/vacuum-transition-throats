# -*- coding: utf-8 -*-
"""Symbolic derivation of the spherically symmetric Einstein-scalar equations in double-null gauge,
    ds^2 = -f(u,v) du dv + r(u,v)^2 dOmega^2,      G_ab = T_ab,   T_ab = d_a phi d_b phi - g_ab [ (1/2)(d phi)^2 + U(phi) ]   (8 pi G = 1).
Prints the components needed by dnull.py. Run once; the result is hard-coded and checked numerically in dnull.py (constraint convergence)."""
import sympy as sp

u, v, th, ph = sp.symbols("u v theta phi_ang")
f = sp.Function("f")(u, v)
r = sp.Function("r")(u, v)
p = sp.Function("p")(u, v)          # scalar field
U = sp.Function("U")
x = [u, v, th, ph]
g = sp.Matrix([[0, -f / 2, 0, 0], [-f / 2, 0, 0, 0], [0, 0, r ** 2, 0], [0, 0, 0, r ** 2 * sp.sin(th) ** 2]])
gi = g.inv()


def Gam(a, b, c):
    return sum(gi[a, d] * (sp.diff(g[d, b], x[c]) + sp.diff(g[d, c], x[b]) - sp.diff(g[b, c], x[d])) for d in range(4)) / 2


G = [[[sp.simplify(Gam(a, b, c)) for c in range(4)] for b in range(4)] for a in range(4)]


def Ric(b, c):
    s = 0
    for a in range(4):
        s += sp.diff(G[a][b][c], x[a]) - sp.diff(G[a][b][a], x[c])
        for d in range(4):
            s += G[a][a][d] * G[d][b][c] - G[a][c][d] * G[d][b][a]
    return sp.simplify(s)


R = sp.Matrix(4, 4, lambda b, c: Ric(b, c))
Rs = sp.simplify(sum(gi[a, b] * R[a, b] for a in range(4) for b in range(4)))
E = sp.simplify(R - Rs * g / 2)
dp = [sp.diff(p, xx) for xx in x]
kin = sum(gi[a, b] * dp[a] * dp[b] for a in range(4) for b in range(4))
T = sp.Matrix(4, 4, lambda a, b: dp[a] * dp[b] - g[a, b] * (kin / 2 + U(p)))
for name, (i, j) in {"uu": (0, 0), "uv": (0, 1), "thth": (2, 2)}.items():
    print(name, "  G - T =", sp.simplify(E[i, j] - T[i, j]))
box = sp.simplify(sum(gi[a, b] * (sp.diff(p, x[a], x[b]) - sum(G[c][a][b] * sp.diff(p, x[c]) for c in range(4))) for a in range(4) for b in range(4)))
print("box phi =", box)
