# -*- coding: utf-8 -*-
import pickle, sympy as sp
D = pickle.load(open("phase2e_perturb_eqs.pkl", "rb"))
t, chi = sp.symbols("t chi")
C, M2, mu, lam = sp.symbols("C M2 mu lam", positive=True)
k2 = sp.Symbol("k2", positive=True)
a = sp.Function("a")(t); p0 = sp.Function("p0")(t)
Psi = sp.Function("Psi")(t, chi); Phi = sp.Function("Phi")(t, chi); df = sp.Function("df")(t, chi); dl = sp.Function("dl")(t, chi); U = sp.Function("U")(t, chi)
ph = sp.Function("ph")(t); f = sp.Function("f")(t); Dd = sp.Function("Dd")(t); v = sp.Function("v")(t)
Qf = sp.Function("Q")(chi)
Q, Qp = sp.symbols("Q Qp")
rep = {Phi: ph * Qf, Psi: ph * Qf, df: f * Qf, dl: Dd * Qf, U: v * sp.diff(Qf, chi)}
def red(e):
    e = e.subs(rep).doit()
    e = e.subs(sp.tan(chi), 1 / sp.Symbol("cot"))     # tan -> 1/cot
    e = e.subs(sp.Derivative(Qf, (chi, 3)), sp.diff(-2 * sp.cos(chi) / sp.sin(chi) * sp.diff(Qf, chi) - k2 * Qf, chi))
    e = e.subs(sp.Derivative(Qf, (chi, 2)), -2 * sp.Symbol("cot") * sp.diff(Qf, chi) - k2 * Qf)
    e = e.subs(sp.Derivative(Qf, chi), Qp).subs(Qf, Q)
    return sp.expand(e)
cot = sp.Symbol("cot")
res = {}
for key in [(0, 0), (0, 1), (1, 1), (2, 2)]:
    e = D["E"][key]
    if key == (2, 2): e = e / sp.sin(chi) ** 2
    e = sp.simplify(e)
    res[key] = red(e)
res["KG"] = red(D["KG1"]); res["NC"] = red(sp.simplify(D["NC1"].subs(sp.sqrt(a ** 6 * sp.sin(chi) ** 4 * sp.sin(sp.Symbol("th")) ** 2), a ** 3 * sp.sin(chi) ** 2 * sp.sin(sp.Symbol("th"))))); res["F"] = red(D["F1"])
for k, e in res.items():
    print(k, ":", sp.collect(e, [Q, Qp, cot]))
    print()
# anisotropy check: (2,2)/sin^2 - (1,1)
an = sp.simplify(res[(2, 2)] - res[(1, 1)])
print("traceless (theta-theta/sin^2 minus chi-chi) with Phi=Psi:", an)
pickle.dump(res, open("phase2e_perturb_reduced.pkl", "wb"))
