# -*- coding: utf-8 -*-
"""Phase 2F: attribute the inhomogeneous growth: scalar-sector-only Floquet exponent (hilltop/tachyonic) vs the full system, at several mu along the fixed-trace family.
Scalar-only: y'' = -3 H y' - (k^2/a^2 + V''(phi0) + rho_*/M^2) y on the background (no gravity back-reaction, no fluid dynamics)."""
import json, os
import numpy as np
from scipy.integrate import solve_ivp
import phase2f_model as F

HERE = os.path.dirname(os.path.abspath(__file__))
fam = json.load(open(os.path.join(HERE, "phase2f_muscan_up_rc0.2_valid.json")))
fam = [r for r in fam if r["mu"] in (fam[0]["mu"],)] if False else fam
pick = [fam[i] for i in (0,) ] + [r for r in fam if 60 < r["mu"] < 70][:1] + [r for r in fam if 120 < r["mu"] < 135][:1] + [r for r in fam if 200 < r["mu"] < 215][:1]
out = []
for r in pick:
    m = F.Model2(r["x"], mu=r["mu"], rc=0.2)
    cd = F.cycle_data(m, r["phi0_over_v"], span=0.04)
    row = dict(mu=r["mu"], T=cd["T"], mu_T=r["mu"] * cd["T"], n={})
    for n in (3, 6, 12, 30, 60):
        k2 = n * n - 1.0
        def rhs(t, Y):
            bg = Y[:4]; a, ad, p, pd = bg; H = ad / a
            d = m.f(t, bg)
            Vpp = 2 * m.a2 + 12 * m.a4 * p ** 2 + 30 * m.a6 * p ** 4
            meff = k2 / a ** 2 + Vpp + (3 * m.C / a ** 3) / m.M2
            M = Y[4:].reshape(2, 2)
            L = np.array([[0, 1], [-meff, -3 * H]])
            return np.concatenate([d, (L @ M).ravel()])
        sol = solve_ivp(rhs, (0, cd["T"]), np.concatenate([cd["bg0"], np.eye(2).ravel()]), method="DOP853", rtol=1e-11, atol=1e-13)
        Mm = sol.y[4:, -1].reshape(2, 2)
        full = F.growth(m, cd, (n,))[n]
        row["n"][str(n)] = dict(scalar_only=float(max(abs(np.linalg.eigvals(Mm)))), full=full)
    out.append(row)
    print("mu=%.1f mu*T=%.1f | " % (r["mu"], row["mu_T"]) + " ".join("n=%s: KG-only %.3g / full %.3g;" % (k, v["scalar_only"], v["full"]) for k, v in row["n"].items()), flush=True)
json.dump(out, open(os.path.join(HERE, "phase2f_attribution.json"), "w"), indent=1)
