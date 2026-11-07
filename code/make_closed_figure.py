# -*- coding: utf-8 -*-
"""Figure of an exact closed cold cycle and of its instability (closed_cycle.py model; parameters from closed_stability.json)."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import closed_cycle as CC, closed_orbits as CO

HERE = os.path.dirname(os.path.abspath(__file__))
ST = json.load(open(os.path.join(HERE, "closed_stability.json")))
pick = [r for r in ST if r["mu"] == 10.0 and r["rc"] == 0.05 and abs(r["C_over_Ccrit"] - 0.9) < 1e-9][0]
m = CC.Model(q=pick["q"], rc=pick["rc"], mu=pick["mu"], C=pick["C_over_Ccrit"] * CC.C_CRIT)
ls = pick["log10_seed_over_v"]
fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.6))
for d, col, lab in ((0.0, "#1f77b4", "exact seed"), (1e-10, "#d62728", "relative seed error 1e-10")):
    s = 10 ** ls * (1 + d) * m.v
    a0 = CC.neck_a(m, s)
    sol = solve_ivp(m.rhs, (0, 14), [a0, 0, s, 0, 0], method="DOP853", rtol=1e-13, atol=1e-15, max_step=0.01, events=lambda t, y: y[0] - 1e-3)
    sol_ev = lambda t, y: y[0] - 1e-3
    t, a, ad, p, pd = sol.t, sol.y[0], sol.y[1], sol.y[2], sol.y[3]
    rho = np.array([m.rho(a[i], p[i], pd[i], 0.0) for i in range(len(t))]) / 3.0
    ax[0].plot(t, a, color=col, label=lab); ax[1].plot(t, p / m.v, color=col); ax[2].plot(t, rho, color=col)
ax[0].set_yscale("log"); ax[0].set_xlabel("t (R_dS)"); ax[0].set_ylabel("curvature radius a (R_dS)"); ax[0].legend(fontsize=7)
ax[1].set_xlabel("t (R_dS)"); ax[1].set_ylabel("phi / v"); ax[2].set_yscale("symlog", linthresh=1.0); ax[2].set_xlabel("t (R_dS)"); ax[2].set_ylabel("rho / V_f")
for a_ in ax: a_.set_xlim(0, 14)
ax[1].set_title("q = 1, rho_crit = 0.05 V_f, mu = 10 H, C = 0.9 C_c", fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "..", "figures", "FigC1_closed_cycle.png"), dpi=150)
print("ok", pick["period"], pick["a_neck"])
