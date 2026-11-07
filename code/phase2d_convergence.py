# -*- coding: utf-8 -*-
"""Phase 2D: convergence and constraint diagnostics for the periodic moving shell (q = 0, m/m_max = 0.10) around the C/Cc = 0.35 cycle."""
import json, os
import numpy as np
from scipy.integrate import solve_ivp
import bridge_moving_shell as MS
from phase2d_shell_search import Interior, setup

HERE = os.path.dirname(os.path.abspath(__file__))
I, yi0, T = setup()
sh = json.load(open(os.path.join(HERE, "phase2d_shells_q0.json")))["shells"][0]
q, w, m, chi0 = sh["q"], sh["w"], sh["m_ext"], sh["chi0_over_pi"] * np.pi
out = dict(chi0_over_pi=sh["chi0_over_pi"], m=m, w=w, convergence=[])
for rtol, ms in ((1e-7, 0.02), (1e-9, 0.01), (1e-10, 0.005), (1e-11, 0.002), (1e-12, 0.002), (1e-13, 0.001)):
    ev = lambda tau, y: y[0] - T
    ev.terminal = True
    sol = solve_ivp(MS.rhs_factory(I, m, q, w), (0.0, 4.0 * T), [0.0, chi0, 0.0] + yi0, method="DOP853", rtol=rtol, atol=rtol * 1e-2, max_step=ms, events=[ev])
    y = sol.y_events[0][0]
    out["convergence"].append(dict(rtol=rtol, max_step=ms, return_dchi=float(y[1] - chi0), return_beta=float(y[2]), return_a=float(y[3] - yi0[0]), return_phi=float(y[5] - yi0[2])))
    print("rtol=%g: return dchi=%.2e beta=%.2e da=%.2e dphi=%.2e" % (rtol, y[1] - chi0, y[2], y[3] - yi0[0], y[5] - yi0[2]), flush=True)
# residuals along the orbit
ev = lambda tau, y: y[0] - T
ev.terminal = True
sol = solve_ivp(MS.rhs_factory(I, m, q, w), (0.0, 4.0 * T), [0.0, chi0, 0.0] + yi0, method="DOP853", rtol=1e-12, atol=1e-14, max_step=0.001, dense_output=True, events=[ev])
tau = np.linspace(0, sol.t_events[0][0], 4001)
Y = sol.sol(tau)
res = dict(friedmann=0.0, gamma_in=0.0, israel_p=0.0, mass=0.0, flux=0.0)
S = [MS.shell_state(I, Y[:, i], m, q, w) for i in range(len(tau))]
sig = np.array([s["sigma"] for s in S]); pr = np.array([s["p"] for s in S]); R = np.array([s["R"] for s in S])
for i, s in enumerate(S):
    t_, chi, beta, a, ad, ph, wd = Y[:7, i]
    res["friedmann"] = max(res["friedmann"], abs(I.c.constraint([a, ad, ph, wd])))
    res["gamma_in"] = max(res["gamma_in"], abs(s["G_in"] ** 2 - (1 + s["Rp"] ** 2 - s["rho"] * s["R"] ** 2 / 3.0)))
    res["israel_p"] = max(res["israel_p"], abs(s["p"] - w * s["sigma"]))
    # mass function: m = R^3 (rho/3 + q) + sigma R^2 Gamma_in - sigma^2 R^3 / 4  must equal the constant exterior mass
    mm = s["R"] ** 3 * (s["rho"] / 3.0 + q) + s["sigma"] * s["R"] ** 2 * s["G_in"] - s["sigma"] ** 2 * s["R"] ** 3 / 4.0
    res["mass"] = max(res["mass"], abs(mm - m) / m)
f = -(R ** 2) * np.array([(s["rho"] + s["P"]) * np.cosh(Y[2, i]) * np.sinh(Y[2, i]) for i, s in enumerate(S)])
lhs = np.gradient(sig * R ** 2, tau) + pr * np.gradient(R ** 2, tau)
res["flux"] = float(np.max(np.abs(lhs[20:-20] - f[20:-20])) / np.max(np.abs(f[20:-20])))
out["residuals"] = {k: float(v) for k, v in res.items()}
print("residuals: Friedmann %.1e | Gamma_in^2 identity %.1e | p - w sigma %.1e | mass function relative %.1e | flux law (relative) %.1e" % tuple(res.values()))
json.dump(out, open(os.path.join(HERE, "phase2d_convergence.json"), "w"), indent=1)
