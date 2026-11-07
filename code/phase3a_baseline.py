# -*- coding: utf-8 -*-
"""Phase 3A: classical collapse baseline.  Units 8 pi G = 1 for the scalar solver (as in dnull.py), G = 1 for the exact dust solutions (stated where used).

(1) Einstein-scalar collapse with a regular centre (dnull.py, double-null): stored 72-run scan (collapse_scan.json, validation_extra.json) + a fresh run with a resolution check.
(2) Oppenheimer-Snyder (closed dust ball + Schwarzschild exterior): numerical ODE against the exact cycloid solution; Misner-Sharp mass, trapped-surface time, density divergence.
(3) Inhomogeneous LTB dust collapse (exact parametric solution): trapped-surface formation, singularity times, shell crossing, absence of any re-expanding region.
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import dnull

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}

# ------------------------------------------------------------------ (1) scalar collapse
VE = json.load(open(os.path.join(HERE, "validation_extra.json")))
SC = json.load(open(os.path.join(HERE, "collapse_scan.json")))
ended = {}
for r in SC:
    ended[r["ended_by"]] = ended.get(r["ended_by"], 0) + 1
out["scalar_scan"] = dict(n_runs=len(SC), n_trapped=sum(1 for r in SC if r["first_trapped"] is not None), creation=VE["creation"], ended_by=ended)
print("scalar scan: %d runs, trapped surface formed in %d; initially normal core %d, anti-trapped region created later in %d; ended by %s" %
      (len(SC), out["scalar_scan"]["n_trapped"], VE["creation"]["n_initially_normal_core"], VE["creation"]["n_created_later"], ended))


def fresh(h):
    o = dnull.run(A=0.05, V_w=0.8, w=0.05, h=h, snap_every=2)
    dg = o["diag"]
    ft = int(np.argmax(dg[:, 5] > 0)) if dg[:, 5].max() > 0 else None
    uft = float(dg[ft, 0])
    ks = [k for k in range(len(o["U"])) if 0 < o["U"][k] <= uft - 0.05]
    k = ks[-1]
    return dict(h=h, rows=int(dg.shape[0]), u_first_trapped=uft, trapped_points_max=int(dg[:, 5].max()), anti_trapped_points_max=int(dg[:, 6].max()),
                u_snapshot=float(o["U"][k]), min_rv_before_trap=float(dg[:ft, 1].min()), max_constraint_Cuu=float(dg[:, 10].max()))


out["scalar_fresh"] = [fresh(0.0125), fresh(0.00625)]
for r in out["scalar_fresh"]:
    print("fresh h=%.5f: trapped from u=%.5f, anti-trapped points max %d, min r_v before trapping %.4f (snapshot u=%.3f)" %
          (r["h"], r["u_first_trapped"], r["anti_trapped_points_max"], r["min_rv_before_trap"], r["u_snapshot"]))

# ------------------------------------------------------------------ (2) Oppenheimer-Snyder (G = 1, c = 1)
M, chi0 = 1.0, 0.9                       # exterior mass M, boundary angle chi0 (boundary radius R = a sin chi0 ; M = (1/2) a_max sin^3 chi0)
amax = 2.0 * M / np.sin(chi0) ** 3
rho_of = lambda a: 3.0 * amax / (8.0 * np.pi * a ** 3)            # a'^2 = amax/a - 1 ; rho = 3 amax / (8 pi a^3)
fo = lambda t, y: [y[1], -amax / (2.0 * y[0] ** 2)]                 # a'' = -amax / (2 a^2)
ev = lambda t, y: y[0] - 1e-4
ev.terminal = True
sol = solve_ivp(fo, (0, 5 * amax), [amax, 0.0], method="DOP853", rtol=1e-13, atol=1e-15, events=ev, dense_output=True)
t_s_num = sol.t_events[0][0]
eta = np.linspace(0, np.pi, 2001)
a_exact = 0.5 * amax * (1 + np.cos(eta)); t_exact = 0.5 * amax * (eta + np.sin(eta))
msk = t_exact < t_s_num * 0.999
err_a = float(np.max(np.abs(sol.sol(t_exact[msk])[0] - a_exact[msk])))
a_h = 2.0 * M / np.sin(chi0)
t_h_num = brentq(lambda t: sol.sol(t)[0] * np.sin(chi0) - 2 * M, 1e-9, t_s_num * 0.999)
eta_h = np.arccos(2 * a_h / amax - 1); t_h_exact = 0.5 * amax * (eta_h + np.sin(eta_h))
tt = np.linspace(0, t_s_num * 0.999, 2000); Y = sol.sol(tt)
m_ms = 4.0 * np.pi / 3.0 * rho_of(Y[0]) * (Y[0] * np.sin(chi0)) ** 3
fried = np.abs(Y[1] ** 2 - (amax / Y[0] - 1.0))
out["oppenheimer_snyder"] = dict(M=M, chi0=chi0, a_max=amax, t_singular_num=float(t_s_num), t_singular_exact=float(np.pi * amax / 2), a_error_vs_exact=err_a,
                                 t_horizon_num=float(t_h_num), t_horizon_exact=float(t_h_exact), mass_spread=float(np.ptp(m_ms) / M), friedmann_residual=float(fried.max()),
                                 density_at_end=float(rho_of(1e-4)), min_a=float(Y[0].min()), reexpands=bool(np.any(np.diff(Y[0]) > 0)))
print("OS: singularity t=%.6f (exact %.6f), horizon t=%.6f (exact %.6f), a error %.1e, mass spread %.1e, Friedmann residual %.1e, re-expands: %s" %
      (t_s_num, np.pi * amax / 2, t_h_num, t_h_exact, err_a, np.ptp(m_ms) / M, fried.max(), out["oppenheimer_snyder"]["reexpands"]))

# ------------------------------------------------------------------ (3) LTB dust collapse (G = 1): R_dot^2 = 2 M(x)/R + 2 E(x),  E < 0
M0, m1, k0, e1 = 1.0, 0.4, 0.9, -0.3
Mf = lambda x: M0 * x ** 3 * (1 + m1 * x ** 2)
Ef = lambda x: -0.5 * k0 * x ** 2 * (1 + e1 * x ** 2)
dM = lambda x: M0 * (3 * x ** 2 + 5 * m1 * x ** 4)
xs = np.linspace(0.05, 1.0, 96)


def R_of(t, x):
    """R at time t after the moment of maximum expansion (R_dot = 0), collapsing branch: R = (M/b)(1 - cos eta), t + t_m = (M/b^{3/2})(eta - sin eta), eta in (pi, 2 pi), b = -2E."""
    Mx, Ex = Mf(x), Ef(x); b = -2 * Ex
    tc = 2 * np.pi * Mx / b ** 1.5; tm = tc / 2.0
    if t >= tc - tm:
        return 0.0
    eta_ = brentq(lambda e: Mx / b ** 1.5 * (e - np.sin(e)) - (t + tm), np.pi - 1e-9, 2 * np.pi - 1e-12)
    return Mx / b * (1 - np.cos(eta_))


res = dict(t_collapse=[], t_horizon=[], trapped=[])
for x in xs:
    Mx, Ex = Mf(x), Ef(x); b = -2 * Ex
    tc = 2 * np.pi * Mx / b ** 1.5; tmax = tc / 2
    th = None
    if 2 * Mx / b > 2 * Mx:
        eta_h = brentq(lambda e: Mx / b * (1 - np.cos(e)) - 2 * Mx, np.pi, 2 * np.pi - 1e-9)
        th = Mx / b ** 1.5 * (eta_h - np.sin(eta_h)) - tmax
    res["t_collapse"].append(float(tc - tmax)); res["t_horizon"].append(th); res["trapped"].append(th is not None and th < tc - tmax)
t_eval = np.linspace(0, 0.98 * min(res["t_collapse"]), 60)
mins_dR, maxs_rho = [], []
for t in t_eval:
    Rr = np.array([R_of(t, x) for x in xs]); dR = np.gradient(Rr, xs)
    mins_dR.append(float(dR.min())); maxs_rho.append(float(np.max(dM(xs) / (4 * np.pi * Rr ** 2 * dR))))
x0 = 0.5; tc0 = 2 * np.pi * Mf(x0) / (-2 * Ef(x0)) ** 1.5
tq = np.linspace(0, 0.9 * tc0 / 2, 400)
Rq = np.array([R_of(t, x0) for t in tq]); dt = tq[1] - tq[0]
Rdot = np.gradient(Rq, dt); Rddot = np.gradient(Rdot, dt)
eom = float(np.max(np.abs(Rddot[5:-5] + Mf(x0) / Rq[5:-5] ** 2) / (Mf(x0) / Rq[5:-5] ** 2)))
fi = float(np.max(np.abs(Rdot[5:-5] ** 2 - (2 * Mf(x0) / Rq[5:-5] + 2 * Ef(x0)))))
out["ltb"] = dict(n_shells=len(xs), first_singularity=float(min(res["t_collapse"])), last_singularity=float(max(res["t_collapse"])), central_first=bool(np.argmin(res["t_collapse"]) == 0),
                  n_trapped=int(sum(res["trapped"])), min_dR_dx=float(min(mins_dR)), shell_crossing_before_collapse=bool(min(mins_dR) <= 0),
                  max_density_grid=float(max(maxs_rho)), eom_residual=eom, first_integral_residual=fi)
print("LTB: %d shells, singular times %.4f..%.4f (centre first: %s), trapped on %d/%d shells, shell crossing before collapse: %s, EOM residual %.1e, first integral %.1e" %
      (len(xs), min(res["t_collapse"]), max(res["t_collapse"]), out["ltb"]["central_first"], out["ltb"]["n_trapped"], len(xs), out["ltb"]["shell_crossing_before_collapse"], eom, fi))
json.dump(out, open(os.path.join(HERE, "phase3a_baseline.json"), "w"), indent=1)
