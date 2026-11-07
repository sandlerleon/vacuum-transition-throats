# -*- coding: utf-8 -*-
"""Phase 3C: growth of radial perturbations of the spin-fluid ball in the REGULAR window n <= f n*, avoiding the singular point n* = 1/(2B) where rho + P = n (1 - 2 B n) vanishes (c_s^2 diverges,
the momentum equation of the fluid degenerates).  Transfer matrix from the post-bounce time at which n = f n* through the turnaround back to n = f n*.  det = 1 is the consistency check."""
import json, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import phase2f_model as F

HERE = os.path.dirname(os.path.abspath(__file__))
C, B = 0.01, 2.0e-7
m = F.Model2(0.35); m.C = C; m.K = -B; m.gam = 2.0; m.a2 = m.a4 = m.a6 = 0.0; m.M2 = 1e12
ab = brentq(lambda a: C * a ** 3 - a ** 4 - 3 * B * C * C, 1e-9 * C, 0.75 * C)
bg0 = np.array([ab, 0.0, 0.0, 0.0])
T = 0.0317018
sol = solve_ivp(m.f, (0, 0.04), bg0, method="DOP853", rtol=1e-13, atol=1e-16, dense_output=True)
res = dict(C=C, B=B, windows=[])
for fr in (0.3, 0.5, 0.7):
    a_f = (3 * C / (fr / (2 * B))) ** (1.0 / 3.0)          # n = fr n*  <=>  a = (3C B 2 / fr)^(1/3)
    tgrid = np.linspace(1e-9, 0.04, 400001)
    A_ = sol.sol(tgrid)[0]
    i_turn = int(np.argmax(A_))
    t1 = float(tgrid[np.argmax(A_ > a_f)]); t2 = float(tgrid[i_turn + np.argmax(A_[i_turn:] < a_f)])
    row = dict(f=fr, a_window=float(a_f), t1=t1, t2=t2, duration=t2 - t1, modes=[])
    nn = 3 * C / A_ ** 3; cs2 = -2 * B * nn / (1 - 2 * B * nn)
    sel = (tgrid >= t1) & (tgrid <= t2)
    for nmode in (3, 6, 12, 30, 60):
        P = F.Pert2(m, nmode)
        Y0 = np.concatenate([sol.sol(t1), np.eye(4).ravel()])
        s = solve_ivp(P.rhs_aug, (t1, t2), Y0, method="DOP853", rtol=1e-10, atol=1e-14)
        M = s.y[4:, -1].reshape(4, 4)
        sv = np.linalg.svd(M, compute_uv=False)
        k = np.sqrt(nmode ** 2 - 1.0)
        S = float(np.trapezoid((k * np.sqrt(np.maximum(-cs2, 0)) / A_)[sel], tgrid[sel]))
        row["modes"].append(dict(n=nmode, growth=float(sv[0]), ln_growth=float(np.log(sv[0])), wkb_S=S, det=float(np.linalg.det(M))))
        print("f=%.1f n=%2d: window [%.5f, %.5f] growth %.4g (ln %.2f), WKB S %.2f, det %.5f" % (fr, nmode, t1, t2, sv[0], np.log(sv[0]), S, np.linalg.det(M)), flush=True)
    res["windows"].append(row)
json.dump(res, open(os.path.join(HERE, "phase3c_ec_window.json"), "w"), indent=1)
