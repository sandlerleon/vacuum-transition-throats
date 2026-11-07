# -*- coding: utf-8 -*-
"""Phase 3C: linear stability of the bouncing spin-fluid (Einstein-Cartan-type) closed universe.  The effective fluid rho = n - B n^2, p = -B n^2 is the polytrope of phase2f_model with gamma = 2 and K = -B
(no scalar field, no coupling; a cosmological constant V = 3 is present in the machinery and is negligible at the scales used, a <= 0.01).  Radial harmonics Q_n, k^2 = n^2 - 1, longitudinal gauge,
equations of phase2f_derive.py.  Monodromy over one bounce-to-bounce cycle; WKB estimate  S(n) = int k sqrt(max(-c_s^2, 0)) / a dt  for the imaginary-sound-speed phase."""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import phase2f_model as F

HERE = os.path.dirname(os.path.abspath(__file__))
C, B = 0.01, 2.0e-7
m = F.Model2(0.35)
m.C = C; m.K = -B; m.gam = 2.0; m.a2 = m.a4 = m.a6 = 0.0; m.M2 = 1e12
f = lambda a: C * a ** 3 - a ** 4 - 3 * B * C * C
ab = brentq(f, 1e-9 * C, 0.75 * C)
bg0 = np.array([ab, 0.0, 0.0, 0.0])
ev = lambda t, y: y[1]
ev.direction = 1
ev.terminal = True
s1 = solve_ivp(m.f, (0, 1e-4), bg0, method="DOP853", rtol=1e-13, atol=1e-16)
s2 = solve_ivp(m.f, (1e-4, 1.0), s1.y[:, -1], method="DOP853", rtol=1e-13, atol=1e-16, events=[ev])
T = float(s2.t_events[0][0])
print("background: bounce a=%.5f, period T=%.5f; with the cosmological constant V=3 present the period differs from the Lambda=0 value 0.03169 by %.1e" % (ab, T, abs(T / 0.0316869 - 1)))
bgT = s2.y_events[0][0]
print("return error of the background: %.1e" % (np.linalg.norm(bgT[:2] - bg0[:2]) / ab))
sol = solve_ivp(m.f, (0, T), bg0, method="DOP853", rtol=1e-13, atol=1e-16, dense_output=True)
tt = np.linspace(0, T, 8001); Y = sol.sol(tt)
a = Y[0]; nn = 3 * C / a ** 3
cs2 = -2 * B * nn / (1 - 2 * B * nn)
frac_imag = float(np.mean(cs2 < 0)); frac_super = float(np.mean(cs2 > 1))
n_star_time = float(np.trapezoid((nn < 1 / (2 * B)).astype(float), tt) / T)
print("fraction of the cycle with c_s^2 < 0: %.3f ; with c_s^2 > 1 (superluminal): %.3f ; max |c_s^2| away from n*: %.3g" % (frac_imag, frac_super, np.max(np.abs(cs2[np.abs(nn - 1 / (2 * B)) > 0.05 / B]))))
res = dict(C=C, B=B, a_bounce=float(ab), T=T, frac_imaginary_cs=frac_imag, frac_superluminal=frac_super, modes=[])
for nmode in (3, 6, 12, 30, 60):
    P = F.Pert2(m, nmode)
    M = P.monodromy(bg0, T, rtol=1e-10)
    ev_ = np.linalg.eigvals(M)
    k = np.sqrt(nmode ** 2 - 1.0)
    S = float(np.trapezoid(k * np.sqrt(np.maximum(-cs2, 0.0)) / a, tt))
    rec = dict(n=nmode, max_multiplier=float(max(abs(ev_))), ln_multiplier=float(np.log(max(abs(ev_)))), wkb_S=S, det=float(np.linalg.det(M)))
    res["modes"].append(rec)
    print("n=%2d: max|multiplier| = %.4g (ln = %.2f), WKB estimate of the imaginary-c_s growth exponent %.2f, det %.5f" % (nmode, max(abs(ev_)), np.log(max(abs(ev_))), S, rec["det"]), flush=True)
json.dump(res, open(os.path.join(HERE, "phase3c_ec_perturb.json"), "w"), indent=1)
