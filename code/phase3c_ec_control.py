# -*- coding: utf-8 -*-
"""Phase 3C control: contraction leg (turnaround -> n = f n*) of the closed ball, with the spin term (B > 0) and without it (B = 0, plain dust).  Same C, same window in a; transfer-matrix growth for n = 3, 6, 12."""
import json, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import phase2f_model as F

HERE = os.path.dirname(os.path.abspath(__file__))
C, B = 0.01, 2.0e-7
out = []
for fr in (0.3, 0.5, 0.7):
    a_f = (3 * C * 2 * B / fr) ** (1.0 / 3.0)
    row = dict(f=fr, a_end=float(a_f), cases=[])
    for Bv in (B, 0.0):
        m = F.Model2(0.35); m.C = C; m.K = -Bv; m.gam = 2.0; m.a2 = m.a4 = m.a6 = 0.0; m.M2 = 1e12
        at = C if Bv == 0.0 else brentq(lambda a: C * a ** 3 - a ** 4 - 3 * Bv * C * C, 0.75 * C, C)
        ev = lambda t, y: y[0] - a_f
        ev.terminal = True
        sol = solve_ivp(m.f, (0, 0.2), [at, 0.0, 0.0, 0.0], method="DOP853", rtol=1e-13, atol=1e-16, events=[ev], dense_output=True)
        t_end = float(sol.t_events[0][0])
        cs2_min = None
        for nmode in (3, 6, 12):
            P = F.Pert2(m, nmode)
            Y0 = np.concatenate([[at, 0.0, 0.0, 0.0], np.eye(4).ravel()])
            s = solve_ivp(P.rhs_aug, (0, t_end), Y0, method="DOP853", rtol=1e-10, atol=1e-14)
            M = s.y[4:, -1].reshape(4, 4)
            Mf = M[2:4, 2:4]; cross = float(np.max(np.abs(np.concatenate([M[2:4, 0:2].ravel(), M[0:2, 2:4].ravel()]))) / np.max(np.abs(Mf)))
            sv = np.linalg.svd(Mf, compute_uv=False)
            row["cases"].append(dict(B=Bv, n=nmode, growth=float(sv[0]), ln=float(np.log(sv[0])), det_fluid_block=float(np.linalg.det(Mf)), cross_coupling=cross, duration=t_end))
            print("f=%.1f  B=%.1e  n=%2d: leg %.5f  fluid-block growth %.4g (ln %.2f)  det %.4f  cross-coupling %.1e" % (fr, Bv, nmode, t_end, sv[0], np.log(sv[0]), np.linalg.det(Mf), cross), flush=True)
    out.append(row)
json.dump(out, open(os.path.join(HERE, "phase3c_ec_control.json"), "w"), indent=1)
