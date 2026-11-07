# -*- coding: utf-8 -*-
"""Nonlinear survival of the elliptic single-field cold cycle (C/Cc = 0.35, mu = 20, rc = 0.2, q = 1): iterate the neck map from perturbed starts (Phase 1 code)."""
import json, os
import numpy as np
from scipy.optimize import brentq
import closed_cycle as CC, closed_orbits as CO, closed_stability as CS

HERE = os.path.dirname(os.path.abspath(__file__))
m = CC.Model(q=1.0, rc=0.2, mu=20.0, C=0.35 * CC.C_CRIT)
ls0 = np.log10(0.22295)
r = brentq(lambda z: CO.g(m, z), ls0 - 0.01, ls0 + 0.01, xtol=1e-14)
s = 10 ** r * m.v
nm = CS.neck_map(m, s, 0.0)
print("orbit phi0/v=%.8f period %.6f return err %.1e" % (10 ** r, nm[2], abs(nm[0] - s) / m.v), flush=True)
res = dict(phi0_over_v=float(10 ** r), period=nm[2], survival=[])
for rel in (1e-6, 1e-4, 1e-3, 1e-2, 3e-2, 1e-1):
    row = []
    for sd in range(3):
        rng = np.random.default_rng(sd); d = rng.normal(size=2); d /= np.linalg.norm(d)
        st = (s + rel * d[0] * m.v, rel * d[1] * m.v); n = 0; dev = 0.0
        for k in range(200):
            q_ = CS.neck_map(m, st[0], st[1], t_max=10.0)
            if q_ is None:
                break
            n += 1; st = (q_[0], q_[1]); dev = max(dev, float(np.hypot(st[0] - s, st[1]) / m.v))
        row.append((n, dev))
    print(" perturbation %.0e: bounces survived (cap 200) %s, max neck deviation/v %s" % (rel, [x[0] for x in row], ["%.2e" % x[1] for x in row]), flush=True)
    res["survival"].append(dict(rel=rel, bounces=[x[0] for x in row], dev=[x[1] for x in row]))
json.dump(res, open(os.path.join(HERE, "single_elliptic_survival.json"), "w"), indent=1)
