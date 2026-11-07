# -*- coding: utf-8 -*-
"""Phase 2F item 3: coupled fluid WITH PRESSURE (polytrope P = K n^gam, gam = 5/3).  Background and perturbation equations re-derived consistently (phase2f_derive.py).
For each K the homogeneous elliptic orbit is continued at fixed neck-map trace (C/Cc adjusted); inhomogeneous growth for n = 3..240; c_s^2 at the neck."""
import json, os, sys
import numpy as np
import phase2f_model as F

HERE = os.path.dirname(os.path.abspath(__file__))
TARGET = -1.704348
GAM = float(sys.argv[1]) if len(sys.argv) > 1 else 5.0 / 3.0
MU = 20.0


def trace_at(K, x, zg):
    m = F.Model2(x, mu=MU, rc=0.2, K=K, gam=GAM)
    for span in (0.04, 0.12, 0.3):
        p = m.orbit(zg, span=span)
        if p is not None:
            J = m.jac_neck(p, 0.0)
            if J is not None:
                return float(np.trace(J)), p / m.v
    return None, None


def hit(K, x0, z0):
    t0, z = trace_at(K, x0, z0)
    if t0 is None: return None
    x1 = x0 + 0.0008 * (1 if t0 > TARGET else -1)
    t1, z1 = trace_at(K, x1, z)
    if t1 is None: return None
    xa, ta, xb, tb, zb = x0, t0, x1, t1, z1
    for it in range(14):
        if abs(tb - TARGET) < 2e-5: return xb, zb
        if tb == ta: return None
        xn = min(max(xb - (tb - TARGET) * (xb - xa) / (tb - ta), xb - 0.01), xb + 0.01)
        tn, zn = trace_at(K, xn, zb)
        if tn is None: return None
        xa, ta, xb, tb, zb = xb, tb, xn, tn, zn
    return (xb, zb) if abs(tb - TARGET) < 5e-3 else None


x, z = 0.35, 0.22294841
out = []
for K in [0.0, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56]:
    r = hit(K, x, z)
    if r is None:
        print("lost at K=%g" % K, flush=True); break
    x, z = r
    m = F.Model2(x, mu=MU, rc=0.2, K=K, gam=GAM)
    cd = F.cycle_data(m, z, span=0.04)
    if cd is None or abs(cd["trace"]) >= 2:
        print("elliptic orbit lost at K=%g" % K, flush=True); break
    a0 = cd["bg0"][0]; n, P, e = m.fluid(a0); A = 1 + cd["phi0"] ** 2 / (2 * m.M2)
    cs2 = GAM * K * n ** (GAM - 1) / (A + GAM * K * n ** (GAM - 1) / (GAM - 1)) if K else 0.0
    gr = F.growth(m, cd, (3, 6, 12, 30, 60, 120, 240))
    rec = dict(K=K, gam=GAM, x=float(x), phi0_over_v=float(z), T=cd["T"], trace=cd["trace"], P_over_rho_dust=float(P / (A * n)), cs2_neck=float(cs2), a_neck=float(a0), growth={str(k): v for k, v in gr.items()})
    out.append(rec)
    print("K=%-5g C/Cc=%.5f P/rho_d=%.3f cs2=%.3f T=%.4f | growth n=3,6,12,30,60,120,240: %s" % (K, x, P / (A * n), cs2, cd["T"], ", ".join("%.3g" % gr[k] for k in (3, 6, 12, 30, 60, 120, 240))), flush=True)
    json.dump(out, open(os.path.join(HERE, "phase2f_pressure.json"), "w"), indent=1)
