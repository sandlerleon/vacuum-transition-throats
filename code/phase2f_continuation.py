# -*- coding: utf-8 -*-
"""Phase 2F item 1-2: continuation of the homogeneous elliptic orbit in the dust-scalar coupling strength rc = rho_crit/V_f (M^2 = rc V_f / mu^2) at FIXED neck-map trace (so the orbit stays elliptic):
at each rc solve phi0 (turnaround condition) and adjust C/Cc by secant to hold the trace; at every point the inhomogeneous growth for radial harmonics n = 3..60.
python phase2f_continuation.py <up|down> [mu]"""
import json, os, sys
import numpy as np
import phase2f_model as F

HERE = os.path.dirname(os.path.abspath(__file__))
direction = sys.argv[1]; mu = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0
TARGET = -1.704348
fac = 1.06 if direction == "up" else 1 / 1.06


def trace_at(rc, x, zg):
    m = F.Model2(x, mu=mu, rc=rc)
    for span in (0.04, 0.12, 0.3):
        p = m.orbit(zg, span=span)
        if p is not None:
            J = m.jac_neck(p, 0.0)
            if J is not None:
                return float(np.trace(J)), p / m.v
    return None, None


def hit(rc, x0, z0):
    t0, z = trace_at(rc, x0, z0)
    if t0 is None: return None
    x1 = x0 + 0.0008 * (1 if t0 > TARGET else -1)
    t1, z1 = trace_at(rc, x1, z)
    if t1 is None: return None
    xa, ta, xb, tb, zb = x0, t0, x1, t1, z1
    for it in range(12):
        if abs(tb - TARGET) < 2e-5: return xb, zb
        if tb == ta: return None
        xn = xb - (tb - TARGET) * (xb - xa) / (tb - ta)
        xn = min(max(xn, xb - 0.01), xb + 0.01)
        tn, zn = trace_at(rc, xn, zb)
        if tn is None: return None
        xa, ta, xb, tb, zb = xb, tb, xn, tn, zn
    return (xb, zb) if abs(tb - TARGET) < 5e-3 else None

rc, x, z = 0.2, 0.35, 0.22294841
out = []; hist = []
for step in range(80):
    if step > 0:
        rcn = rc * fac
        xg = x + 1.4 * (rcn - rc)
        if len(hist) >= 2:
            (r1, x1), (r2, x2) = hist[-2], hist[-1]
            xg = x2 + (np.log(rcn) - np.log(r2)) / (np.log(r2) - np.log(r1)) * (x2 - x1)
        r = hit(rcn, xg, z) or hit(rcn, x, z)
        if r is None:
            print("continuation lost at rc=%.4f" % rcn, flush=True); break
        rc, (x, z) = rcn, r
    hist.append((rc, x))
    m = F.Model2(x, mu=mu, rc=rc)
    cd = F.cycle_data(m, z, span=0.04)
    if cd is None or abs(cd["trace"]) >= 2:
        print("orbit/elliptic lost at rc=%.4f" % rc, flush=True); break
    gr = F.growth(m, cd, (3, 6, 12, 30, 60))
    g2 = 6 * m.C * cd["phi0"] ** 2 / (m.M2 * cd["bg0"][0] ** 3 * (2 * m.M2 + cd["phi0"] ** 2))
    rec = dict(rc=rc, x=float(x), phi0_over_v=float(z), T=cd["T"], trace=cd["trace"], growth={str(k): v for k, v in gr.items()}, fifth_rate_neck=float(np.sqrt(g2)), a_neck=float(cd["bg0"][0]))
    out.append(rec)
    print("rc=%.4f C/Cc=%.5f phi0/v=%.4f T=%.4f trace=%+.4f | growth n=3,6,12,30,60: %s | g_neck=%.2f" % (rc, x, z, cd["T"], cd["trace"], ", ".join("%.3g" % gr[k] for k in (3, 6, 12, 30, 60)), np.sqrt(g2)), flush=True)
    json.dump(out, open(os.path.join(HERE, "phase2f_continuation_%s_mu%g.json" % (direction, mu)), "w"), indent=1)
