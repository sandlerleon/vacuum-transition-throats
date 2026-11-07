# -*- coding: utf-8 -*-
"""Phase 2F item 4: reshape the hilltop.  Sextic family V = 3 + a2 phi^2 + a4 phi^4 + a6 phi^6 with the SAME curvature at the top (V'' = -mu^2) and the SAME true-vacuum depth (-q V_f) but a free well position v
(v/v0 = 1 is the quartic).  For each v the elliptic orbit is continued at fixed neck-map trace (C/Cc adjusted); inhomogeneous growth for n = 3..60."""
import json, os, sys
import numpy as np
import phase2f_model as F

HERE = os.path.dirname(os.path.abspath(__file__))
TARGET = -1.704348
MU = 20.0
V0 = np.sqrt(4 * 2 * 3.0 / MU ** 2)


def trace_at(vr, x, zabs):
    m = F.Model2(x, mu=MU, rc=0.2, v=vr * V0)
    for span in (0.04, 0.12, 0.3):
        p = m.orbit(zabs / m.v, span=span)
        if p is not None:
            J = m.jac_neck(p, 0.0)
            if J is not None:
                return float(np.trace(J)), p
    return None, None


def hit(vr, x0, z0):
    t0, z = trace_at(vr, x0, z0)
    if t0 is None: return None
    x1 = x0 + 0.0008 * (1 if t0 > TARGET else -1)
    t1, z1 = trace_at(vr, x1, z)
    if t1 is None: return None
    xa, ta, xb, tb, zb = x0, t0, x1, t1, z1
    for it in range(14):
        if abs(tb - TARGET) < 2e-5: return xb, zb
        if tb == ta: return None
        xn = min(max(xb - (tb - TARGET) * (xb - xa) / (tb - ta), xb - 0.01), xb + 0.01)
        tn, zn = trace_at(vr, xn, zb)
        if tn is None: return None
        xa, ta, xb, tb, zb = xb, tb, xn, tn, zn
    return (xb, zb) if abs(tb - TARGET) < 5e-3 else None


res = {}
for name, seq in (("wider", [round(1.0 * 1.02 ** k, 4) for k in range(0, 56)]), ("narrower", [round(0.98 ** k, 4) for k in range(0, 45)])):
    x, z = 0.35, 0.22294841 * V0
    out = []
    for vr in seq:
        r = hit(vr, x, z)
        if r is None:
            print("%s: lost at v/v0=%g" % (name, vr), flush=True); break
        x, z = r
        m = F.Model2(x, mu=MU, rc=0.2, v=vr * V0)
        cd = F.cycle_data(m, z / m.v, span=0.04)
        if cd is None or abs(cd["trace"]) >= 2: print("%s: elliptic lost at %g" % (name, vr)); break
        gr = F.growth(m, cd, (3, 6, 12, 30, 60))
        rec = dict(v_over_v0=vr, a4=float(m.a4), a6=float(m.a6), x=float(x), phi0=float(z), T=cd["T"], trace=cd["trace"], growth={str(k): v for k, v in gr.items()})
        out.append(rec)
        print("%-8s v/v0=%.2f a4=%.1f a6=%.1f C/Cc=%.4f phi0=%.4f T=%.4f | growth n=3,6,12,30,60: %s" % (name, vr, m.a4, m.a6, x, z, cd["T"], ", ".join("%.3g" % gr[k] for k in (3, 6, 12, 30, 60))), flush=True)
    res[name] = out
    json.dump(res, open(os.path.join(HERE, "phase2f_potential.json"), "w"), indent=1)
