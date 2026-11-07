# -*- coding: utf-8 -*-
"""Phase 2D step 2b: robust test for a local 2-cycle near the multiplier = -1 point.  d(z) = P^2(z) - z around the fixed point; winding number of d on circles of radius r.
Index of the fixed point itself: +1 (elliptic or inverse-reflection...), -1 (hyperbolic).  A change of winding between radii reveals other zeros of d (period-2 points) inside."""
import json, os
import numpy as np
from phase2d_cycle_validation import Cycle

HERE = os.path.dirname(os.path.abspath(__file__))
xs = 0.34584411
out = []
for dx in (-0.004, -0.001, 0.001, 0.004, 0.01):
    x = xs + dx
    c = Cycle(x)
    g = 0.19931 + (0.22295 - 0.19931) / 0.005 * (x - 0.345)
    p = c.orbit(g, span=0.2)
    J = c.jac_neck(p, 0.0)
    rec = dict(x=float(x), dx=dx, trace=float(np.trace(J)), windings={})
    for r in (0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3):
        th = np.linspace(0, 2 * np.pi, 97)[:-1]
        d = []
        ok = True
        for t in th:
            z = r * np.array([np.cos(t), np.sin(t)])
            r1 = c.P(p + z[0] * c.v, z[1] * c.v)
            r2 = None if r1 is None else c.P(r1[0], r1[1])
            if r2 is None:
                ok = False; break
            d.append([(r2[0] - p) / c.v - z[0], r2[1] / c.v - z[1]])
        if not ok:
            rec["windings"][str(r)] = None; continue
        d = np.array(d); ang = np.unwrap(np.arctan2(d[:, 1], d[:, 0]))
        wn = (ang[-1] - ang[0] + (np.arctan2(d[0, 1], d[0, 0]) - np.arctan2(d[-1, 1], d[-1, 0]) if False else 0)) / (2 * np.pi)
        # close the loop
        closing = np.arctan2(d[0, 1], d[0, 0]) - np.arctan2(d[-1, 1], d[-1, 0])
        closing = (closing + np.pi) % (2 * np.pi) - np.pi
        wn = (ang[-1] - ang[0] + closing) / (2 * np.pi)
        rec["windings"][str(r)] = float(round(wn, 3))
    out.append(rec)
    print("C/Cc=%.5f (dx=%+.3f) trace %+.4f windings(r): %s" % (x, dx, np.trace(J), rec["windings"]), flush=True)
json.dump(out, open(os.path.join(HERE, "phase2d_winding.json"), "w"), indent=1)
