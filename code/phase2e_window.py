# -*- coding: utf-8 -*-
"""Phase 2E, item 5: recheck the elliptic window: nonlinear survival of the neck map near rational rotation numbers (resonances) and convergence.  (independent Cycle code)"""
import json, os
import numpy as np
from scipy.optimize import brentq
from phase2d_cycle_validation import Cycle

HERE = os.path.dirname(os.path.abspath(__file__))


def trace_of(x, guess):
    c = Cycle(x); p = c.orbit(guess, span=0.1)
    return float(np.trace(c.jac_neck(p, 0.0))), p / c.v


# locate C where the rotation number is exactly 2/5, 3/8, 1/3? (trace = 2 cos(2 pi r))
targets = {"2/5": 0.4, "3/8": 0.375, "0.36": 0.36}
guess_of = lambda x: 0.19931 + (0.2229 - 0.19931) / 0.005 * (x - 0.345) if x < 0.3515 else 0.2229 + (0.2499 - 0.2229) / 0.005 * (x - 0.35)
out = dict(points=[])
Cs = {}
for name, r in (("2/5", 0.4), ("3/8", 0.375)):
    tt = 2 * np.cos(2 * np.pi * r)
    f = lambda x: trace_of(x, guess_of(x))[0] - tt
    Cs[name] = brentq(f, 0.349, 0.3598, xtol=1e-9)
Cs["near_flip"] = 0.3470
Cs["centre"] = 0.3500
Cs["near_end"] = 0.3598
print("resonance locations:", {k: round(v, 6) for k, v in Cs.items()}, flush=True)
for name, x in Cs.items():
    c = Cycle(x)
    tr, g = trace_of(x, guess_of(x))
    p = g * c.v
    rot = float(np.arccos(tr / 2) / (2 * np.pi))
    rec = dict(name=name, x=float(x), phi0_over_v=float(g), trace=tr, rotation=rot, survival=[])
    for rel in (1e-3, 1e-2, 3e-2):
        row = []
        for sd in range(2):
            rng = np.random.default_rng(sd); d = rng.normal(size=2); d /= np.linalg.norm(d)
            z = np.array([p + rel * d[0] * c.v, rel * d[1] * c.v]); n = 0; dev = 0.0
            for k in range(150):
                r_ = c.P(z[0], z[1])
                if r_ is None: break
                z = np.array([r_[0], r_[1]]); n += 1; dev = max(dev, float(np.hypot(z[0] - p, z[1]) / c.v))
            row.append((n, dev))
        rec["survival"].append(dict(rel=rel, bounces=[a for a, b in row], dev=[b for a, b in row]))
    out["points"].append(rec)
    print("%-9s C/Cc=%.6f phi0/v=%.5f trace=%+.5f rot=%.5f | survival(150 cap): %s" % (name, x, g, tr, rot, [(s["rel"], s["bounces"], ["%.1e" % d for d in s["dev"]]) for s in rec["survival"]]), flush=True)
json.dump(out, open(os.path.join(HERE, "phase2e_window.json"), "w"), indent=1)
