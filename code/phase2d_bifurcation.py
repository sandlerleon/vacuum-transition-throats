# -*- coding: utf-8 -*-
"""Phase 2D step 2: classify the bifurcation of the main orbit branch near C/Cc = 0.3458 and find the end of the elliptic window (independent code, phase2d_cycle_validation.Cycle)."""
import json, os
import numpy as np
from scipy.optimize import fsolve, brentq
from phase2d_cycle_validation import Cycle

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}

# ---- A. extend the branch beyond 0.36 with secant-predicted guesses, reporting where it stops
print("A. branch beyond 0.36")
pts = [(0.3575, 0.26508), (0.36, 0.28182)]
out["branch_ext"] = []
x = 0.36
while x < 0.40:
    x = round(x + 0.0025, 5)
    g = pts[-1][1] + (pts[-1][1] - pts[-2][1]) / (pts[-1][0] - pts[-2][0]) * (x - pts[-1][0])
    c = Cycle(x)
    p = None
    for span in (0.04, 0.1, 0.25):
        p = c.orbit(g, span=span)
        if p is not None:
            break
    if p is None:
        print("  x=%.4f: orbit not found near guess %.4f" % (x, g)); out["branch_ext"].append(dict(x=x, found=False, guess=float(g))); break
    J = c.jac_neck(p, 0.0)
    tr = float(np.trace(J))
    pts.append((x, p / c.v))
    out["branch_ext"].append(dict(x=x, found=True, phi0_over_v=float(p / c.v), trace=tr, period=float(c.flow_to_bounce(p, 0.0)[0]), det=float(np.linalg.det(J))))
    print("  x=%.4f phi0/v=%.5f T=%.5f trace=%+.5f det=%.6f rot=%s" % (x, p / c.v, out["branch_ext"][-1]["period"], tr, np.linalg.det(J), "%.4f" % (np.arccos(tr / 2) / 2 / np.pi) if abs(tr) < 2 else "-"), flush=True)

# ---- B. flip: displacement of P^2 along the eigenvector at C/Cc = x* + dx
print("B. P^2 displacement along the -1 eigenvector")
xs = 0.34584411
out["flip_scan"] = []
for dx in (-0.006, -0.003, 0.003, 0.006, 0.012):
    x = xs + dx
    c = Cycle(x)
    g = 0.19931 + (0.22295 - 0.19931) / 0.005 * (x - 0.345)
    p = c.orbit(g, span=0.2)
    J = c.jac_neck(p, 0.0)
    w, V = np.linalg.eig(J)
    k = int(np.argmin(np.abs(w + 1)))
    e = np.real(V[:, k]); e /= np.linalg.norm(e)
    row = []
    for s in (-0.3, -0.2, -0.1, -0.05, -0.02, -0.01, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3):
        x0 = np.array([p, 0.0]) + s * e * c.v
        r1 = c.P(*x0)
        r2 = None if r1 is None else c.P(r1[0], r1[1])
        if r2 is None:
            row.append((s, None)); continue
        d = np.array([r2[0], r2[1]]) - x0
        row.append((s, float(d @ e / c.v), float(np.linalg.norm(d) / c.v)))
    out["flip_scan"].append(dict(dx=dx, x=float(x), trace=float(np.trace(J)), eig=[complex(z).real for z in w], rows=row))
    print(" C/Cc=%.5f trace %+.4f mult %s" % (x, np.trace(J), np.round(w, 4)))
    print("    s:    " + " ".join("%7.2f" % r[0] for r in row))
    print("    dP2.e:" + " ".join("%7.4f" % r[1] if r[1] is not None else "   None" for r in row), flush=True)

# ---- C. 2-D multi-start search for fixed points of P^2 that are not fixed points of P, on both sides
print("C. multi-start period-2 search")
out["period2"] = []
rng = np.random.default_rng(1)
for dx in (-0.003, 0.003, 0.006, 0.012):
    x = xs + dx
    c = Cycle(x)
    g = 0.19931 + (0.22295 - 0.19931) / 0.005 * (x - 0.345)
    p = c.orbit(g, span=0.2)

    def P2(z):
        r1 = c.P(p + z[0] * c.v, z[1] * c.v)
        if r1 is None: return np.array([9.0, 9.0])
        r2 = c.P(r1[0], r1[1])
        if r2 is None: return np.array([9.0, 9.0])
        return np.array([(r2[0] - p) / c.v - z[0], r2[1] / c.v - z[1]])
    found = []
    for k in range(36):
        amp = 10 ** rng.uniform(-2.3, -0.5); th = rng.uniform(0, 2 * np.pi)
        z0 = amp * np.array([np.cos(th), np.sin(th)])
        try:
            z, info, ier, msg = fsolve(P2, z0, full_output=True, xtol=1e-12)
        except Exception:
            continue
        if ier == 1 and np.linalg.norm(P2(z)) < 1e-9 and np.linalg.norm(z) > 3e-4 and not any(np.linalg.norm(z - f) < 1e-5 for f in found):
            found.append(z)
    out["period2"].append(dict(dx=dx, x=float(x), points=[[float(a) for a in f] for f in found]))
    print(" C/Cc=%.5f: %d distinct non-trivial fixed points of P^2: %s" % (x, len(found), [np.round(f, 4).tolist() for f in found]), flush=True)
json.dump(out, open(os.path.join(HERE, "phase2d_bifurcation.json"), "w"), indent=1)
