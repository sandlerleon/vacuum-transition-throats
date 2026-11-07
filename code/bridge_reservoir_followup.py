# -*- coding: utf-8 -*-
"""Follow-up on the stable-looking reservoir orbit: step convergence, exact multipliers, long-run survival, continuation in (g, m_chi), and the g = 0 control."""
import json, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import bridge_reservoir as BR
import closed_cycle as CC

HERE = os.path.dirname(os.path.abspath(__file__))
B = json.load(open(os.path.join(HERE, "closed_stability_breather.json")))
O = [b for b in B if b["mu"] == 20.0 and b["rc"] == 0.2][0]


def model(mchi, g):
    return BR.TwoField(O["q"], O["rc"], O["mu"], O["C_over_Ccrit"] * CC.C_CRIT, mchi, g)


def solve_orbit(m, guess):
    F = lambda x: (lambda r: np.array([10.0, 10.0]) if r is None else r)(BR.residual(m, x[0] * m.v, x[1]))
    x, info, ier, msg = fsolve(F, guess, full_output=True, xtol=1e-13)
    return (x, np.linalg.norm(F(x))) if ier == 1 else (None, None)


def survive(m, s, rel, max_bounces=400, seedrng=0):
    """Perturb the scalar start state by relative amount rel (random direction in the 4 scalar variables), count bounces until singular/runaway; return count and max |phi|/v deviation of the neck states."""
    rng = np.random.default_rng(seedrng)
    sc = np.array([m.v, m.v, 1.0, 1.0])
    d = rng.normal(size=4); d /= np.linalg.norm(d)
    s2 = np.array(s) + rel * d * sc
    a0 = m.neck_a(*s2)
    if a0 is None:
        return 0, None
    y0 = [a0, 0.0, s2[0], s2[1], s2[2], s2[3]]
    ev = lambda t, y: y[1]; ev.direction = 0
    small = lambda t, y: y[0] - 1e-3; small.terminal = True
    big = lambda t, y: y[0] - 1e6; big.terminal = True
    T = 0.62 * (max_bounces + 2)
    sol = solve_ivp(m.rhs, (0.0, T), y0, method="DOP853", rtol=1e-11, atol=1e-14, events=[ev, small, big], max_step=0.05)
    n = 0; dev = 0.0
    for t, y in zip(sol.t_events[0], sol.y_events[0]):
        if t < 1e-6: continue
        if m.rhs(t, y)[1] > 0:
            n += 1; dev = max(dev, float(np.linalg.norm((np.array([y[2], y[3], y[4], y[5]]) - np.array(s)) / sc)))
    ended = "singular" if sol.t_events[1].size else ("runaway" if sol.t_events[2].size else "t_max")
    return n, dev, ended, float(sol.t[-1])


if __name__ == "__main__":
    res = {}
    m = model(20.0, 1000.0)
    x, r = solve_orbit(m, [0.2318, 0.0452])
    s = (x[0] * m.v, 0.0, x[1], 0.0)
    print("orbit phi0/v=%.8f chi0=%.8f residual %.1e" % (x[0], x[1], r))
    res["orbit"] = dict(phi0_over_v=float(x[0]), chi0=float(x[1]), residual=float(r))
    for h in (1e-5, 1e-6, 1e-7, 1e-8):
        fl = BR.floquet(m, s, h=h)
        print(" h=%g multipliers %s det %.9f" % (h, ["%.6f%+.6fj" % (a, b) for a, b in fl["multipliers"]], fl["det"]))
    fl = BR.floquet(m, s, h=1e-6)
    res["floquet"] = fl
    ang = [float(np.angle(complex(a, b))) for a, b in fl["multipliers"]]
    print(" arg(mult) =", ["%.5f" % x for x in ang], " |mu|-1 =", ["%.1e" % (abs(complex(a, b)) - 1) for a, b in fl["multipliers"]])
    res["survival"] = []
    for rel in (1e-8, 1e-6, 1e-4, 1e-3, 1e-2):
        out = [survive(m, s, rel, 400, sd) for sd in range(3)]
        print(" perturbation %.0e: bounces survived %s, max deviation %s, end %s" % (rel, [o[0] for o in out], ["%.2e" % o[1] if o[1] is not None else None for o in out], [o[2] for o in out]))
        res["survival"].append(dict(rel=rel, bounces=[o[0] for o in out], dev=[o[1] for o in out], end=[o[2] for o in out]))
    # continuation in g and m_chi from the stable orbit
    res["map"] = []
    for mchi in (16.0, 18.0, 20.0, 22.0, 25.0):
        guess = [x[0], x[1] * 20.0 / mchi]
        for g in (300.0, 600.0, 800.0, 1000.0, 1200.0, 1500.0, 2000.0):
            mm = model(mchi, g)
            xs, rr = solve_orbit(mm, guess)
            if xs is None or xs[1] < 0.005:
                print(" m=%g g=%g: no orbit from guess" % (mchi, g)); continue
            guess = list(xs)
            f2 = BR.floquet(mm, (xs[0] * mm.v, 0.0, xs[1], 0.0), h=1e-6)
            if f2 is None:
                continue
            print(" m=%g g=%g: phi0/v=%.4f chi0=%.4f T=%.4f max|mu|=%.4f  |mu|=%s" % (mchi, g, xs[0], xs[1], f2["period"], f2["max_abs"], ",".join("%.4f" % q for q in f2["multipliers_abs"])), flush=True)
            res["map"].append(dict(mchi=mchi, g=g, phi0_over_v=float(xs[0]), chi0=float(xs[1]), period=f2["period"], max_abs=f2["max_abs"], multipliers_abs=f2["multipliers_abs"]))
    json.dump(res, open(os.path.join(HERE, "bridge_reservoir_followup.json"), "w"), indent=1)
