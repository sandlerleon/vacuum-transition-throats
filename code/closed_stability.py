# -*- coding: utf-8 -*-
"""Stability of the exact time-symmetric cold cycles of closed_orbits.py.

Poincare map on the neck: P(phi, phi_dot) = (phi, phi_dot) at the next bounce (a_dot = 0, a_ddot > 0), with a fixed by the Friedmann constraint a^2 rho = 3.
A periodic orbit is a fixed point; its Floquet multipliers are the eigenvalues of the Jacobian of P (central differences).  |lambda| > 1: unstable (hyperbolic);
|lambda| = 1 (complex pair): linearly stable (elliptic).  The map is also iterated from perturbed starts to see how many bounces survive.

    python closed_stability.py          (writes closed_stability.json)
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import closed_cycle as CC
import closed_orbits as CO

HERE = os.path.dirname(os.path.abspath(__file__))


def neck_radius(m, p, pd):
    """Smallest-branch-free choice: the root a of a^2 rho(a; p, pd) = 3 closest to the vacuum-only neck."""
    return CC.neck_a(m, p, pd)


def neck_map(m, p, pd, t_max=80.0, rtol=1e-13):
    a0 = neck_radius(m, p, pd)
    if a0 is None:
        return None
    y0 = [a0, 0.0, p, pd, 0.0]
    ev = lambda t, y: y[1]
    ev.direction = 0
    ev2 = lambda t, y: y[0] - 1e-3
    ev2.terminal = True
    sol = solve_ivp(m.rhs, (0.0, t_max), y0, method="DOP853", rtol=rtol, atol=1e-15, events=[ev, ev2], max_step=0.02)
    for k, tt in enumerate(sol.t_events[0]):
        if tt < 1e-6:
            continue
        ye = sol.y_events[0][k]
        if m.rhs(tt, ye)[1] > 0:
            return float(ye[2]), float(ye[3]), float(tt), float(ye[0])
    return None


def jacobian(m, p, pd, h=1e-7):
    ps = m.v
    cols = []
    for dp, dpd in ((h * ps, 0.0), (0.0, h * ps)):
        a = neck_map(m, p + dp, pd + dpd); b = neck_map(m, p - dp, pd - dpd)
        if a is None or b is None:
            return None
        cols.append([(a[0] - b[0]) / (2 * h * ps), (a[1] - b[1]) / (2 * h * ps)])
    return np.array(cols).T


def analyse(q, rc, mu, x, ls):
    m = CC.Model(q=q, rc=rc, mu=mu, C=x * CC.C_CRIT)
    s = 10.0 ** ls * m.v
    r = neck_map(m, s, 0.0)
    rec = dict(q=q, rc=rc, mu=mu, C_over_Ccrit=x, log10_seed_over_v=ls, period=None, a_neck=None, a_turn=None, fixed_point_error=None, multipliers=None, determinant=None, trace=None, survive=None)
    if r is None:
        return rec
    rec["period"] = r[2]; rec["a_neck"] = r[3]
    rec["fixed_point_error"] = [(r[0] - s) / m.v, r[1] / m.v]
    J = jacobian(m, s, 0.0)
    if J is not None:
        ev = np.linalg.eigvals(J)
        rec["multipliers"] = [[float(z.real), float(z.imag)] for z in ev]
        rec["determinant"] = float(np.linalg.det(J)); rec["trace"] = float(np.trace(J))
    # survival: iterate the map from the exact seed and from perturbed seeds
    surv = {}
    for d in (0.0, 1e-12, 1e-9, 1e-6):
        p, pd = s * (1 + d), 0.0
        n = 0
        for _ in range(30):
            nxt = neck_map(m, p, pd)
            if nxt is None or abs(nxt[0] - s) > 0.5 * m.v or abs(nxt[1]) > 5 * m.v:
                break
            p, pd = nxt[0], nxt[1]; n += 1
        surv[str(d)] = n
    rec["bounces_survived_by_relative_seed_error"] = surv
    return rec


if __name__ == "__main__":
    scan = json.load(open(os.path.join(HERE, "closed_orbit_scan.json")))
    cases = []
    for row in scan:
        if row["closing"]:
            best = max(row["closing"], key=lambda g: g["excursion"])
            if best["excursion"] > 3.0:
                cases.append((row["q"], row["rc"], row["mu"], row["C_over_Ccrit"], best["log10_seed"], best["excursion"]))
    out = []
    for q, rc, mu, x, ls, exc in cases:
        m = CC.Model(q=q, rc=rc, mu=mu, C=x * CC.C_CRIT)
        try:
            ls = brentq(lambda z: CO.g(m, z), ls - 0.02, ls + 0.02, xtol=1e-14)
        except Exception:
            pass
        r = analyse(q, rc, mu, x, ls)
        r["excursion"] = exc
        out.append(r)
        mult = r["multipliers"]
        print("mu=%g rc=%g C/Cc=%.1f seed=%.4f exc=%.1f | period %s | fp err %s | |multipliers| %s | det %s | survived %s" %
              (mu, rc, x, ls, exc, None if r["period"] is None else round(r["period"], 3), ["%.0e" % e for e in (r["fixed_point_error"] or [])],
               None if not mult else ["%.3g" % np.hypot(*z) for z in mult], None if r["determinant"] is None else "%.3g" % r["determinant"], r.get("bounces_survived_by_relative_seed_error")), flush=True)
    json.dump(out, open(os.path.join(HERE, "closed_stability.json"), "w"), indent=1)
