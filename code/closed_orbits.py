# -*- coding: utf-8 -*-
"""Exact time-symmetric cold cycles of the closed matter-coupled scalar model (no dissipation, no assumed bounce).

The equations of closed_cycle.py are invariant under time reversal.  A trajectory that starts at the classical neck (a_dot = 0, phi_dot = 0, phi = seed s) is mapped to its own
time reverse if it later reaches a turnaround with phi_dot = 0 as well; the contraction then retraces the expansion in reverse, restoration returns the scalar to the seed, and the
universe returns to the neck: an exact periodic orbit with a nonsingular bounce, a held epoch, a release, a turnaround and a restoration.  That is one condition, g(s) = phi_dot at the
first turnaround = 0, on the one free parameter s, so such orbits are isolated (a coherent seed, as in the earlier cyclic papers), and their stability has to be computed.

    python closed_orbits.py          (writes closed_orbits.json)
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import closed_cycle as CC

HERE = os.path.dirname(os.path.abspath(__file__))


def first_turnaround(m, s, rtol=1e-12, t_max=60.0):
    """Integrate from the neck with seed s (in units of v) to the first turnaround; return (t, state) or None."""
    a0 = CC.neck_a(m, s * m.v)
    if a0 is None:
        return None
    y0 = [a0, 0.0, s * m.v, 0.0, 0.0]
    ev = lambda t, y: y[1]
    ev.direction = -1                                    # a_dot crossing zero from above: a maximum
    ev2 = lambda t, y: y[0] - 1e-3
    ev2.terminal = True
    ev3 = lambda t, y: y[0] - 1e6
    ev3.terminal = True
    sol = solve_ivp(m.rhs, (0.0, t_max), y0, method="DOP853", rtol=rtol, atol=1e-15, events=[ev, ev2, ev3], max_step=0.02)
    if sol.t_events[0].size:
        return float(sol.t_events[0][0]), sol.y_events[0][0]
    return None


def g(m, ls):
    r = first_turnaround(m, 10.0 ** ls)
    return None if r is None else float(r[1][3])


def find_orbits(m, lo=-6.0, hi=-0.5, n=140):
    grid = np.linspace(lo, hi, n)
    vals = [g(m, x) for x in grid]
    roots = []
    for i in range(n - 1):
        if vals[i] is not None and vals[i + 1] is not None and vals[i] * vals[i + 1] < 0:
            try:
                r = brentq(lambda x: g(m, x), grid[i], grid[i + 1], xtol=1e-13)
                roots.append(float(r))
            except Exception:
                pass
    return roots, grid, vals


def period_check(m, ls, n_turns=4):
    """Integrate through several turning points from the symmetric start; report the sequence and how well the second neck matches the first."""
    o = CC.run(m, seed=10.0 ** ls * m.v, t_max=120.0, max_turns=n_turns, rtol=1e-12)
    return o


def instability(m, ls, delta=1e-6):
    """Growth of a perturbation of the seed over one full cycle: compare the second-neck state of the exact seed and of a perturbed seed."""
    res = []
    for d in (0.0, delta):
        o = CC.run(m, seed=10.0 ** ls * (1 + d) * m.v, t_max=120.0, max_turns=2, rtol=1e-13)
        res.append(o)
    return res


if __name__ == "__main__":
    out = []
    for q in (1.0,):
        for rc in (0.05, 0.2):
            for mu in (20.0,):
                for x in (0.3, 0.6):
                    m = CC.Model(q=q, rc=rc, mu=mu, C=x * CC.C_CRIT)
                    roots, grid, vals = find_orbits(m)
                    rec = dict(q=q, rc=rc, mu=mu, C_over_Ccrit=x, log10_seed_over_v_roots=roots, orbits=[])
                    print("q=%g rc=%g mu=%g C/Cc=%.2f: %d symmetric-turnaround roots in log10(s/v) in [-6,-0.5]: %s" % (q, rc, mu, x, len(roots), ["%.4f" % r for r in roots]), flush=True)
                    for r in roots[:3]:
                        o = period_check(m, r, n_turns=4)
                        seq = "".join("B" if t["kind"] == "bounce" else "T" for t in o["turns"])
                        a_turns = [round(t["a"], 4) for t in o["turns"]]
                        phis = [round(t["phi_over_v"], 6) for t in o["turns"]]
                        print("   root %.4f -> turning points %s a=%s phi/v=%s | %s" % (r, seq, a_turns, phis, o["end"]), flush=True)
                        rec["orbits"].append(dict(log10_seed=r, sequence=seq, a=a_turns, phi_over_v=phis, end=o["end"], turns=o["turns"]))
                    out.append(rec)
    json.dump(out, open(os.path.join(HERE, "closed_orbits.json"), "w"), indent=1)
