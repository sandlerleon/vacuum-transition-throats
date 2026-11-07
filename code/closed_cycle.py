# -*- coding: utf-8 -*-
"""Phase 1 of the pre-existing-bridge programme: does the matter-coupled scalar model complete cold cycles in a CLOSED universe with no assumed bounce?

Units: 8 pi G = 1, R_dS = 1 (false-vacuum density V_f = 3), curvature K = +1 with a the curvature radius, c = 1.  Dust parameter C = rho_* a^3 / 3 (conserved, rho_* the coupled
dark-matter 'bare' density).  Scalar: V_0 = V_c - mu^2 phi^2/2 + lam phi^4/4 with V_c = V_f, coupling A = 1 + phi^2/(2 M^2), rho_crit = mu^2 M^2, true vacuum V_0(v) = -q V_f.

Exact equations (no bounce prescription):
    a_ddot = -(a/6)(rho + 3P),  rho = phi_dot^2/2 + V_0 + A rho_*,  P = phi_dot^2/2 - V_0,     phi_ddot + 3 (a_dot/a) phi_dot + V_0' + rho_* A' = 0,
with the Friedmann constraint  a_dot^2 = a^2 rho / 3 - 1  monitored.  A bounce is a minimum of a that the equations produce (a_dot = 0, a_ddot > 0), a turnaround a maximum.

    python closed_cycle.py            (scan -> closed_cycle.json)
"""
import json
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
VF = 3.0
C_CRIT = 2.0 / (3.0 * np.sqrt(3.0))


class Model:
    def __init__(self, q=1.0, rc=0.05, mu=30.0, C=0.2, ups=0.0, ups_fn=None):
        self.q, self.rc, self.mu, self.C, self.ups, self.ups_fn = q, rc, mu, C, ups, ups_fn
        self.rho_crit = rc * VF
        self.M2 = self.rho_crit / mu ** 2
        self.v2 = 4.0 * (1.0 + q) * VF / mu ** 2
        self.lam = mu ** 2 / self.v2
        self.v = np.sqrt(self.v2)

    def V0(self, p):
        return VF - 0.5 * self.mu ** 2 * p * p + 0.25 * self.lam * p ** 4

    def rho_star(self, a):
        return 3.0 * self.C / a ** 3

    def rho(self, a, p, pd, rr=0.0):
        return 0.5 * pd * pd + self.V0(p) + (1.0 + p * p / (2.0 * self.M2)) * self.rho_star(a) + rr

    def rhs(self, t, y):
        a, ad, p, pd, rr = y
        rs = self.rho_star(a)
        A = 1.0 + p * p / (2.0 * self.M2)
        add = -(a / 6.0) * (2.0 * pd * pd - 2.0 * self.V0(p) + A * rs + 2.0 * rr)
        u = self.ups if self.ups_fn is None else self.ups_fn(p, pd, rr)
        pdd = -(3.0 * ad / a + u) * pd - (-self.mu ** 2 * p + self.lam * p ** 3 + rs * p / self.M2)
        rrd = -4.0 * ad / a * rr + u * pd * pd
        return [ad, add, pd, pdd, rrd]

    def constraint(self, y):
        a, ad, p, pd, rr = y
        return ad * ad - (a * a * self.rho(a, p, pd, rr) / 3.0 - 1.0)


def neck_a(m, p, pd=0.0, S=0.0):
    """Radius a at which a_dot = 0 satisfies the Friedmann constraint exactly for the scalar state (p, pd): root of a^2 rho/3 = 1 (outer branch)."""
    f = lambda a: a * a * m.rho(a, p, pd, S / a ** 4) / 3.0 - 1.0
    lo, hi = (m.C / 2.0) ** (1.0 / 3.0) * 1.0001, 10.0
    if f(lo) * f(hi) > 0:
        return None
    return brentq(f, lo, hi)


def bounce_start(C):
    """Outer root a2 of a^2 + C/a - 1 = 0 (R_dS = 1)."""
    return brentq(lambda a: a * a + C / a - 1.0, (C / 2.0) ** (1.0 / 3.0), 10.0)


def run(m, seed=1e-12, t_max=400.0, max_turns=6, a_min=1e-3, a_max=1e7, rtol=1e-11, y0=None, S0=0.0):
    """Start at the classical neck (a_dot = 0, scalar at the held point with a small seed) and integrate; list the turning points of a."""
    a0 = neck_a(m, seed, 0.0, S0)
    if a0 is None:
        return dict(params=dict(q=m.q, rc=m.rc, mu=m.mu, C=m.C, C_over_Ccrit=m.C / C_CRIT, seed=seed, ups=m.ups), turns=[], end="no classical neck for this seed", a_neck_initial=None, constraint_initial=None)
    y = [a0, 0.0, seed, 0.0, S0 / a0 ** 4] if y0 is None else list(y0)
    # constraint correction: a_dot = 0 requires a^2 rho/3 = 1; with phi = seed the residual is O(seed^2); report it
    res0 = m.constraint(y)
    turns = []

    def ev_turn(t, y):
        return y[1]
    ev_turn.direction = 0

    def ev_small(t, y):
        return y[0] - a_min
    ev_small.terminal = True

    def ev_big(t, y):
        return y[0] - a_max
    ev_big.terminal = True

    out = dict(params=dict(q=m.q, rc=m.rc, mu=m.mu, C=m.C, C_over_Ccrit=m.C / C_CRIT, seed=seed, ups=m.ups), a_neck_initial=a0, constraint_initial=float(res0), turns=[], end=None)
    sol = solve_ivp(m.rhs, (0.0, t_max), y, method="DOP853", rtol=rtol, atol=1e-14, events=[ev_turn, ev_small, ev_big], max_step=0.05)
    for k, tt in enumerate(sol.t_events[0]):
        if tt < 1e-6:
            continue                                       # the starting neck itself
        ye = sol.y_events[0][k]
        if len(out["turns"]) >= max_turns:
            break
        add = m.rhs(tt, ye)[1]
        rho = m.rho(ye[0], ye[2], ye[3], ye[4])
        vac = m.V0(ye[2])
        rec = dict(t=float(tt), a=float(ye[0]), kind="bounce" if add > 0 else "turnaround", phi_over_v=float(ye[2] / m.v), kinetic=float(0.5 * ye[3] ** 2),
                   rho_over_Vf=float(rho / VF), V0_over_Vf=float(vac / VF), rho_dm=float(m.rho_star(ye[0]) * (1 + ye[2] ** 2 / (2 * m.M2))), constraint=float(m.constraint(ye)))
        rec["rho_plus_3P_over_Vf"] = float((2 * ye[3] ** 2 - 2 * vac + rec["rho_dm"] + 2 * ye[4]) / VF)
        rec["radiation_over_Vf"] = float(ye[4] / VF)
        out["turns"].append(rec)
    if sol.t_events[1].size:
        out["end"] = "singular (a < %g) at t = %.3f" % (a_min, sol.t_events[1][0])
    elif sol.t_events[2].size:
        out["end"] = "a > %g at t = %.3f (expands without turning)" % (a_max, sol.t_events[2][0])
    elif len(out["turns"]) >= max_turns:
        out["end"] = "max turning points reached"
    else:
        out["end"] = "t_max reached"
    return out


def summarize(o):
    p = o["params"]
    seq = "".join("B" if t["kind"] == "bounce" else "T" for t in o["turns"])
    return "C/Cc=%.2f rc=%.3f mu=%g q=%g ups=%g seed=%.0e | %s | %s" % (p["C_over_Ccrit"], p["rc"], p["mu"], p["q"], p["ups"], p["seed"], seq, o["end"])


if __name__ == "__main__":
    res = []
    for ups in (0.0, 0.3, 3.0):
        for rc in (0.02, 0.05, 0.2):
            for mu in (20.0, 30.0):
                for x in (0.1, 0.3, 0.6, 0.9):
                    for seed in (1e-3, 1e-8, 1e-20):
                        m = Model(q=1.0, rc=rc, mu=mu, C=x * C_CRIT, ups=ups)
                        o = run(m, seed=seed * m.v)
                        res.append(o)
                        print(summarize(o), flush=True)
    json.dump(res, open(os.path.join(HERE, "closed_cycle.json"), "w"), indent=1)
    n2 = sum(1 for o in res if sum(1 for t in o["turns"] if t["kind"] == "bounce") >= 1)
    print("runs with at least one further bounce after the starting neck:", n2, "of", len(res))
