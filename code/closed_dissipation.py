# -*- coding: utf-8 -*-
"""Phase 1B: can self-consistent dissipation stabilize the exact cold cycles of the closed matter-coupled model?

(a) Explicit periodicity check of the exact (non-dissipative) orbits: the full state (a, a_dot, phi, phi_dot, rho_r) after one period against the start.
(b) Dissipation with the energy going into radiation, never deleted:  phi_ddot + (3H + G) phi_dot + V' = 0,  rho_r_dot + 4 H rho_r = G phi_dot^2,  with rho_r in the Friedmann and
    acceleration equations.  G is (1) constant, (2) gated to rapid field motion (|phi_dot| above half the peak of the undamped cycle), (3) thermal, G = C_T rho_r^(3/4)
    (Paper 1's Upsilon = C_T T^3 form) from a small radiation seed.
For each case: bounces completed from the exact orbit's start, radiation accumulated at the necks, entropy parameter S = rho_r a^4 per bounce, the stretching of a perturbation over the
first return (the Floquet proxy when G = 0), constraint residuals, and the smallest scale factor.

    python closed_dissipation.py          (writes closed_dissipation.json)
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp

import closed_cycle as CC

HERE = os.path.dirname(os.path.abspath(__file__))
VF = CC.VF


def necks(m, y0, t_max, rtol=1e-12):
    """Integrate from y0 and return the list of bounce states (a_dot = 0, a_ddot > 0), the turnarounds, the end status, and the max |constraint| at the events."""
    ev = lambda t, y: y[1]
    ev.direction = 0
    small = lambda t, y: y[0] - 1e-3
    small.terminal = True
    big = lambda t, y: y[0] - 1e6
    big.terminal = True
    sol = solve_ivp(m.rhs, (0.0, t_max), y0, method="DOP853", rtol=rtol, atol=1e-15, events=[ev, small, big], max_step=0.02)
    B, T, res = [], [], []
    for k, tt in enumerate(sol.t_events[0]):
        if tt < 1e-6:
            continue
        ye = sol.y_events[0][k]
        res.append(abs(m.constraint(ye)))
        (B if m.rhs(tt, ye)[1] > 0 else T).append((float(tt), ye))
    end = "singular" if sol.t_events[1].size else ("runaway" if sol.t_events[2].size else "t_max")
    return B, T, end, (max(res) if res else None), sol


def start_state(m, seed, S0=0.0):
    a0 = CC.neck_a(m, seed, 0.0, S0)
    if a0 is None:
        return None
    return [a0, 0.0, seed, 0.0, S0 / a0 ** 4]


def case_runs(m, y0, t_max):
    B, T, end, res, sol = necks(m, y0, t_max)
    rad = [float(b[1][4] / VF) for b in B]
    S = [float(b[1][4] * b[1][0] ** 4) for b in B]
    return dict(bounces=len(B), turnarounds=len(T), end=end, a_neck=[float(b[1][0]) for b in B], radiation_over_Vf=rad, S=S, max_constraint=res,
                a_min=float(min([y0[0]] + [b[1][0] for b in B])), first_return=None if not B else [float(B[0][1][2]), float(B[0][1][3])]), sol


def stretch(m, seed, S0, t_max, d=1e-8):
    """Stretching of a seed perturbation over the first return: |P(x+d) - P(x)| / (d * seed) in (phi, phi_dot)."""
    out = []
    for dd in (0.0, d):
        y0 = start_state(m, seed * (1 + dd), S0)
        B, T, end, res, sol = necks(m, y0, t_max)
        if not B:
            return None
        out.append(np.array([B[0][1][2], B[0][1][3]]))
    return float(np.linalg.norm(out[1] - out[0]) / (d * seed))


if __name__ == "__main__":
    B1 = json.load(open(os.path.join(HERE, "closed_stability_breather.json")))
    W = [w for w in json.load(open(os.path.join(HERE, "closed_stability.json"))) if w["mu"] == 10.0 and w["rc"] == 0.05 and abs(w["C_over_Ccrit"] - 0.9) < 1e-9]
    orbits = [("breather", b) for b in B1] + [("wide", w) for w in W]
    out = dict(periodicity=[], dissipation=[])
    for kind, o in orbits:
        m0 = CC.Model(q=o["q"], rc=o["rc"], mu=o["mu"], C=o["C_over_Ccrit"] * CC.C_CRIT)
        seed = 10 ** o["log10_seed_over_v"] * m0.v
        y0 = start_state(m0, seed)
        period = o["period"]
        # (a) periodicity of the full state after one period
        B, T, end, res, sol = necks(m0, y0, t_max=period * 1.5)
        yR = B[0][1]
        per = dict(kind=kind, mu=o["mu"], rc=o["rc"], C_over_Ccrit=o["C_over_Ccrit"], period=period, start=[float(v) for v in y0],
                   after_one_period=[float(v) for v in yR], delta=[float(yR[i] - y0[i]) for i in range(5)], delta_over_v=[float((yR[2] - y0[2]) / m0.v), float((yR[3] - y0[3]) / m0.v)],
                   constraint_start=float(m0.constraint(y0)), constraint_end=float(m0.constraint(yR)))
        out["periodicity"].append(per)
        print("%s mu=%g rc=%g C/Cc=%.1f: |a-a0|=%.1e |a_dot|=%.1e |dphi|/v=%.1e |dphi_dot|/v=%.1e constraint %.1e / %.1e" %
              (kind, o["mu"], o["rc"], o["C_over_Ccrit"], abs(per["delta"][0]), abs(per["delta"][1]), abs(per["delta_over_v"][0]), abs(per["delta_over_v"][1]), per["constraint_start"], per["constraint_end"]), flush=True)
        t_max = 40.0 if kind == "wide" else 60.0
        base, solb = case_runs(m0, y0, t_max)
        c_gate = 0.5 * float(np.max(np.abs(solb.y[3][: max(5, int(len(solb.t) * min(1.0, 1.5 * period / solb.t[-1])))])))
        cases = [("none", 0.0, None, 0.0)]
        cases += [("constant", g, None, 0.0) for g in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 5.0)]
        cases += [("gated", g, (lambda g=g, c=c_gate: (lambda p, pd, rr: g * 0.5 * (1.0 + np.tanh((abs(pd) - c) / (0.2 * c)))))(), 0.0) for g in (1e-2, 1e-1, 1.0, 5.0)]
        cases += [("thermal", ct, (lambda ct=ct: (lambda p, pd, rr: ct * max(rr, 0.0) ** 0.75))(), 1e-6) for ct in (1e-1, 1.0, 10.0, 100.0)]
        for name, g, fn, rad_seed in cases:
            m = CC.Model(q=o["q"], rc=o["rc"], mu=o["mu"], C=o["C_over_Ccrit"] * CC.C_CRIT, ups=g if fn is None else 0.0, ups_fn=fn)
            S0 = rad_seed * VF * 0.5 ** 4          # radiation seed: rho_r a^4 = rad_seed * V_f * a0^4 with a0 ~ 0.5-1 (order of magnitude)
            y0d = start_state(m, seed, S0)
            r, sol = case_runs(m, y0d, t_max)
            r["stretch"] = stretch(m, seed, S0, t_max) if r["bounces"] >= 1 else None
            r.update(kind=kind, mu=o["mu"], rc=o["rc"], C_over_Ccrit=o["C_over_Ccrit"], mode=name, gamma=g)
            out["dissipation"].append(r)
            print("   %-8s G=%-7g bounces %2d | end %-8s | radiation at last neck %s | stretch %s" % (name, g, r["bounces"], r["end"], None if not r["radiation_over_Vf"] else "%.2e V_f" % r["radiation_over_Vf"][-1], None if r["stretch"] is None else "%.3g" % r["stretch"]), flush=True)
    json.dump(out, open(os.path.join(HERE, "closed_dissipation.json"), "w"), indent=1)
