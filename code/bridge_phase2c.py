# -*- coding: utf-8 -*-
"""Phase 2C driver: (i) the stable reservoir orbit, (ii) junction consistency with a MOVING shell around it (total rho, P), (iii) controls: g = 0, coupling switched off, no reservoir,
(iv) fine continuation of the stability window in m_chi at g = 1000.   python bridge_phase2c.py  (writes bridge_phase2c.json)"""
import json, os
import numpy as np
from scipy.integrate import solve_ivp
import bridge_reservoir as BR
import bridge_reservoir_followup as BF
import bridge_moving_shell as MS
import closed_cycle as CC

HERE = os.path.dirname(os.path.abspath(__file__))


def count_bounces(m, y0, t_max):
    ev = lambda t, y: y[1]; ev.direction = 0
    small = lambda t, y: y[0] - 1e-3; small.terminal = True
    big = lambda t, y: y[0] - 1e6; big.terminal = True
    sol = solve_ivp(m.rhs, (0.0, t_max), y0, method="DOP853", rtol=1e-11, atol=1e-14, events=[ev, small, big], max_step=0.05)
    n = sum(1 for t, y in zip(sol.t_events[0], sol.y_events[0]) if t > 1e-6 and m.rhs(t, y)[1] > 0)
    return n, ("singular" if sol.t_events[1].size else ("runaway" if sol.t_events[2].size else "t_max")), float(sol.t[-1])


if __name__ == "__main__":
    out = {}
    m = BF.model(20.0, 1000.0)
    x, r = BF.solve_orbit(m, [0.2318, 0.0452])
    s = (x[0] * m.v, 0.0, x[1], 0.0)
    fl = BR.floquet(m, s, h=1e-6)
    fr = BR.reservoir_fraction(m, s)
    out["orbit"] = dict(mchi=20.0, g=1000.0, phi0_over_v=float(x[0]), chi0=float(x[1]), period=fl["period"], multipliers=fl["multipliers"], max_abs=fl["max_abs"], reservoir_frac=[fr[0], fr[1]], constraint_max=fr[2], a_range=[fr[3], fr[4]])
    print("stable orbit: T=%.5f max|mu|-1=%.1e reservoir energy fraction %.4f-%.4f a in [%.4f,%.4f] constraint %.1e" % (fl["period"], fl["max_abs"] - 1, fr[0], fr[1], fr[3], fr[4], fr[2]))
    # (ii) moving shell around the stable orbit
    a0 = m.neck_a(*s)
    yi0 = [a0, 0.0, s[0], 0.0, s[2], 0.0]
    out["shells"] = []
    for q in (1.0, 0.0):
        for m_frac in (0.2, 0.5, 0.8):
            for w in (-1.0, -0.5, 0.0, 1.0 / 3.0, 1.0):
                for rr in MS.periodic_shell_generic(m, yi0, fl["period"], (0.52, 0.9), m_frac, q, w):
                    rr.update(q=q, m_frac=m_frac, w=w)
                    out["shells"].append(rr)
                    d = rr["diag"]
                    print("shell q=%g m/mmax=%.1f w=%+.2f chi0=%.4fpi | |mult|=%s | sigma>0 %.2f NEC %.2f | flux-identity %.1e | R[%.3f,%.3f]" % (q, m_frac, w, rr["chi0_over_pi"],
                          ",".join("%.4f" % z for z in rr["multipliers_abs"]), d["sigma_pos"], d["nec"], d["flux_identity_resid"], d["R_min"], d["R_max"]), flush=True)
    # (iii) controls
    ctrl = {}
    m0 = BF.model(20.0, 0.0)
    x0, r0 = BF.solve_orbit(m0, [0.07, x[1]])
    if x0 is not None:
        fl0 = BR.floquet(m0, (x0[0] * m0.v, 0.0, x0[1], 0.0), h=1e-6)
        ctrl["g0_orbit"] = dict(phi0_over_v=float(x0[0]), chi0=float(x0[1]), max_abs=fl0["max_abs"], multipliers_abs=fl0["multipliers_abs"])
        print("control g=0 (decoupled reservoir, same chi0): phi0/v=%.4f max|mu|=%.4g |mu|=%s" % (x0[0], fl0["max_abs"], fl0["multipliers_abs"]))
    ctrl["coupling_off"] = {}
    for mchi, g, label in ((20.0, 0.0, "g=0 from the stable orbit's initial data"), (20.0, 1000.0, "g=1000 (stable orbit)")):
        mm = BF.model(mchi, g)
        a0 = mm.neck_a(*s)
        n, end, tend = count_bounces(mm, [a0, 0.0, s[0], 0.0, s[2], 0.0], 0.62 * 402)
        print("  %-42s: %d bounces, end %s at t=%.1f" % (label, n, end, tend))
        ctrl["coupling_off"][label] = dict(bounces=n, end=end, t_end=tend)
    ctrl["no_reservoir"] = dict(phi0_over_v=10 ** BF.O["log10_seed_over_v"], max_abs=12.667543)
    out["controls"] = ctrl
    # (iv) fine continuation in m_chi at g = 1000
    out["m_window"] = []
    for direction in (-1, 1):
        guess = [x[0], x[1]]
        for k in range(1, 41):
            mchi = 20.0 + direction * 0.25 * k
            mm = BF.model(mchi, 1000.0)
            xs, rr = BF.solve_orbit(mm, guess)
            if xs is None or xs[1] < 0.003:
                print(" m=%.2f: orbit lost" % mchi); break
            guess = list(xs)
            f2 = BR.floquet(mm, (xs[0] * mm.v, 0.0, xs[1], 0.0), h=1e-6)
            if f2 is None:
                break
            out["m_window"].append(dict(mchi=mchi, phi0_over_v=float(xs[0]), chi0=float(xs[1]), period=f2["period"], max_abs=f2["max_abs"]))
            print(" m=%.2f phi0/v=%.4f chi0=%.4f T=%.4f max|mu|=%.6f" % (mchi, xs[0], xs[1], f2["period"], f2["max_abs"]), flush=True)
            if f2["max_abs"] > 3.0:
                break
    json.dump(out, open(os.path.join(HERE, "bridge_phase2c.json"), "w"), indent=1)
