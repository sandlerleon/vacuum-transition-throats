# -*- coding: utf-8 -*-
"""Numerical experiment: collapse of a false-vacuum core bounded by an ingoing wall, with regular centre (dnull.py).

Units: 8 pi G = 1, false-vacuum density rho_in = 3 (so the core's de Sitter radius is R_dS = sqrt(3/rho_in) = 1), true-vacuum density rho_out = -q rho_in.
Initial data on the outgoing cone u = 0 from the centre (r = v/2): psi(v) = false vacuum for v < V_w, true vacuum beyond, with a tanh transition of width w.
The field amplitude A sets the physical wall: thickness ~ A sqrt(2/kappa), tension ~ A sqrt(kappa/2) 4/3, gravitational parameter 3 sigma^2 / (2 Delta) ~ A^2.

For every run: the exterior mass parameter m_ext = m_MS - rho_out r^3 / 3 on the initial cone beyond the wall, R_M = (3 m_ext / Delta)^(1/3) (radius at which the vacuum energy
deficit equals m_ext), the first trapped point (r_v < 0, r_u < 0), the apparent-horizon radius against the Schwarzschild-AdS value r_h (r_h^3 + r_h = m_ext for q = 1),
and the existence of anti-trapped (expanding) points inside the wall.

    python collapse_scan.py          (writes collapse_scan.json)
"""
import json
import os
import sys
import time

import numpy as np

import dnull as D

HERE = os.path.dirname(os.path.abspath(__file__))


def m_ext_on_cone(o, q=1.0, rho_in=3.0):
    """Exterior (Schwarzschild-AdS) mass parameter from the Misner-Sharp mass on the initial cone beyond the wall."""
    R, S, Ru, h = o["R"][0], o["S"][0], o["Ru"][0], o["h"]
    Vw, w = o["V_w"], o["w"]
    rho_out = -q * rho_in
    vals = []
    for vv in np.linspace(Vw + 6 * w, Vw + 9 * w, 7):
        i = int(round(vv / h))
        if i + 2 >= len(R):
            continue
        rv = (R[i + 1] - R[i - 1]) / (2 * h)
        f = np.exp(S[i])
        m = R[i] * (1 + 4 * Ru[i] * rv / f)
        vals.append(m - rho_out * R[i] ** 3 / 3.0)
    return float(np.mean(vals)) if vals else float("nan"), float(np.std(vals)) if vals else float("nan")


def horizon_root(m, q):
    """Black-hole horizon radius (R_dS = 1) of the exterior f(r) = 1 - m/r + q r^2: smallest positive real root of q r^3 + r - m = 0 (r = m for q = 0)."""
    if m <= 0:
        return None
    if q == 0:
        return float(m)
    rts = [x.real for x in np.roots([q, 0.0, 1.0, -m]) if abs(x.imag) < 1e-12 and x.real > 0]
    return float(min(rts)) if rts else None


def analyse(A, Vw, w, h, kappa=24.0, q=1.0, rho_in=3.0, smax=12.0):
    o = D.run(A=A, V_w=Vw, w=w, h=h, kappa=kappa, q=q, rho_in=rho_in, snap_every=1, smax=smax, stop_after_trap=40)
    R, Ru, S, N = o["R"], o["Ru"], o["S"], o["N"]
    pot = o["pot"]
    mx, mx_sd = m_ext_on_cone(o, q, rho_in)
    Delta = (1.0 + q) * rho_in
    RM = (3.0 * mx / Delta) ** (1.0 / 3.0) if mx > 0 else 0.0
    first, last_ok = None, 0
    at_inside = 0; at_last = 0; at_later = 0
    ended = "domain"
    min_rv_last = None; smax_last = None
    iw = int(round((Vw - 4 * w) / o["h"]))                                      # points well inside the wall (false-vacuum core)
    for k in range(R.shape[0]):
        j = k
        r = R[k]
        if not np.isfinite(r[j + 1:]).all():
            ended = "non-finite values"
            break
        if np.nanmax(S[k, j + 1:]) > 12.0:
            ended = "sigma > 12"
            break
        last_ok = k
        rv = np.zeros_like(r); rv[1:-1] = (r[2:] - r[:-2]) / (2 * o["h"])
        at_row = 0
        for i in range(j + 4, N - 8):
            if rv[i] < 0 and Ru[k, i] < 0 and first is None:
                first = dict(u=k * o["h"], v=i * o["h"], r=float(r[i]))
            if i < iw and rv[i] > 0 and Ru[k, i] > 0:
                at_inside += 1; at_row += 1
        if k >= 1:
            at_later += at_row
        at_last = at_row
        min_rv_last = float(rv[j + 4:N - 8].min()); smax_last = float(np.nanmax(S[k, j + 1:]))
    if first is not None and ended == "domain":
        ended = "stopped after the trapped surface formed"
    # anti-trapped points present in the core at the initial cone
    at0 = int(sum(1 for i in range(5, iw) if Ru[0, i] > 0))
    return dict(A=A, V_w=Vw, w=w, h=h, kappa=kappa, q=q, R_w=Vw / 2.0, m_ext=mx, m_ext_sd=mx_sd, R_M=RM, r_h_SAdS=horizon_root(mx, q), first_trapped=first,
                at_points_in_core=int(at_inside), at_points_last_row=int(at_last), at_points_initial_core=at0, at_points_later_rows=int(at_later),
                rows_ok=last_ok, rows=int(R.shape[0]), u_end=float(last_ok * o["h"]), ended_by=ended, min_rv_last_row=min_rv_last, max_sigma_last_row=smax_last,
                thickness_over_RdS=float(A * np.sqrt(2.0 / kappa)), sigma_est=float(A * np.sqrt(kappa / 2.0) * 4.0 / 3.0))


if __name__ == "__main__":
    t0 = time.time()
    out = []
    for A in (0.05, 0.1, 0.2):
        for w in (0.05, 0.1, 0.2):
            for Vw in (0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 4.0, 6.0):
                h = min(0.01, w / 8)
                r = analyse(A, Vw, w, h)
                out.append(r)
                ft = r["first_trapped"]
                print("A %.2f w %.2f Rw %.2f | m_ext %.3f R_M %.3f | trap %s | AT(core) %d | r_h %s" %
                      (A, w, Vw / 2, r["m_ext"], r["R_M"], "none" if ft is None else "u=%.2f r=%.3f" % (ft["u"], ft["r"]), r["at_points_in_core"], r["r_h_SAdS"]), flush=True)
    json.dump(out, open(os.path.join(HERE, "collapse_scan.json"), "w"), indent=1)
    # exteriors with other vacuum depths: q = -0.5 (lower de Sitter exterior), 0 (Minkowski exterior), 2 (deeper anti-de Sitter)
    outq = []
    for q in (-0.5, 0.0, 2.0):
        for w in (0.05, 0.1):
            for Vw in (1.2, 1.6, 2.0, 2.4, 3.0, 4.0, 6.0):
                r = analyse(0.1, Vw, w, min(0.01, w / 8), q=q)
                outq.append(r)
                ft = r["first_trapped"]
                print("q %.1f w %.2f Rw %.2f | m_ext %.3f R_M %.3f | trap %s | AT(core) %d | r_h %s | %s" %
                      (q, w, Vw / 2, r["m_ext"], r["R_M"], "none" if ft is None else "u=%.2f r=%.3f" % (ft["u"], ft["r"]), r["at_points_last_row"], r["r_h_SAdS"], r["ended_by"]), flush=True)
    json.dump(outq, open(os.path.join(HERE, "collapse_scan_q.json"), "w"), indent=1)
    print("saved [%.0fs]" % (time.time() - t0))
