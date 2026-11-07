# -*- coding: utf-8 -*-
"""Phase 2, step 1: junction feasibility of a pre-existing bridge around the exact cold-cycle orbits.

Geometry: the closed region of closed_cycle.py (FRW, K = +1, curvature radius a(t)) is cut at polar angle chi0 > pi/2 (beyond the equator of its S^3, the baby-universe branch of
the thin-wall analysis) and glued across a spherical shell, comoving with the interior, to a static Schwarzschild-(A)dS exterior of constant mass parameter m and vacuum density
rho_out = -q V_f.  The interior is exactly homogeneous, so the wall is at rest in the interior frame (tau = t, R = a sin chi0, R_dot = a_dot sin chi0).

Israel junction conditions (8 pi G = 1, R_dS = 1, V_f = 3):
    Gamma_in = cos(chi0)  (exact for the FRW interior: Gamma_in^2 = 1 + R_dot^2 - rho R^2/3 = 1 - sin^2 chi0),
    Gamma_out = -sqrt(F),  F = 1 + R_dot^2 - m/R + q R^2   (the sign that keeps sigma > 0 for chi0 > pi/2),
    sigma = 2 (Gamma_in - Gamma_out) / R,
    p     = -(sigma R^2)_dot / (2 R R_dot) = -(Gamma_in + sqrt F)/R - (2 R_ddot + m/R^2 + 2 q R)/(2 sqrt F),
    energy conservation on the shell, (sigma R^2)_dot + p (R^2)_dot = 0 (no flux across a comoving wall).
Requirements for a physical shell over the whole cycle: sigma > 0, sigma + p >= 0 (null energy), and a constant exterior mass m > 0.  sigma > 0 for all t is equivalent to
m < m_max = min_t R (sin^2 chi0 + R_dot^2 + q R^2).

Structural fact (not a computation): because the wall is comoving and the interior stays homogeneous, the interior equations are those of closed_cycle.py with no change; the shell
absorbs whatever the exterior requires.  The return map of the interior, and therefore its Floquet multipliers, are identical with and without the bridge.

    python bridge_junction.py          (writes bridge_junction.json)
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp

import closed_cycle as CC

HERE = os.path.dirname(os.path.abspath(__file__))


def orbit_history(o, n=6000):
    m = CC.Model(q=o["q"], rc=o["rc"], mu=o["mu"], C=o["C_over_Ccrit"] * CC.C_CRIT)
    s = 10 ** o["log10_seed_over_v"] * m.v
    a0 = CC.neck_a(m, s)
    T = o["period"]
    sol = solve_ivp(m.rhs, (0.0, T), [a0, 0.0, s, 0.0, 0.0], method="DOP853", rtol=1e-13, atol=1e-15, dense_output=True)
    t = np.linspace(0.0, T, n)
    Y = sol.sol(t)
    acc = np.array([m.rhs(t[i], Y[:, i])[1] for i in range(n)])
    return m, t, Y, acc


def shell(m_ext, q, chi0, a, ad, add):
    sn, cs = np.sin(chi0), np.cos(chi0)
    R, Rd, Rdd = a * sn, ad * sn, add * sn
    F = 1.0 + Rd ** 2 - m_ext / R + q * R ** 2
    if np.any(F <= 0):
        return None
    sq = np.sqrt(F)
    sigma = 2.0 * (cs + sq) / R
    p = -(cs + sq) / R - (2.0 * Rdd + m_ext / R ** 2 + 2.0 * q * R) / (2.0 * sq)
    return sigma, p, R


def analyse_orbit(o, label):
    m, t, Y, acc = orbit_history(o)
    a, ad = Y[0], Y[1]
    rows = []
    for q in (1.0, 0.0, -0.5):
        for chi_frac in (0.55, 0.7, 0.85, 0.95):
            chi0 = chi_frac * np.pi
            sn = np.sin(chi0)
            R = a * sn
            m_max = float(np.min(R * (sn ** 2 + (ad * sn) ** 2 + q * R ** 2)))
            rec = dict(orbit=label, q=q, chi0_over_pi=chi_frac, m_max=m_max, n_m=0, frac_nec=0.0, frac_sigma=0.0, best=None)
            if m_max > 0:
                grid = np.linspace(0.02 * m_max, 0.98 * m_max, 60)
                ok_sigma, ok_nec = 0, 0
                best = None
                for mm in grid:
                    res = shell(mm, q, chi0, a, ad, acc)
                    if res is None:
                        continue
                    sg, p, Rr = res
                    s_ok = bool(np.all(sg > 0)); n_ok = bool(np.all(sg + p >= 0))
                    ok_sigma += s_ok; ok_nec += (s_ok and n_ok)
                    margin = float(np.min(sg + p) / np.max(sg))
                    if s_ok and (best is None or margin > best["nec_margin"]):
                        best = dict(m=float(mm), sigma_min=float(sg.min()), sigma_max=float(sg.max()), p_min=float(p.min()), p_max=float(p.max()), nec_margin=margin,
                                    E_shell_range=[float((4 * np.pi * Rr ** 2 * sg).min()), float((4 * np.pi * Rr ** 2 * sg).max())])
                rec.update(n_m=len(grid), frac_sigma=ok_sigma / len(grid), frac_nec=ok_nec / len(grid), best=best)
            rows.append(rec)
    return rows


if __name__ == "__main__":
    B = json.load(open(os.path.join(HERE, "closed_stability_breather.json")))
    W = [w for w in json.load(open(os.path.join(HERE, "closed_stability.json"))) if w["mu"] == 10.0 and w["rc"] == 0.05 and abs(w["C_over_Ccrit"] - 0.9) < 1e-9]
    orbits = [("breather mu=%g rc=%g C=%.1fCc" % (b["mu"], b["rc"], b["C_over_Ccrit"]), b) for b in B] + [("wide mu=%g rc=%g C=%.1fCc" % (w["mu"], w["rc"], w["C_over_Ccrit"]), w) for w in W]
    out = []
    for label, o in orbits:
        rows = analyse_orbit(o, label)
        out += rows
        for r in rows:
            b = r["best"]
            print("%-34s q=%4.1f chi0=%.2fpi | m_max %8.4f | sigma>0 for %3.0f%% of m, NEC too for %3.0f%% | best m=%s sigma in [%s, %s] p in [%s, %s]" %
                  (label, r["q"], r["chi0_over_pi"], r["m_max"], 100 * r["frac_sigma"], 100 * r["frac_nec"], None if not b else "%.3f" % b["m"],
                   None if not b else "%.3g" % b["sigma_min"], None if not b else "%.3g" % b["sigma_max"], None if not b else "%.3g" % b["p_min"], None if not b else "%.3g" % b["p_max"]), flush=True)
    json.dump(out, open(os.path.join(HERE, "bridge_junction.json"), "w"), indent=1)
