# -*- coding: utf-8 -*-
"""Phase 2D, steps 3-8: the two-ended geometry around the stable cold cycle, and its causal structure.

GEOMETRY (specified, not inferred).  Parent spacetime: the maximally extended Schwarzschild solution of mass parameter m (f(R) = 1 - m/R, 2GM = m), four Kruskal regions
    I (right exterior), II (black-hole interior), III (left exterior), IV (white-hole interior)    with  U V = -e^{2 kappa r*}  in I and III.
In region III (future-directed Killing time s) U = e^{kappa (r* + s)}, V = -e^{kappa (r* - s)}.  Cut region III along a timelike spherical shell S = {R = R_s(tau)} and DISCARD the part of
III with R > R_s (the left asymptotic end).  Glue the discarded side to a closed FRW region (the stable cold cycle, a cap of S^3 of polar angle chi_0 > pi/2) across S.  The kept exterior is
    E = I  U  II  U  IV  U  {points of III with R < R_s},        the junction orientation is Gamma_out = -sqrt(F) < 0 (area decreases away from the cap, toward the horizon).
Along a spatial slice through the cap the areal radius R rises from 0 to a maximum a (the equator of S^3), decreases through the shell to the throat at r_h = m (the bifurcation 2-sphere, a minimum),
and increases again to infinity in region I.  So: two ends (the closed cyclic region and the asymptotically flat region I), joined by a Schwarzschild throat that is hidden behind horizons.

Israel conditions are those of bridge_moving_shell.py.  Interior = the INDEPENDENT Cycle model (phase2d_cycle_validation.py) at C/Cc = 0.35.

    python phase2d_bridge.py          (writes phase2d_bridge.json and figure FigD1.png)
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

import bridge_moving_shell as MS
from phase2d_shell_search import Interior, setup

HERE = os.path.dirname(os.path.abspath(__file__))


def ext_geometry(m, q=0.0):
    f = lambda R: 1.0 - m / R + q * R * R
    rh = brentq(f, 1e-9, 50.0)
    kappa = 0.5 * (m / rh ** 2 + 2.0 * q * rh)
    Rref = 3.0 * rh

    def rstar(R):
        if abs(q) < 1e-14:
            return R + m * np.log(abs(R / m - 1.0))
        return quad(lambda z: 1.0 / f(z), Rref, R, limit=200)[0]
    return f, rh, kappa, rstar


def shell_history(I, yi0, shell, N, npts_per=800):
    """Evolve interior + shell for N interior periods starting from the periodic-shell data; return dense arrays."""
    q, w, m_ext, chi0 = shell["q"], shell["w"], shell["m_ext"], shell["chi0_over_pi"] * np.pi
    T = shell["T"]
    ev = lambda tau, y: y[0] - N * T
    ev.terminal = True
    evs = lambda tau, y: y[3] - 1e-3
    evs.terminal = True
    sol = solve_ivp(MS.rhs_factory(I, m_ext, q, w), (0.0, 4.0 * N * T), [0.0, chi0, 0.0] + list(yi0), method="DOP853", rtol=1e-12, atol=1e-14, max_step=0.002, dense_output=True, events=[evs, ev])
    tau_end = sol.t_events[1][0] if sol.t_events[1].size else sol.t[-1]
    tau = np.linspace(0.0, tau_end, N * npts_per + 1)
    Y = sol.sol(tau)
    S = [MS.shell_state(I, Y[:, i], m_ext, q, w) for i in range(len(tau))]
    return sol, tau, Y, S, bool(sol.t_events[1].size)


def analyse_shell(I, yi0, shell, N=8):
    q, w, m = shell["q"], shell["w"], shell["m_ext"]
    f, rh, kappa, rstar = ext_geometry(m, q)
    sol, tau, Y, S, done = shell_history(I, yi0, shell, N)
    t, chi, beta = Y[0], Y[1], Y[2]
    a = Y[3]
    R = np.array([s["R"] for s in S]); Rp = np.array([s["Rp"] for s in S]); Gin = np.array([s["G_in"] for s in S]); Gout = np.array([s["G_out"] for s in S])
    sig = np.array([s["sigma"] for s in S]); pr = np.array([s["p"] for s in S]); F = np.array([s["F"] for s in S]); rho = np.array([s["rho"] for s in S])
    rec = dict(q=q, w=w, m=m, m_frac=shell["m_frac"], chi0_over_pi=shell["chi0_over_pi"], completed_periods=bool(done), r_h=float(rh), kappa=float(kappa), throat_area_over_4pi=float(rh ** 2),
               R_shell_min=float(R.min()), R_shell_max=float(R.max()), equator_a_min=float(a.min()), equator_a_max=float(a.max()),
               ordering_equator_gt_shell_gt_throat=bool(np.all(a > R) and np.all(R > rh)), shell_over_throat_min=float(R.min() / rh),
               sigma_pos=float(np.mean(sig > 0)), nec_fraction=float(np.mean((sig > 0) & (sig + pr >= 0))), sigma_range=[float(sig.min()), float(sig.max())], p_range=[float(pr.min()), float(pr.max())])
    # null expansions on both sides: l_pm(R) = R' +- Gamma
    rec["interior_trapped_fraction"] = float(np.mean(Gin ** 2 < Rp ** 2))        # (R'+Gamma)(R'-Gamma) > 0  <=>  trapped
    rec["exterior_trapped_fraction"] = float(np.mean(Gout ** 2 < Rp ** 2))
    rec["interior_Gamma_sign"] = [float(Gin.min()), float(Gin.max())]
    rec["exterior_Gamma_sign"] = [float(Gout.min()), float(Gout.max())]
    # equator is trapped or marginal: R = a, (grad R)^2 = 1 - rho a^2/3 = -a_dot^2 <= 0
    rec["max_equator_ad2"] = float(np.max(Y[4] ** 2))
    # exterior Killing time and Kruskal coordinates of the shell worldline
    dsdtau = np.sqrt(F) / np.array([f(r) for r in R])
    s = np.concatenate([[0.0], np.cumsum(0.5 * (dsdtau[1:] + dsdtau[:-1]) * np.diff(tau))])
    rs = np.array([rstar(r) for r in R])
    lnU = kappa * (rs + s); lnmV = kappa * (rs - s)          # ln U and ln(-V): Kruskal coordinates in log form (kappa is ~1/(2m), so U, V overflow)
    rec["shell_outside_horizon_always"] = bool(np.all(R > rh))      # <=> the shell worldline stays in region III
    rec["killing_time_per_interior_period"] = float(s[-1] / N)
    # radial null rays from the cap centre (chi = 0): chi_ray(t) = Phi(t) - Phi(t_emit), Phi = int dt / a along the interior history
    Phi = np.concatenate([[0.0], np.cumsum(0.5 * ((np.cosh(beta[1:]) / a[1:]) + (np.cosh(beta[:-1]) / a[:-1])) * np.diff(tau))])
    rays = []
    nem = 24
    temis = np.linspace(0.0, shell["T"], nem, endpoint=False)
    for te in temis:
        Phi_te = np.interp(te, t, Phi)
        h = chi + Phi_te - Phi                     # chi_shell - chi_ray
        ic = None
        for i in range(len(tau) - 1):
            if t[i] > te and h[i] > 0 and h[i + 1] <= 0:
                ic = i; break
        if ic is None:
            rays.append(dict(t_emit=float(te), crossed=False)); continue
        fr = h[ic] / (h[ic] - h[ic + 1])
        lerp = lambda arr: float(arr[ic] + fr * (arr[ic + 1] - arr[ic]))
        rays.append(dict(t_emit=float(te), crossed=True, t_cross=lerp(t), R_cross=lerp(R), s_cross=lerp(s), rstar_cross=lerp(rs), ln_U=lerp(lnU), ln_minus_V=lerp(lnmV), delay_in_periods=float((lerp(t) - te) / shell["T"])))
    rec["rays"] = rays
    rec["rays_crossed"] = sum(1 for r in rays if r["crossed"])
    # Causal statement (exact, from the Kruskal structure): in the kept exterior every future-directed causal curve from a point of the shell (U > 0) has U non-decreasing, so U >= U_shell > 0;
    # region I is {U < 0}, hence J^+(cap) and region I are disjoint, and region I cannot send a signal to the cap (J^-(cap) lies in U > 0 as well? no: J^-(cap) is in regions III/IV).
    rec["min_ln_U_on_shell"] = float(lnU.min())
    rec["every_ray_has_V_negative_at_crossing"] = bool(all(r["crossed"] for r in rays))
    rec["Y_diag"] = dict(tau=tau[::10].tolist(), t=t[::10].tolist(), R=R[::10].tolist(), a=a[::10].tolist(), lnU=lnU[::10].tolist(), lnmV=lnmV[::10].tolist(), s=s[::10].tolist(), rstar=rs[::10].tolist(), sigma=sig[::10].tolist(), p=pr[::10].tolist())
    return rec


def long_run(I, yi0, shell, N, pert, seed=0):
    """Combined evolution with a random perturbation of size `pert` (relative) of the interior neck state and of the shell data (chi0, beta0); compares the interior against the interior alone."""
    rng = np.random.default_rng(seed)
    q, w, m_ext = shell["q"], shell["w"], shell["m_ext"]
    T = shell["T"]
    chi0 = shell["chi0_over_pi"] * np.pi
    d = rng.normal(size=4); d /= np.linalg.norm(d)
    p0 = yi0[2] * (1.0 + pert * d[0])
    a0 = I.c.neck_a(p0, 0.0)
    yi = [a0, 0.0, p0, 0.0]
    chi_p, beta_p = chi0 + pert * d[1], pert * d[2]
    ev = lambda tau, y: y[0] - N * T
    ev.terminal = True
    evs = lambda tau, y: y[3] - 1e-3
    evs.terminal = True
    sol = solve_ivp(MS.rhs_factory(I, m_ext, q, w), (0.0, 4.0 * N * T), [0.0, chi_p, beta_p] + yi, method="DOP853", rtol=1e-11, atol=1e-13, max_step=0.004, dense_output=True, events=[evs, ev])
    reached = bool(sol.t_events[1].size)
    tau_end = sol.t_events[1][0] if reached else sol.t[-1]
    tau = np.linspace(0.0, tau_end, 40 * N + 1)
    Y = sol.sol(tau)
    S = [MS.shell_state(I, Y[:, i], m_ext, q, w) for i in range(len(tau))]
    bad = [i for i, s in enumerate(S) if s is None]
    good = [s for s in S if s is not None]
    f, rh, kappa, rstar = ext_geometry(m_ext, q)
    R = np.array([s["R"] for s in good]); sig = np.array([s["sigma"] for s in good]); pr = np.array([s["p"] for s in good])
    # interior alone with the same data
    ev2 = lambda t_, y: y[0] - 1e-3
    ev2.terminal = True
    base = solve_ivp(I.rhs, (0.0, N * T), yi, method="DOP853", rtol=1e-11, atol=1e-13, events=[ev2])
    t_final = Y[0, -1]
    ref = solve_ivp(I.rhs, (0.0, t_final), yi, method="DOP853", rtol=1e-11, atol=1e-13)
    interior_diff = float(np.max(np.abs(Y[3:, -1] - ref.y[:, -1])))
    return dict(pert=pert, seed=seed, periods_reached=float(Y[0, -1] / T), completed=reached and not bad, invalid_points=len(bad), R_min=float(R.min()), R_max=float(R.max()), R_over_rh_min=float(R.min() / rh),
                sigma_min=float(sig.min()), nec_fraction=float(np.mean((sig > 0) & (sig + pr >= 0))), max_abs_beta=float(np.max(np.abs(Y[2]))), interior_vs_alone_diff=interior_diff,
                a_min=float(Y[3].min()), a_max=float(Y[3].max()))


def neck_portrait(c, p_fix, offsets, nit=150):
    """Iterates of the neck map around the elliptic fixed point (for the figure)."""
    out = []
    for off in offsets:
        z = np.array([p_fix + off[0] * c.v, off[1] * c.v]); pts = []
        for k in range(nit):
            r = c.P(z[0], z[1])
            if r is None:
                break
            z = np.array([r[0], r[1]]); pts.append([(z[0] - p_fix) / c.v, z[1] / c.v])
        out.append(pts)
    return out


if __name__ == "__main__":
    I, yi0, T = setup()
    shells = json.load(open(os.path.join(HERE, "phase2d_shells_q0.json")))["shells"]
    for s in shells:
        s["T"] = T
    results = dict(T=T, yi0=yi0, analysis=[], long=[])
    for s in shells:
        r = analyse_shell(I, yi0, s, N=8)
        results["analysis"].append(r)
        print("m/mmax=%.2f chi0=%.4fpi m=%.5f r_h=%.5f | shell R [%.3f,%.3f] (>= %.0f r_h), equator a [%.3f,%.3f] | order a>R>r_h: %s | sigma>0 %.2f NEC %.2f | trapped int/ext %.2f/%.2f | rays crossed %d/24 | shell outside horizon: %s" %
              (r["m_frac"], r["chi0_over_pi"], r["m"], r["r_h"], r["R_shell_min"], r["R_shell_max"], r["shell_over_throat_min"], r["equator_a_min"], r["equator_a_max"], r["ordering_equator_gt_shell_gt_throat"],
               r["sigma_pos"], r["nec_fraction"], r["interior_trapped_fraction"], r["exterior_trapped_fraction"], r["rays_crossed"], r["shell_outside_horizon_always"]), flush=True)
    # long runs on the most conservative shell (smallest |beta|) and the largest one
    for idx in (0, len(shells) - 1):
        for pert in (1e-6, 1e-4, 1e-3, 1e-2):
            for seed in range(2):
                lr = long_run(I, yi0, shells[idx], 60, pert, seed)
                lr["shell_index"] = idx
                results["long"].append(lr)
                print("  shell %d pert %.0e seed %d: periods %.1f completed %s invalid %d R/r_h>=%.1f sigma_min %.2e NEC %.3f |beta|max %.2f interior-vs-alone diff %.1e a in [%.3f,%.3f]" %
                      (idx, pert, seed, lr["periods_reached"], lr["completed"], lr["invalid_points"], lr["R_over_rh_min"], lr["sigma_min"], lr["nec_fraction"], lr["max_abs_beta"], lr["interior_vs_alone_diff"], lr["a_min"], lr["a_max"]), flush=True)
    # neck-map portrait
    c = I.c
    p_fix = yi0[2]
    results["portrait"] = neck_portrait(c, p_fix, [(0.01, 0.0), (0.0, 0.01), (0.02, 0.0), (0.0, 0.03), (0.04, 0.0), (0.0, 0.06), (-0.05, 0.0), (0.06, 0.02), (0.0, 0.09)], nit=120)
    json.dump(results, open(os.path.join(HERE, "phase2d_bridge.json"), "w"))
    print("saved")
