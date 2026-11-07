# -*- coding: utf-8 -*-
"""Phase 2C, step 1: a MOVING thin shell (not comoving) between the closed matter-coupled FRW region of closed_cycle.py and a Schwarzschild-(A)dS exterior.

Units 8 pi G = 1, R_dS = 1.  Interior (exactly homogeneous): ds^2 = -dt^2 + a(t)^2 (dchi^2 + sin^2 chi dOmega^2), a(t), phi(t) from closed_cycle.Model.  Shell worldline (t(tau), chi(tau)) with
rapidity beta relative to the comoving fluid:   dt/dtau = cosh beta,  a dchi/dtau = sinh beta.  Then
    R = a sin(chi),  R' = a_dot sin(chi) cosh(beta) + cos(chi) sinh(beta),
    Gamma_in = n.grad R = cos(chi) cosh(beta) + a_dot sin(chi) sinh(beta)      [Gamma_in^2 = 1 + R'^2 - rho R^2/3 exactly, by the Friedmann constraint],
    K_in^tau_tau = beta' + H sinh(beta)   (proper acceleration of the worldline along the outward normal; geodesic: beta' = -H sinh beta),
    Gamma_out = -sqrt(F),  F = 1 + R'^2 - m/R + q R^2,   K_out^tau_tau = Gamma_out'/R' = -(2 R'' + m/R^2 + 2 q R)/(2 sqrt F).
Israel conditions (S^a_b = -([K^a_b] - delta^a_b [K])):
    sigma = 2 (Gamma_in - Gamma_out)/R,        p = [K^tau_tau] + [K^theta_theta] = K_out^tau_tau - K_in^tau_tau + (Gamma_out - Gamma_in)/R.
With a shell equation of state p = w sigma this is a linear equation for beta' (R'' = R0 + Gamma_in beta'), so the shell obeys a closed second-order ODE; nothing is added by hand.
The shell conservation law with flux (Gauss-Codazzi) is  (sigma R^2)' + p (R^2)' = -R^2 (rho + P) cosh(beta) sinh(beta)  (interior fluid, exterior vacuum; sign fixed by the numerical check below);
the code checks this identity independently of the equation of motion (agreement ~1e-6 with the opposite-sign form disagreeing at O(1)).

STRUCTURAL FACT (proved by inspection of the equations): the interior equations for (a, phi, rho_r) contain no shell variable, so for an exactly homogeneous interior the scalar return map
and its Floquet multipliers are IDENTICAL for every shell trajectory.  A moving gravitational shell can exchange energy with the shell/exterior only; it cannot act on phi.

    python bridge_moving_shell.py          (writes bridge_moving_shell.json)
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp

import closed_cycle as CC
import bridge_junction as BJ

HERE = os.path.dirname(os.path.abspath(__file__))
EPS = -1.0


class Ph1(CC.Model):
    """Phase 1 single-field interior with the generic interface used by the shell: state [a, ad, phi, phid, rr]."""
    def rho_P(self, yi):
        a, ad, p, pd, rr = yi
        return self.rho(a, p, pd, rr), 0.5 * pd * pd - self.V0(p)


def shell_state(model, y, m_ext, q, w):
    """Quantities of the shell at state y = [t, chi, beta, *interior]; the interior object supplies rhs(t, yi) (yi[0] = a, yi[1] = a_dot) and rho_P(yi); returns a dict (including beta')."""
    t, chi, beta = y[:3]
    yi = list(y[3:])
    a, ad = yi[0], yi[1]
    add = model.rhs(t, yi)[1]
    s, c, ch, sh = np.sin(chi), np.cos(chi), np.cosh(beta), np.sinh(beta)
    H = ad / a
    R = a * s
    Rp = ad * s * ch + c * sh
    G_in = c * ch + ad * s * sh
    R0 = add * s * ch ** 2 + H * c * sh * ch - (s / a) * sh ** 2
    F = 1.0 + Rp ** 2 - m_ext / R + q * R ** 2
    if F <= 0:
        return None
    sq = np.sqrt(F)
    G_out = EPS * sq
    sigma = 2.0 * (G_in - G_out) / R
    A0 = EPS * (2.0 * R0 + m_ext / R ** 2 + 2.0 * q * R) / (2.0 * sq)
    A1 = EPS * G_in / sq - 1.0
    bp = ((w + 0.5) * sigma + H * sh - A0) / A1
    Rpp = R0 + G_in * bp
    K_out = EPS * (2.0 * Rpp + m_ext / R ** 2 + 2.0 * q * R) / (2.0 * sq)
    K_in = bp + H * sh
    pres = K_out - K_in + (G_out - G_in) / R
    rho, P = model.rho_P(yi)
    return dict(R=R, Rp=Rp, G_in=G_in, G_out=G_out, sigma=sigma, p=pres, bp=bp, Rpp=Rpp, rho=rho, P=P, F=F, add=add)


def rhs_factory(model, m_ext, q, w):
    def rhs(tau, y):
        st = shell_state(model, y, m_ext, q, w)
        if st is None:
            return [np.nan] * len(y)
        t, chi, beta = y[:3]
        ch, sh = np.cosh(beta), np.sinh(beta)
        d = model.rhs(t, list(y[3:]))
        return [ch, sh / y[3], st["bp"]] + [ch * x for x in d]
    return rhs


def evolve(o, chi0, m_frac, q, w, n_per=1.0, npts=4000):
    model = Ph1(q=o["q"], rc=o["rc"], mu=o["mu"], C=o["C_over_Ccrit"] * CC.C_CRIT)
    seed = 10 ** o["log10_seed_over_v"] * model.v
    a0 = CC.neck_a(model, seed)
    m, t, Y, acc = BJ.orbit_history(o, n=4000)
    sn = np.sin(chi0)
    R = Y[0] * sn
    m_max = float(np.min(R * (sn ** 2 + (Y[1] * sn) ** 2 + q * R ** 2)))
    m_ext = m_frac * m_max
    y0 = [0.0, chi0, 0.0, a0, 0.0, seed, 0.0, 0.0]
    T = o["period"] * n_per
    ev = lambda tau, y: y[3] - 1e-3
    ev.terminal = True
    sol = solve_ivp(rhs_factory(model, m_ext, q, w), (0.0, 3.0 * T), y0, method="DOP853", rtol=1e-11, atol=1e-13, max_step=0.002, dense_output=True,
                    events=[ev, lambda tau, y: y[0] - T])
    return model, sol, m_ext, m_max, T


def diagnostics(model, sol, m_ext, q, w, T):
    """Shell history over one interior period (t from 0 to T): sigma, p, NEC, beta, flux identity residual."""
    tau_end = sol.t_events[1][0] if sol.t_events[1].size else sol.t[-1]
    tau = np.linspace(0.0, tau_end, 6000)
    Y = sol.sol(tau)
    S = [shell_state(model, Y[:, i], m_ext, q, w) for i in range(len(tau))]
    if any(s is None for s in S):
        k = next(i for i, s in enumerate(S) if s is None)
        tau, Y, S = tau[:k], Y[:, :k], S[:k]
    sig = np.array([s["sigma"] for s in S]); pr = np.array([s["p"] for s in S]); R = np.array([s["R"] for s in S])
    f = np.array([(s["rho"] + s["P"]) * np.cosh(Y[2, i]) * np.sinh(Y[2, i]) for i, s in enumerate(S)])
    lhs = np.gradient(sig * R ** 2, tau) + pr * np.gradient(R ** 2, tau)
    rhsv = -R ** 2 * f
    sel = slice(20, -20)
    scale = max(np.max(np.abs(rhsv[sel])), np.max(np.abs(lhs[sel])), 1e-12)
    resid = float(np.max(np.abs(lhs[sel] - rhsv[sel])) / scale)
    resid_flip = float(np.max(np.abs(lhs[sel] + rhsv[sel])) / scale)   # wrong-sign comparison
    return dict(reached_period=bool(sol.t_events[1].size), t_end=float(Y[0, -1]), sigma_pos=float(np.mean(sig > 0)), nec=float(np.mean((sig > 0) & (sig + pr >= 0))),
                max_abs_beta=float(np.max(np.abs(Y[2]))), R_min=float(R.min()), R_max=float(R.max()), flux_identity_resid=resid, flux_identity_resid_wrong_sign=resid_flip,
                dE_shell=float(sig[-1] * R[-1] ** 2 - sig[0] * R[0] ** 2), sigmaR2_range=[float((sig * R ** 2).min()), float((sig * R ** 2).max())])


def one_period_map(model, yi0, chi0, beta0, m_ext, q, w, T):
    y0 = [0.0, chi0, beta0] + list(yi0)
    ev = lambda tau, y: y[0] - T
    ev.terminal = True
    bad = lambda tau, y: y[3] - 1e-3
    bad.terminal = True
    sol = solve_ivp(rhs_factory(model, m_ext, q, w), (0.0, 4.0 * T), y0, method="DOP853", rtol=1e-12, atol=1e-14, max_step=0.002, events=[ev, bad])
    if not sol.t_events[0].size:
        return None
    return sol.y_events[0][0]


def periodic_shell_generic(model, yi0, T, chi_frac_range, m_frac, q, w, ngrid=39):
    """Time-symmetric periodic shell around a time-symmetric periodic interior (neck at t = 0): beta = 0 at t = 0 and t = T/2; shoot on chi0.  Returns the periodic shells found, their
    full-period 2-D return-map Jacobian and multipliers, and the shell history diagnostics."""
    from scipy.optimize import brentq
    sol = solve_ivp(model.rhs, (0.0, T), list(yi0), method="DOP853", rtol=1e-13, atol=1e-15, dense_output=True)
    tt = np.linspace(0.0, T, 4000)
    Yi = sol.sol(tt)

    def m_ext_of(chi0):
        sn = np.sin(chi0)
        R = Yi[0] * sn
        return m_frac * float(np.min(R * (sn ** 2 + (Yi[1] * sn) ** 2 + q * R ** 2)))

    def beta_half(chi0):
        r = one_period_map(model, yi0, chi0, 0.0, m_ext_of(chi0), q, w, 0.5 * T)
        return np.nan if r is None else float(r[2])

    grid = np.linspace(chi_frac_range[0], chi_frac_range[1], ngrid) * np.pi
    vals = [beta_half(c) for c in grid]
    out = []
    for i in range(len(grid) - 1):
        if not (np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] * vals[i + 1] < 0):
            continue
        try:
            chi0 = brentq(beta_half, grid[i], grid[i + 1], xtol=1e-12)
        except Exception:           # noqa: BLE001
            continue
        mext = m_ext_of(chi0)

        def full(c, b):
            r = one_period_map(model, yi0, c, b, mext, q, w, T)
            return None if r is None else np.array([r[1], r[2]])
        base = full(chi0, 0.0)
        if base is None:
            continue
        h = 1e-6
        J = np.zeros((2, 2))
        ok = True
        for k, (dc, db) in enumerate(((h, 0.0), (0.0, h))):
            fp, fm = full(chi0 + dc, db), full(chi0 - dc, -db)
            if fp is None or fm is None:
                ok = False
                break
            J[:, k] = (fp - fm) / (2 * h)
        # history of this periodic shell
        ev = lambda tau, y: y[0] - T
        ev.terminal = True
        evs = lambda tau, y: y[3] - 1e-3
        evs.terminal = True
        so = solve_ivp(rhs_factory(model, mext, q, w), (0.0, 4.0 * T), [0.0, chi0, 0.0] + list(yi0), method="DOP853", rtol=1e-12, atol=1e-14, max_step=0.002, dense_output=True, events=[evs, ev])
        d = diagnostics(model, so, mext, q, w, T)
        out.append(dict(chi0_over_pi=float(chi0 / np.pi), m_ext=float(mext), return_chi_minus_chi0=float(base[0] - chi0), return_beta=float(base[1]),
                        jac=J.tolist() if ok else None, multipliers_abs=[float(abs(x)) for x in np.linalg.eigvals(J)] if ok else None, det=float(np.linalg.det(J)) if ok else None, diag=d))
    return out


def periodic_shell(o, chi_frac, m_frac, q, w):
    model = Ph1(q=o["q"], rc=o["rc"], mu=o["mu"], C=o["C_over_Ccrit"] * CC.C_CRIT)
    seed = 10 ** o["log10_seed_over_v"] * model.v
    a0 = CC.neck_a(model, seed)
    return periodic_shell_generic(model, [a0, 0.0, seed, 0.0, 0.0], o["period"], (0.52, 0.9), m_frac, q, w)


if __name__ == "__main__":
    B = json.load(open(os.path.join(HERE, "closed_stability_breather.json")))
    o = [b for b in B if b["mu"] == 20.0 and b["rc"] == 0.2][0]
    out = []
    for q in (1.0, 0.0):
        for chi_frac in (0.55, 0.7, 0.85):
            for m_frac in (0.2, 0.5, 0.8):
                for w in (-1.0, -0.5, 0.0, 1.0 / 3.0, 1.0):
                    model, sol, m_ext, m_max, T = evolve(o, chi_frac * np.pi, m_frac, q, w)
                    d = diagnostics(model, sol, m_ext, q, w, T)
                    d.update(q=q, chi0_over_pi=chi_frac, m_frac=m_frac, w=w, m_ext=m_ext)
                    out.append(d)
                    print("q=%g chi0=%.2fpi m/mmax=%.1f w=%+.2f | period %s t_end=%.3f | sigma>0 %.2f NEC %.2f | |beta|max %.2f | R [%.3f,%.3f] | flux-identity resid %.1e (wrong sign %.1e)" %
                          (q, chi_frac, m_frac, w, "yes" if d["reached_period"] else "NO ", d["t_end"], d["sigma_pos"], d["nec"], d["max_abs_beta"], d["R_min"], d["R_max"],
                           d["flux_identity_resid"], d["flux_identity_resid_wrong_sign"]), flush=True)
    per = []
    for q in (1.0, 0.0):
        for m_frac in (0.2, 0.5, 0.8):
            for w in (-1.0, -0.5, 0.0, 1.0 / 3.0, 1.0):
                for r in periodic_shell(o, 0.7, m_frac, q, w):
                    r.update(q=q, m_frac=m_frac, w=w)
                    per.append(r)
                    print("periodic shell q=%g m/mmax=%.1f w=%+.2f chi0=%.4fpi | return dchi=%.1e beta=%.1e | |mult|=%s det=%.4f" % (q, m_frac, w, r["chi0_over_pi"], r["return_chi_minus_chi0"], r["return_beta"],
                          None if r["multipliers_abs"] is None else ",".join("%.3f" % x for x in r["multipliers_abs"]), r["det"] if r["det"] is not None else float("nan")), flush=True)
    json.dump(dict(scan=out, periodic=per), open(os.path.join(HERE, "bridge_moving_shell.json"), "w"), indent=1)
