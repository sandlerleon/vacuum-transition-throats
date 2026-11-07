# -*- coding: utf-8 -*-
"""Phase 2C, step 2: reversible energy exchange between the scalar sector and a bridge-associated reservoir, derived from an ACTION (nothing added by hand).

Because an exactly homogeneous interior is blind to a gravitational shell (bridge_moving_shell.py), the only homogeneous reduction in which the bridge can act on phi is a bulk interaction
with additional bridge-localised degrees of freedom.  Minimal choice: one real field chi (the reservoir) with the most general renormalisable Z2-symmetric action
    L = -1/2 (d phi)^2 - V_0(phi) - A(phi) rho_*  - 1/2 (d chi)^2 - 1/2 m^2 chi^2 - 1/2 g phi^2 chi^2 ,
a quadratic reservoir potential and the portal coupling g.  Then, with no free function:
    phi_ddot + 3H phi_dot + V_0' + rho_* A' + g phi chi^2 = 0,        (J_phi = -g phi chi^2)
    chi_ddot + 3H chi_dot + (m^2 + g phi^2) chi = 0,
    rho_chi + 3 ... :  the reservoir energy exchange is d/dt(chi_dot^2/2 + m^2 chi^2/2) = -3H chi_dot^2 - g phi^2 chi chi_dot ... (total energy conserved, Friedmann constraint exact),
    a_ddot = -(a/6)(rho + 3P),  a_dot^2 = a^2 rho/3 - 1,  rho = (phi_dot^2 + chi_dot^2)/2 + V_0 + m^2 chi^2/2 + g phi^2 chi^2/2 + A rho_*,
i.e. total energy-momentum is conserved and the constraint is satisfied identically.  This is the generic conservative exchange; it contains no wormhole-specific ingredient.

Experiments:
  2C-1  time-symmetric periodic orbits with nonzero reservoir amplitude (phi_dot = chi_dot = 0 at neck AND turnaround), found by 2-D shooting on (phi_0, chi_0): repeated nonsingular cycles by construction,
        reservoir energy bounded (reported as the fraction of the total energy it holds over the cycle);
  2C-2  the full 4-D return map (phi, phi_dot, chi, chi_dot at the next bounce, a fixed by the constraint): Floquet multipliers (all |mu| = 1 would be stable);
  2C-4  control: the same system with g = 0 and with chi = 0.  The coupled orbit is compared against the uncoupled one.

    python bridge_reservoir.py          (writes bridge_reservoir.json)
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve

import closed_cycle as CC

HERE = os.path.dirname(os.path.abspath(__file__))
VF = CC.VF


class TwoField:
    def __init__(self, q, rc, mu, C, mchi, g):
        self.base = CC.Model(q=q, rc=rc, mu=mu, C=C)
        self.mchi, self.g = mchi, g
        self.v = self.base.v

    def pot(self, p, c):
        return self.base.V0(p) + 0.5 * self.mchi ** 2 * c * c + 0.5 * self.g * p * p * c * c

    def rho(self, a, p, pd, c, cd):
        b = self.base
        return 0.5 * (pd * pd + cd * cd) + self.pot(p, c) + (1.0 + p * p / (2.0 * b.M2)) * b.rho_star(a)

    def rhs(self, t, y):
        a, ad, p, pd, c, cd = y
        b = self.base
        rs = b.rho_star(a)
        A = 1.0 + p * p / (2.0 * b.M2)
        K = 0.5 * (pd * pd + cd * cd)
        add = -(a / 6.0) * (4.0 * K - 2.0 * self.pot(p, c) + A * rs)
        H = ad / a
        pdd = -3.0 * H * pd - (-b.mu ** 2 * p + b.lam * p ** 3 + rs * p / b.M2 + self.g * p * c * c)
        cdd = -3.0 * H * cd - (self.mchi ** 2 + self.g * p * p) * c
        return [ad, add, pd, pdd, cd, cdd]

    def constraint(self, y):
        a, ad, p, pd, c, cd = y
        return ad * ad - (a * a * self.rho(a, p, pd, c, cd) / 3.0 - 1.0)

    def rho_P(self, yi):
        a, ad, p, pd, c, cd = yi
        K = 0.5 * (pd * pd + cd * cd)
        return self.rho(a, p, pd, c, cd), K - self.pot(p, c)

    def neck_a(self, p, pd, c, cd):
        f = lambda a: a * a * self.rho(a, p, pd, c, cd) / 3.0 - 1.0
        lo, hi = (self.base.C / 2.0) ** (1.0 / 3.0) * 1.0001, 10.0
        if f(lo) * f(hi) > 0:
            return None
        return brentq(f, lo, hi)


def integrate(m, s, t_max=6.0, stop_after=1):
    """Start from scalar state s = (p, pd, c, cd) on the constraint surface a_dot = 0; integrate to the first turnaround (stop_after = 1) or the first bounce after it (stop_after = 2).
    Returns (status, state at the stop, t)."""
    a0 = m.neck_a(*s)
    if a0 is None:
        return "no-neck", None, None
    y0 = [a0, 0.0, s[0], s[1], s[2], s[3]]
    if m.rhs(0.0, y0)[1] <= 0:
        return "not-a-bounce", None, None
    ev = lambda t, y: y[1]
    ev.direction = 0
    small = lambda t, y: y[0] - 1e-3
    small.terminal = True
    cnt = [0]

    def stop(t, y):
        return y[1]
    stop.terminal = False
    sol = solve_ivp(m.rhs, (0.0, t_max), y0, method="DOP853", rtol=1e-12, atol=1e-15, events=[ev, small], max_step=0.02)
    evs = [(t, y) for t, y in zip(sol.t_events[0], sol.y_events[0]) if t > 1e-6]
    if sol.t_events[1].size and (len(evs) < stop_after or sol.t_events[1][0] < evs[stop_after - 1][0]):
        return "singular", None, None
    if len(evs) < stop_after:
        return "no-return", None, None
    t, y = evs[stop_after - 1]
    return "ok", y, t


def residual(m, p0, c0):
    st, y, t = integrate(m, (p0, 0.0, c0, 0.0), stop_after=1)
    if st != "ok":
        return None
    return np.array([y[3], y[5]])


def find_orbits(m, n_phi=7, n_chi=7):
    """2-D shooting on (phi_0, chi_0): phi_dot = chi_dot = 0 at the first turnaround."""
    sols = []
    v = m.v
    for pf in np.geomspace(0.005, 0.5, n_phi):
        for cf in np.linspace(0.05, 1.2, n_chi):
            def F(x):
                r = residual(m, x[0] * v, x[1])
                return np.array([10.0, 10.0]) if r is None else r
            try:
                x, info, ier, msg = fsolve(F, [pf, cf], full_output=True, xtol=1e-12)
            except Exception:           # noqa: BLE001
                continue
            if ier == 1 and np.linalg.norm(F(x)) < 1e-9 and x[1] > 0.01 and x[0] > 1e-4:
                if not any(abs(x[0] - s[0]) < 1e-6 and abs(x[1] - s[1]) < 1e-5 for s in sols):
                    sols.append((float(x[0]), float(x[1])))
    return sols


def return_map(m, s):
    st, y, t = integrate(m, tuple(s), stop_after=2)
    return None if st != "ok" else (np.array([y[2], y[3], y[4], y[5]]), t)


def floquet(m, s, h=1e-7):
    base = return_map(m, s)
    if base is None:
        return None
    J = np.zeros((4, 4))
    sc = np.array([m.v, m.v, 1.0, 1.0])
    for k in range(4):
        d = np.zeros(4); d[k] = h * sc[k]
        fp, fm = return_map(m, np.array(s) + d), return_map(m, np.array(s) - d)
        if fp is None or fm is None:
            return None
        J[:, k] = (fp[0] - fm[0]) / (2 * d[k])
    mult = np.linalg.eigvals(J)
    return dict(period=float(base[1]), state_return_error=float(np.linalg.norm((base[0] - np.array(s)) / sc)), multipliers_abs=sorted([float(abs(x)) for x in mult], reverse=True),
                multipliers=[[float(x.real), float(x.imag)] for x in mult], det=float(np.linalg.det(J)), max_abs=float(np.max(np.abs(mult))))


def reservoir_fraction(m, s):
    """Fraction of the total energy density held by the reservoir (chi_dot^2/2 + m^2 chi^2/2 + g phi^2 chi^2/2) over one period; min, max."""
    st, y, t = integrate(m, tuple(s), stop_after=2)
    a0 = m.neck_a(*s)
    y0 = [a0, 0.0, s[0], s[1], s[2], s[3]]
    sol = solve_ivp(m.rhs, (0.0, t), y0, method="DOP853", rtol=1e-12, atol=1e-15, dense_output=True)
    tt = np.linspace(0.0, t, 800)
    Y = sol.sol(tt)
    fr = []
    for i in range(len(tt)):
        a, ad, p, pd, c, cd = Y[:, i]
        er = 0.5 * cd * cd + 0.5 * m.mchi ** 2 * c * c + 0.5 * m.g * p * p * c * c
        fr.append(er / m.rho(a, p, pd, c, cd))
    cons = max(abs(m.constraint(Y[:, i])) for i in range(0, len(tt), 20))
    return float(min(fr)), float(max(fr)), float(cons), float(Y[0].min()), float(Y[0].max())


if __name__ == "__main__":
    import sys
    B = json.load(open(os.path.join(HERE, "closed_stability_breather.json")))
    o = [b for b in B if b["mu"] == 20.0 and b["rc"] == 0.2][0]
    masses = [float(x) for x in sys.argv[1:]] or [3.0, 6.0, 9.0, 13.0, 20.0]
    out = []
    for mchi in masses:
        for g in (30.0, 100.0, 300.0, 1000.0, 3000.0):
            m = TwoField(o["q"], o["rc"], o["mu"], o["C_over_Ccrit"] * CC.C_CRIT, mchi, g)
            orbs = find_orbits(m, n_phi=5, n_chi=5)
            print("m_chi=%g g=%g: %d orbit(s) with chi0>0" % (mchi, g, len(orbs)), flush=True)
            for (pf, c0) in orbs:
                s = (pf * m.v, 0.0, c0, 0.0)
                fl = floquet(m, s)
                if fl is None:
                    continue
                fr = reservoir_fraction(m, s)
                rec = dict(mchi=mchi, g=g, phi0_over_v=pf, chi0=c0, reservoir_frac_min=fr[0], reservoir_frac_max=fr[1], constraint_max=fr[2], a_min=fr[3], a_max=fr[4], **fl)
                out.append(rec)
                print("m_chi=%-4g g=%-6g phi0/v=%.4f chi0=%.4f | T=%.4f ret.err %.1e | |mult| %s det %.3f | reservoir %.3f-%.3f | constraint %.1e" %
                      (mchi, g, pf, c0, fl["period"], fl["state_return_error"], ",".join("%.3g" % x for x in fl["multipliers_abs"]), fl["det"], fr[0], fr[1], fr[2]), flush=True)
    tag = "_".join("%g" % x for x in masses)
    json.dump(out, open(os.path.join(HERE, "bridge_reservoir_%s.json" % tag), "w"), indent=1)
