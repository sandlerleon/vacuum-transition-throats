# -*- coding: utf-8 -*-
"""Phase 2D, step 1-2: INDEPENDENT reproduction and classification of the stable cold cycle (does not import closed_cycle.py, closed_orbits.py or closed_stability.py).

Model (units 8 pi G = 1, R_dS = 1, V_f = 3): closed FRW, curvature +1, scale factor a; scalar V(phi) = V_f - mu^2 phi^2/2 + lam phi^4/4, v^2 = 4(1+q)V_f/mu^2, lam = mu^2/v^2;
dust rho_* = 3 C / a^3 with the coupling A(phi) = 1 + phi^2/(2 M^2), M^2 = rho_crit/mu^2, rho_crit = rc V_f; C = x C_c, C_c = 2/(3 sqrt 3).
Equations written independently here: (a, b = a_dot, phi, w = phi_dot)
    b' = -(a/6)(2 w^2 - 2 V + A rho_*),    w' = -3 (b/a) w - V'(phi) - rho_* phi/M^2,   Friedmann: b^2 = a^2 rho/3 - 1.
Different from Phase 1: (i) integrator LSODA/Radau/DOP853 chosen per test, (ii) the clock is cosmic time for the flow and CONFORMAL time for an independent orbit solve, (iii) the stability
is computed three ways: finite-difference Jacobian of the neck map, the 4-D monodromy of the variational equations, and a conformal-time neck map.

    python phase2d_cycle_validation.py          (writes phase2d_cycle_validation.json)
"""
import json
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve

HERE = os.path.dirname(os.path.abspath(__file__))
VF = 3.0
CC = 2.0 / (3.0 * np.sqrt(3.0))


class Cycle:
    def __init__(self, x, q=1.0, rc=0.2, mu=20.0):
        self.x, self.q, self.rc, self.mu = x, q, rc, mu
        self.C = x * CC
        self.M2 = rc * VF / mu ** 2
        self.v2 = 4.0 * (1.0 + q) * VF / mu ** 2
        self.lam = mu ** 2 / self.v2
        self.v = np.sqrt(self.v2)

    def V(self, p):
        return VF - 0.5 * self.mu ** 2 * p * p + 0.25 * self.lam * p ** 4

    def Vp(self, p):
        return -self.mu ** 2 * p + self.lam * p ** 3

    def rho(self, a, p, w):
        return 0.5 * w * w + self.V(p) + (1.0 + p * p / (2.0 * self.M2)) * 3.0 * self.C / a ** 3

    def f(self, t, y):
        a, b, p, w = y
        Q = 3.0 * self.C / a ** 3
        A = 1.0 + p * p / (2.0 * self.M2)
        return [b, -(a / 6.0) * (2.0 * w * w - 2.0 * self.V(p) + A * Q), w, -3.0 * b / a * w - self.Vp(p) - Q * p / self.M2]

    def jac(self, y):
        a, b, p, w = y
        Q = 3.0 * self.C / a ** 3
        A = 1.0 + p * p / (2.0 * self.M2)
        X = 2.0 * w * w - 2.0 * self.V(p) + A * Q
        Vpp = -self.mu ** 2 + 3.0 * self.lam * p * p
        return np.array([
            [0.0, 1.0, 0.0, 0.0],
            [-X / 6.0 + A * Q / 2.0, 0.0, -(a / 6.0) * (-2.0 * self.Vp(p) + p * Q / self.M2), -2.0 * a * w / 3.0],
            [0.0, 0.0, 0.0, 1.0],
            [3.0 * b * w / a ** 2 + 3.0 * Q * p / (self.M2 * a), -3.0 * w / a, -Vpp - Q / self.M2, -3.0 * b / a]])

    def constraint(self, y):
        a, b, p, w = y
        return b * b - (a * a * self.rho(a, p, w) / 3.0 - 1.0)

    def neck_a(self, p, w):
        g = lambda a: a * a * self.rho(a, p, w) / 3.0 - 1.0
        lo, hi = (self.C / 2.0) ** (1.0 / 3.0) * 1.0001, 10.0
        return brentq(g, lo, hi, xtol=1e-15) if g(lo) * g(hi) < 0 else None

    # ---- neck map: bounce to next bounce (a_dot crosses zero upward)
    def flow_to_bounce(self, p, w, rtol=1e-13, method="DOP853", want_turn=False):
        a0 = self.neck_a(p, w)
        if a0 is None:
            return None
        y0 = [a0, 0.0, p, w]
        T1 = 0.05
        s1 = solve_ivp(self.f, (0.0, T1), y0, method=method, rtol=rtol, atol=1e-15)
        ev = lambda t, y: y[1]
        ev.direction = 1
        ev.terminal = True
        tr = lambda t, y: y[1]
        tr.direction = -1
        small = lambda t, y: y[0] - 1e-3
        small.terminal = True
        s2 = solve_ivp(self.f, (T1, 40.0), s1.y[:, -1], method=method, rtol=rtol, atol=1e-15, events=[ev, tr, small], max_step=0.05)
        if not s2.t_events[0].size:
            return None
        if want_turn:
            return float(s2.t_events[0][0]), s2.y_events[0][0], (s2.t_events[1][0] if s2.t_events[1].size else None, s2.y_events[1][0] if s2.t_events[1].size else None)
        return float(s2.t_events[0][0]), s2.y_events[0][0]

    def P(self, p, w, **k):
        r = self.flow_to_bounce(p, w, **k)
        return None if r is None else (r[1][2], r[1][3], r[0])

    def g_turn(self, p, rtol=1e-13):
        """phi_dot at the first turnaround for a start at the neck with phi_dot = 0 (time-symmetric orbit iff zero)."""
        a0 = self.neck_a(p, 0.0)
        if a0 is None:
            return None
        ev = lambda t, y: y[1]
        ev.direction = -1
        ev.terminal = True
        small = lambda t, y: y[0] - 1e-3
        small.terminal = True
        s1 = solve_ivp(self.f, (0.0, 0.05), [a0, 0.0, p, 0.0], method="DOP853", rtol=rtol, atol=1e-15)
        s2 = solve_ivp(self.f, (0.05, 40.0), s1.y[:, -1], method="DOP853", rtol=rtol, atol=1e-15, events=[ev, small], max_step=0.05)
        return float(s2.y_events[0][0][3]) if s2.t_events[0].size else None

    def orbit(self, guess_over_v, span=0.15):
        """Root of g_turn near guess (multiplicative bracket search)."""
        lo, hi = guess_over_v * (1 - span), guess_over_v * (1 + span)
        grid = np.linspace(lo, hi, 13)
        vals = [self.g_turn(z * self.v) for z in grid]
        for i in range(len(grid) - 1):
            if vals[i] is not None and vals[i + 1] is not None and vals[i] * vals[i + 1] < 0:
                return brentq(lambda z: self.g_turn(z * self.v), grid[i], grid[i + 1], xtol=1e-14) * self.v
        return None

    def jac_neck(self, p, w, h=1e-7):
        cols = []
        for dp, dw in ((h * self.v, 0.0), (0.0, h * self.v)):
            a, b = self.P(p + dp, w + dw), self.P(p - dp, w - dw)
            if a is None or b is None:
                return None
            cols.append([(a[0] - b[0]) / (2 * h * self.v), (a[1] - b[1]) / (2 * h * self.v)])
        return np.array(cols).T

    def monodromy(self, p, rtol=1e-13):
        """4-D monodromy of the variational equations over one period of the orbit through the neck state (phi = p, phi_dot = 0)."""
        a0 = self.neck_a(p, 0.0)
        T = self.flow_to_bounce(p, 0.0)[0]

        def aug(t, Y):
            y = Y[:4]
            Phi = Y[4:].reshape(4, 4)
            return np.concatenate([self.f(t, y), (self.jac(y) @ Phi).ravel()])
        Y0 = np.concatenate([[a0, 0.0, p, 0.0], np.eye(4).ravel()])
        sol = solve_ivp(aug, (0.0, T), Y0, method="DOP853", rtol=rtol, atol=1e-15)
        return T, sol.y[:4, -1], sol.y[4:, -1].reshape(4, 4)

    def conformal_orbit_check(self, p, rtol=1e-13):
        """Independent clock: integrate the same orbit in CONFORMAL time eta (dt = a d eta) and return (period in cosmic time, return error)."""
        a0 = self.neck_a(p, 0.0)

        def fe(eta, Y):
            a, b, ph, w, t = Y
            d = self.f(t, [a, b, ph, w])
            return [a * d[0], a * d[1], a * d[2], a * d[3], a]
        ev = lambda eta, Y: Y[1]
        ev.direction = 1
        ev.terminal = True
        s1 = solve_ivp(fe, (0.0, 0.03), [a0, 0.0, p, 0.0, 0.0], method="Radau", rtol=rtol, atol=1e-15)
        s2 = solve_ivp(fe, (0.03, 40.0), s1.y[:, -1], method="Radau", rtol=rtol, atol=1e-15, events=[ev], max_step=0.05)
        Y = s2.y_events[0][0]
        return float(Y[4]), np.array([Y[0] - a0, Y[1], Y[2] - p, Y[3]])


def branch(x_list, guess):
    out, g = [], guess
    for x in x_list:
        c = Cycle(x)
        p = c.orbit(g, span=0.12)
        if p is None:
            out.append(dict(x=float(x), found=False))
            continue
        J = c.jac_neck(p, 0.0)
        if J is None:
            out.append(dict(x=float(x), found=False)); continue
        g = p / c.v
        tr = float(np.trace(J))
        det = float(np.linalg.det(J))
        ev = np.linalg.eigvals(J)
        rot = None
        if abs(tr) < 2.0:
            rot = float(np.arccos(np.clip(tr / 2.0, -1, 1)) / (2 * np.pi))
        out.append(dict(x=float(x), found=True, phi0_over_v=float(g), period=float(c.flow_to_bounce(p, 0.0)[0]), trace=tr, det=det, multipliers=[[float(e.real), float(e.imag)] for e in ev], rotation_number=rot))
        print("  x=%.4f phi0/v=%.5f trace=%+.6f det=%.7f rot=%s" % (x, g, tr, det, "-" if rot is None else "%.5f" % rot), flush=True)
    return out


def period2_search(c, p_fix, J, delta_list=(0.003, 0.01, 0.03, 0.1)):
    """Fixed points of P^2 that are not fixed points of P, started along the eigenvector of the multiplier near -1."""
    w_, V_ = np.linalg.eig(J)
    k = int(np.argmin(np.abs(w_ + 1.0)))
    e = np.real(V_[:, k]); e = e / np.linalg.norm(e)
    sols = []

    def P2(x):
        r1 = c.P(p_fix + x[0] * c.v, x[1] * c.v)
        if r1 is None:
            return np.array([9.0, 9.0])
        r2 = c.P(r1[0], r1[1])
        if r2 is None:
            return np.array([9.0, 9.0])
        return np.array([(r2[0] - p_fix) / c.v - x[0], r2[1] / c.v - x[1]])
    for d in delta_list:
        for sgn in (1, -1):
            x0 = sgn * d * e
            try:
                xs, info, ier, msg = fsolve(P2, x0, full_output=True, xtol=1e-12)
            except Exception:           # noqa: BLE001
                continue
            if ier == 1 and np.linalg.norm(P2(xs)) < 1e-9 and np.linalg.norm(xs) > 2e-4:
                if not any(np.linalg.norm(xs - s) < 1e-6 for s in sols):
                    sols.append(xs)
    return [[float(s[0]), float(s[1])] for s in sols]


if __name__ == "__main__":
    out = {}
    # ---- 1. independent reproduction at C/Cc = 0.35, 0.355, 0.36
    out["reproduction"] = []
    for x, guess in ((0.35, 0.2229), (0.355, 0.2499), (0.36, 0.2818)):
        c = Cycle(x)
        p = c.orbit(guess, span=0.05)
        T, yT, M = c.monodromy(p)
        J = c.jac_neck(p, 0.0)
        Tc, errc = c.conformal_orbit_check(p)
        mono_ev = np.linalg.eigvals(M)
        rec = dict(x=x, phi0_over_v=float(p / c.v), period=float(T), return_err_full_state=float(np.linalg.norm(yT - np.array([c.neck_a(p, 0.0), 0.0, p, 0.0]))),
                   constraint_start=float(c.constraint([c.neck_a(p, 0.0), 0.0, p, 0.0])), constraint_end=float(c.constraint(yT)),
                   neck_map_trace=float(np.trace(J)), neck_map_det=float(np.linalg.det(J)), neck_map_mult=[[float(e.real), float(e.imag)] for e in np.linalg.eigvals(J)],
                   monodromy_eigs=sorted([[float(e.real), float(e.imag)] for e in mono_ev]), monodromy_det=float(np.linalg.det(M)),
                   conformal_period=Tc, conformal_return_err=float(np.linalg.norm(errc)))
        out["reproduction"].append(rec)
        print("x=%.3f phi0/v=%.8f T=%.8f (conformal-time clock %.8f) return err %.1e constraint %.1e/%.1e" % (x, rec["phi0_over_v"], T, Tc, rec["return_err_full_state"], rec["constraint_start"], rec["constraint_end"]))
        print("   neck map trace %.8f det %.8f | monodromy eigenvalues %s | det %.8f" % (rec["neck_map_trace"], rec["neck_map_det"], [complex(*e) for e in rec["monodromy_eigs"]], rec["monodromy_det"]), flush=True)
    # ---- 2. convergence in tolerance at C/Cc = 0.35
    out["convergence"] = []
    c = Cycle(0.35)
    for rt in (1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13):
        a0 = None
        # re-find the orbit with this tolerance
        gt = lambda z: c.g_turn(z * c.v, rtol=rt)
        z0 = brentq(gt, 0.2229 * 0.97, 0.2229 * 1.03, xtol=1e-14)
        p = z0 * c.v
        Pm = c.P(p, 0.0, rtol=rt)
        Jt = c.jac_neck(p, 0.0, h=1e-6)
        out["convergence"].append(dict(rtol=rt, phi0_over_v=float(z0), period=float(Pm[2]), return_err=float(np.hypot(Pm[0] - p, Pm[1]) / c.v)))
        print("  rtol=%g phi0/v=%.12f T=%.12f return err %.2e" % (rt, z0, Pm[2], out["convergence"][-1]["return_err"]), flush=True)
    for method in ("LSODA", "Radau", "DOP853"):
        p = c.orbit(0.2229, span=0.05)
        r = c.P(p, 0.0, rtol=1e-11, method=method)
        out["convergence"].append(dict(method=method, return_err=float(np.hypot(r[0] - p, r[1]) / c.v), period=float(r[2])))
        print("  method %-7s T=%.10f return err %.2e" % (method, r[2], out["convergence"][-1]["return_err"]), flush=True)
    # ---- 3. branch continuation across the bifurcation
    print("branch continuation")
    xs = list(np.round(np.arange(0.330, 0.3605, 0.0025), 5))
    out["branch"] = branch(xs, 0.1421)
    # ---- 4. flip location and period-2 orbit
    tr_of_x = {}

    def trace_at(x, guess):
        c_ = Cycle(x)
        p_ = c_.orbit(guess, span=0.2)
        return float(np.trace(c_.jac_neck(p_, 0.0))), p_ / c_.v
    g0 = [b for b in out["branch"] if b["found"] and abs(b["x"] - 0.345) < 1e-9][0]["phi0_over_v"]
    xl, xr = 0.345, 0.350
    f_ = lambda x: trace_at(x, g0 * (1 + (x - 0.345) * 4.7))[0] + 2.0
    xstar = brentq(f_, xl, xr, xtol=1e-10)
    out["flip"] = dict(C_over_Ccrit_star=float(xstar))
    print("flip (trace = -2) at C/Cc = %.8f" % xstar, flush=True)
    for dx in (-0.004, -0.0015, 0.0015, 0.004):
        x = xstar + dx
        tr, gz = trace_at(x, g0 * (1 + (x - 0.345) * 4.7))
        c = Cycle(x)
        p = gz * c.v
        J = c.jac_neck(p, 0.0)
        sols = period2_search(c, p, J)
        rec = dict(dx=dx, x=float(x), trace=float(np.trace(J)), period2_points=sols)
        out["flip"].setdefault("scan", []).append(rec)
        print("  C/Cc=%.5f (dx=%+.4f) trace %+.5f: period-2 points (offsets in v) %s" % (x, dx, np.trace(J), sols), flush=True)
    json.dump(out, open(os.path.join(HERE, "phase2d_cycle_validation.json"), "w"), indent=1)
