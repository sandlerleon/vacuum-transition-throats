# -*- coding: utf-8 -*-
"""Phase 2, Stage A: reduced leak-model feasibility test for a pre-existing bridge (Model B, energy flux through the bridge).  NOT a wormhole solver.

Question: can a bridge that carries energy out of the closed cold cycle (so the dissipative heat never accumulates as radiation) keep the exact cold cycles of closed_cycle.py bouncing,
and is a bridge needed, or would any radiating boundary do?

Interior: the closed FRW model of closed_cycle.py with dissipation G (phi_ddot + (3H + G) phi_dot + V' = 0), radiation rho_r (rho_r_dot + 4H rho_r = G phi_dot^2 - Q), and a leak Q (energy
density removed per unit time, rho_r only):
    'diffuse' : Q = eps rho_r                   (volume-proportional emission: any ordinary radiating boundary / transparent medium)
    'area'    : Q = eps rho_r / a               (flux through a fixed-angular-size mouth: shell area / volume ~ 1/a; the bridge-specific law)
    'sink'    : Q = G phi_dot^2, rho_r = 0      (ideal sink: ALL dissipated heat leaves instantly; the best any radiation leak could do)
Energy bookkeeping:  E_leak = int Q a^3 dt (per 2 pi^2) is the energy received by the exterior; m_leak = (pi/2) E_leak is the corresponding exterior mass parameter (8 pi G = 1).

Approximation, stated plainly: a homogeneous interior cannot lose energy and keep the Friedmann constraint a_dot^2 = a^2 rho/3 - 1 (at a_dot = 0 the density is pinned by a).  We keep the
acceleration equation and the scalar/radiation equations as written and TRACK the effective curvature K_eff = a^2 rho/3 - a_dot^2 (K_eff = 1 initially); d K_eff/dt = -a^2 Q/3 exactly.  The model is
trusted only while |K_eff - 1| << 1; a drift means the leak has changed the interior's own global constraint, i.e. the homogeneous reduced model is no longer self-consistent and the
inhomogeneous (Stage B) problem would have to be solved.  Shell stresses are evaluated quasi-statically with the growing exterior mass (bridge_junction.shell).

    python bridge_leak.py          (writes bridge_leak.json)
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp

import closed_cycle as CC
import bridge_junction as BJ

HERE = os.path.dirname(os.path.abspath(__file__))
VF = CC.VF
CHI0 = 0.7 * np.pi          # wall position used for the quasi-static shell stresses
QEXT = 1.0


class LeakModel(CC.Model):
    def __init__(self, q, rc, mu, C, gamma=0.0, eps=0.0, law="none"):
        super().__init__(q=q, rc=rc, mu=mu, C=C, ups=gamma)
        self.gamma, self.eps, self.law = gamma, eps, law

    def Qleak(self, a, pd, rr):
        if self.law == "diffuse":
            return self.eps * rr
        if self.law == "area":
            return self.eps * rr / a
        if self.law == "sink":
            return self.gamma * pd * pd
        return 0.0

    def rhs(self, t, y):
        a, ad, p, pd, rr, El = y
        d = super().rhs(t, [a, ad, p, pd, rr])
        Q = self.Qleak(a, pd, rr)
        if self.law == "sink":
            d[4] = 0.0                               # radiation never accumulates; the heat leaves at once
        else:
            d[4] -= Q
        return d + [Q * a ** 3]

    def keff(self, y):
        a, ad, p, pd, rr = y[:5]
        return a * a * self.rho(a, p, pd, rr) / 3.0 - ad * ad


def run_leak(m, y0, t_max, rtol=1e-11):
    ev = lambda t, y: y[1]
    ev.direction = 0
    small = lambda t, y: y[0] - 1e-3
    small.terminal = True
    big = lambda t, y: y[0] - 1e6
    big.terminal = True
    sol = solve_ivp(m.rhs, (0.0, t_max), y0, method="DOP853", rtol=rtol, atol=1e-15, events=[ev, small, big], max_step=0.02, dense_output=True)
    B, T = [], []
    for k, tt in enumerate(sol.t_events[0]):
        if tt < 1e-6:
            continue
        ye = sol.y_events[0][k]
        (B if m.rhs(tt, ye)[1] > 0 else T).append((float(tt), ye))
    end = "singular" if sol.t_events[1].size else ("runaway" if sol.t_events[2].size else "t_max")
    return sol, B, T, end


def start(m, seed):
    a0 = CC.neck_a(m, seed, 0.0, 0.0)
    return [a0, 0.0, seed, 0.0, 0.0, 0.0]


def shell_stress(m, sol, m_max_ref):
    """Quasi-static shell stresses at chi0 along the trajectory, exterior mass m(t) = (pi/2) E_leak(t); returns fractions of the time with sigma > 0 and with sigma + p >= 0."""
    t = np.linspace(0.0, sol.t[-1], 800)
    Y = sol.sol(t)
    n_ok = s_ok = valid = 0
    for i in range(len(t)):
        y = Y[:, i]
        add = m.rhs(t[i], y)[1]
        mm = 0.5 * np.pi * y[5] + 1e-6
        res = BJ.shell(mm, QEXT, CHI0, y[0], y[1], add)
        if res is None:
            continue
        valid += 1
        sg, p, R = res
        s_ok += bool(sg > 0)
        n_ok += bool(sg > 0 and sg + p >= 0)
    return dict(frac_valid=valid / len(t), frac_sigma_pos=s_ok / len(t), frac_nec=n_ok / len(t))


def lyapunov(m, seed, t_max, d=1e-9):
    """Finite-time exponent: growth of the separation of two nearby cold-cycle starts (state a, a_dot, phi, phi_dot) until it reaches 1e-3 or the run ends."""
    ys = [start(m, seed * (1 + dd)) for dd in (0.0, d)]
    sols = [run_leak(m, y, t_max)[0] for y in ys]
    tend = min(s.t[-1] for s in sols)
    t = np.linspace(0.0, tend, 4000)
    sep = np.linalg.norm(sols[0].sol(t)[:4] - sols[1].sol(t)[:4], axis=0)
    s0 = sep[0]
    hit = np.where(sep > 1e-3)[0]
    if hit.size:
        return float(np.log(sep[hit[0]] / s0) / t[hit[0]]), True
    return float(np.log(max(sep[-1], 1e-300) / s0) / max(tend, 1e-9)), False


def one_case(o, gamma, law, eps, orbit_label):
    m = LeakModel(q=o["q"], rc=o["rc"], mu=o["mu"], C=o["C_over_Ccrit"] * CC.C_CRIT, gamma=gamma, eps=eps, law=law)
    seed = 10 ** o["log10_seed_over_v"] * m.v
    T = o["period"]
    t_max = 40.0 * T
    y0 = start(m, seed)
    sol, B, T_, end = run_leak(m, y0, t_max)
    ks = [m.keff(b[1]) for b in B]
    Yall = sol.y
    kall = np.array([m.keff(Yall[:, i]) for i in range(0, Yall.shape[1], max(1, Yall.shape[1] // 400))])
    E_int = [m.rho(b[1][0], b[1][2], b[1][3], b[1][4]) * b[1][0] ** 3 for b in B]
    rec = dict(orbit=orbit_label, gamma=gamma, law=law, eps=eps, bounces=len(B), turnarounds=len(T_), end=end, t_end=float(sol.t[-1]), period=T,
               a_neck=[float(b[1][0]) for b in B], radiation_over_Vf=[float(b[1][4] / VF) for b in B], keff_at_necks=[float(k) for k in ks],
               keff_min=float(kall.min()), keff_max=float(kall.max()), E_leak=float(sol.y[5][-1]), m_leak=float(0.5 * np.pi * sol.y[5][-1]),
               E_int_first_neck=float(E_int[0]) if E_int else None)
    rec["shell"] = shell_stress(m, sol, None) if sol.t[-1] > 0 else None
    if len(B) >= 1 and law in ("none", "diffuse", "area", "sink") and gamma >= 0:
        try:
            lam, hit = lyapunov(m, seed, min(t_max, 12.0 * T))
            rec["lyapunov"], rec["lyapunov_reached"] = lam, hit
        except Exception as e:                     # noqa: BLE001
            rec["lyapunov"], rec["lyapunov_reached"] = None, False
    return rec


if __name__ == "__main__":
    breathers = json.load(open(os.path.join(HERE, "closed_stability_breather.json")))
    pick = [("mu20_rc0.2_C0.3", breathers[3]), ("mu10_rc0.4_C0.6", breathers[2])]
    out = []
    for label, o in pick:
        for gamma in (0.0, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0):
            cases = [("none", 0.0)]
            if gamma > 0:
                cases += [("diffuse", e) for e in (0.1, 1.0, 10.0, 100.0)] + [("area", e) for e in (0.1, 1.0, 10.0, 100.0)] + [("sink", 0.0)]
            for law, eps in cases:
                r = one_case(o, gamma, law, eps, label)
                out.append(r)
                print("%-16s G=%-5g %-8s eps=%-6g | bounces %2d turn %2d %-8s t=%5.1f | Keff[%.3f,%.3f] | m_leak %.2e | lam %s | shell sigma>0 %.2f NEC %.2f" %
                      (label, gamma, law, eps, r["bounces"], r["turnarounds"], r["end"], r["t_end"], r["keff_min"], r["keff_max"], r["m_leak"],
                       "-" if r.get("lyapunov") is None else "%.1f%s" % (r["lyapunov"], "" if r["lyapunov_reached"] else "*"), r["shell"]["frac_sigma_pos"], r["shell"]["frac_nec"]), flush=True)
    json.dump(out, open(os.path.join(HERE, "bridge_leak.json"), "w"), indent=1)
