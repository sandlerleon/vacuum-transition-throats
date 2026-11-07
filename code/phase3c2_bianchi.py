# -*- coding: utf-8 -*-
"""Phase 3C.2/3C.3: exact axisymmetric Bianchi-I minisuperspace with an Einstein-Cartan spin fluid, a rotation proxy and a uniform magnetic field.  Units 8 pi G = 1.

Metric  ds^2 = -dt^2 + e^{2 alpha} dx^2 + e^{2 beta} (dy^2 + dz^2);  H_x = alpha_dot, H_y = H_z = beta_dot, theta = alpha_dot + 2 beta_dot, V = e^{alpha + 2 beta}, a_perp = e^{beta}.
Einstein equations G_{mu nu} = T_{mu nu}:
    constraint        2 alpha_dot beta_dot + beta_dot^2 = rho                      (= theta^2/3 - sigma^2)
    x equation        2 beta_ddot + 3 beta_dot^2 = - p_x
    y equation        alpha_ddot + beta_ddot + alpha_dot^2 + beta_dot^2 + alpha_dot beta_dot = - p_y
Matter (all conserved separately; checked numerically):
  dust            n = N0 / V,  rho_d = n, p = 0
  spin fluid      random spin orientation (Pasmatsiou, Tsagas, Barrow 2016: <s_ab> = 0, <s^2> != 0): rho_s = - B n^2, p_s = - B n^2 (isotropic), B = kappa sigma^2 / 4
  rotation proxy  azimuthally averaged orbital motion in the equatorial plane with conserved specific angular momentum (circulation law omega ~ a_perp^{-2}, v ~ 1/a_perp, from the vorticity equation):
                  rho_r = (1/2) n v^2, v^2 = v0^2 e^{-2 beta};  p_x = 0, p_y = p_z = rho_r  (kinetic stress of orbiting dust).  Reproduces the centrifugal stress and its anisotropy (flattening),
                  NOT frame dragging or a true global vorticity: a controlled slow-rotation (O(v^2)) proxy.
  magnetic field  along x, flux conserved: B_x ~ e^{-2 beta}, rho_B = b0^2 e^{-4 beta} / 2;  p_x = - rho_B, p_y = p_z = + rho_B.
Outcomes: 'bounce' (theta crosses zero while contracting), 'singular' (V < V_min) .  Critical spin strength B_c = smallest B that gives a bounce (bisection in log B).

    python phase3c2_bianchi.py          (writes phase3c2_bianchi.json)
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
VMIN = 1e-9


class Model:
    def __init__(self, B=0.0, v0sq=0.0, rhoB0=0.0, N0=1.0):
        self.B, self.v0sq, self.b0sq, self.N0 = B, v0sq, 2.0 * rhoB0, N0

    def parts(self, al, be):
        V = np.exp(al + 2 * be)
        n = self.N0 / V
        rho_d = n
        rho_s = -self.B * n * n; p_s = rho_s
        rho_r = 0.5 * n * self.v0sq * np.exp(-2 * be)
        rho_B = 0.5 * self.b0sq * np.exp(-4 * be)
        rho = rho_d + rho_s + rho_r + rho_B
        px = p_s - rho_B
        py = p_s + rho_r + rho_B
        return rho, px, py, dict(dust=rho_d, spin=rho_s, rot=rho_r, mag=rho_B)

    def rhs(self, t, y):
        al, be, ad, bd = y
        rho, px, py, _ = self.parts(al, be)
        bdd = 0.5 * (-px - 3 * bd * bd)
        add = -py - bdd - ad * ad - bd * bd - ad * bd
        return [ad, bd, add, bdd]

    def constraint(self, y):
        al, be, ad, bd = y
        return 2 * ad * bd + bd * bd - self.parts(al, be)[0]


def initial(model, shear):
    """contracting start at V = 1 with shear parameter s: alpha_dot = h (1 + 2 s), beta_dot = h (1 - s); constraint fixes h."""
    rho = model.parts(0.0, 0.0)[0]
    h = -np.sqrt(rho / (3.0 * (1.0 - shear ** 2)))
    return [0.0, 0.0, h * (1 + 2 * shear), h * (1 - shear)]


def run(model, shear=0.0, tmax=2000.0):
    y0 = initial(model, shear)
    bounce = lambda t, y: y[2] + 2 * y[3]
    bounce.direction = 1
    bounce.terminal = True
    small = lambda t, y: y[0] + 2 * y[1] - np.log(VMIN)
    small.terminal = True
    sol = solve_ivp(model.rhs, (0, tmax), y0, method="DOP853", rtol=1e-12, atol=1e-14, events=[bounce, small], dense_output=True)
    if sol.t_events[0].size:
        ye = sol.y_events[0][0]
        Vb = float(np.exp(ye[0] + 2 * ye[1]))
        return dict(outcome="bounce", V_min=Vb, t=float(sol.t_events[0][0]), constraint=float(abs(model.constraint(ye))), sol=sol)
    if sol.t_events[1].size:
        return dict(outcome="singular", V_min=VMIN, t=float(sol.t_events[1][0]), constraint=float(abs(model.constraint(sol.y_events[1][0]))), sol=sol)
    return dict(outcome="none", V_min=None, t=float(sol.t[-1]), constraint=None, sol=sol)


def critical_B(v0sq=0.0, rhoB0=0.0, shear=0.0, lo=1e-8, hi=0.5):
    f = lambda B: run(Model(B=B, v0sq=v0sq, rhoB0=rhoB0), shear)["outcome"] == "bounce"
    if not f(hi):
        return None
    if f(lo):
        return lo
    a, b = np.log(lo), np.log(hi)
    for _ in range(40):
        m = 0.5 * (a + b)
        if f(np.exp(m)): b = m
        else: a = m
    return float(np.exp(b))


if __name__ == "__main__":
    out = {}
    # ---------------------------------------------------------------- checks
    chk = {}
    m = Model(B=0.05, v0sq=0.04, rhoB0=0.03)
    r = run(m, shear=0.02)
    sol = r["sol"]
    tt = np.linspace(0, r["t"] * 0.999, 4000); Y = sol.sol(tt)
    chk["constraint_max"] = float(max(abs(m.constraint(Y[:, i])) / max(abs(m.parts(Y[0, i], Y[1, i])[0]), 1e-12) for i in range(0, len(tt), 10)))
    # separate conservation: rho_dot + sum_i H_i (rho + p_i) = 0 along the solution for the total stress
    rho_t = np.array([m.parts(Y[0, i], Y[1, i])[0] for i in range(len(tt))]); px = np.array([m.parts(Y[0, i], Y[1, i])[1] for i in range(len(tt))]); py = np.array([m.parts(Y[0, i], Y[1, i])[2] for i in range(len(tt))])
    rd = np.gradient(rho_t, tt); cons = rd + Y[2] * (rho_t + px) + 2 * Y[3] * (rho_t + py)
    chk["conservation_rel"] = float(np.max(np.abs(cons[5:-5]) / np.max(np.abs(rd[5:-5]))))
    # GR limit: dust only, small shear: sigma^2 V^2 = const (Kasner-like) and no bounce
    g = run(Model(), shear=0.1)
    ts = np.linspace(0, g["t"] * 0.99, 200); Yg = g["sol"].sol(ts)
    Vg = np.exp(Yg[0] + 2 * Yg[1]); th = Yg[2] + 2 * Yg[3]; sg2 = (th ** 2 / 3 - np.array([Model().parts(Yg[0, i], Yg[1, i])[0] for i in range(len(ts))]))
    chk["GR_limit_outcome"] = g["outcome"]; chk["sigma2_V2_spread"] = float(np.ptp(sg2 * Vg ** 2) / np.mean(sg2 * Vg ** 2))
    # analytic bounce volume: isotropic-asymptotics V_b = B N0 - Sigma  (N0 = 1) for dust + spin + constant shear Sigma = sigma^2 V^2
    s0 = 0.05; Bt = 0.02
    mm = Model(B=Bt); rr = run(mm, shear=s0)
    y0 = initial(mm, s0); Vs = 1.0; th0 = y0[2] + 2 * y0[3]; rho0 = mm.parts(0, 0)[0]; Sig = (th0 ** 2 / 3 - rho0) * Vs ** 2
    chk["analytic_Vb"] = float(Bt - Sig); chk["numeric_Vb"] = rr["V_min"]; chk["Sigma_initial"] = float(Sig)
    out["checks"] = chk
    print("checks:", {k: (v if isinstance(v, str) else "%.3e" % v) for k, v in chk.items()})
    # ---------------------------------------------------------------- critical spin strength for the models
    rows = []
    cases = [("A  nonrotating EC (isotropic start)", 0.0, 0.0, 0.0), ("A' nonrotating EC, initial shear 0.05", 0.0, 0.0, 0.05),
             ("C  rotating EC v0^2 = 0.01", 0.01, 0.0, 0.0), ("C  rotating EC v0^2 = 0.05", 0.05, 0.0, 0.0), ("C  rotating EC v0^2 = 0.2", 0.2, 0.0, 0.0),
             ("D  magnetized nonrotating, rho_B0 = 0.01", 0.0, 0.01, 0.0), ("D  magnetized nonrotating, rho_B0 = 0.05", 0.0, 0.05, 0.0), ("D  magnetized nonrotating, rho_B0 = 0.2", 0.0, 0.2, 0.0),
             ("E  rotating + magnetized, 0.05 / 0.05", 0.05, 0.05, 0.0), ("E  rotating + magnetized, 0.2 / 0.2", 0.2, 0.2, 0.0)]
    for name, v0, rb, sh in cases:
        Bc = critical_B(v0sq=v0, rhoB0=rb, shear=sh)
        # properties at a supercritical B (10 x B_c): minimum volume, density at the bounce
        sup = None
        if Bc is not None:
            m2 = Model(B=10 * Bc if Bc > 1e-8 else 1e-3, v0sq=v0, rhoB0=rb)
            r2 = run(m2, sh)
            if r2["outcome"] == "bounce":
                ye = r2["sol"].y_events[0][0]
                tt2 = np.linspace(0, r2["t"], 2000); Y2 = r2["sol"].sol(tt2)
                th2 = Y2[2] + 2 * Y2[3]; rho2 = np.array([m2.parts(Y2[0, i], Y2[1, i])[0] for i in range(len(tt2))])
                sig2 = th2 ** 2 / 3 - rho2
                sup = dict(B=float(m2.B), V_bounce=r2["V_min"], max_shear_sq=float(sig2.max()), constraint=r2["constraint"])
        rows.append(dict(case=name, v0sq=v0, rhoB0=rb, shear0=sh, B_critical=Bc, supercritical=sup))
        print("%-46s B_c = %s %s" % (name, "none (<= 0.5)" if Bc is None else "%.4g" % Bc, "" if sup is None else "| at 10 B_c: V_bounce %.3g, max sigma^2 %.3g" % (sup["V_bounce"], sup["max_shear_sq"])), flush=True)
    out["critical"] = rows
    # model B: rotating GR without torsion (B = 0): always singular?
    mb = [dict(v0sq=v0, rhoB0=rb, outcome=run(Model(B=0.0, v0sq=v0, rhoB0=rb), 0.0)["outcome"]) for v0, rb in ((0.0, 0.0), (0.05, 0.0), (0.2, 0.0), (0.0, 0.2), (0.2, 0.2))]
    out["model_B_no_torsion"] = mb
    print("rotating/magnetized GR without torsion (B = 0):", [x["outcome"] for x in mb])
    json.dump(out, open(os.path.join(HERE, "phase3c2_bianchi.json"), "w"), indent=1)
