# -*- coding: utf-8 -*-
"""Phase 3H.0: does causal shear viscosity relax the anisotropy obstruction of the spin-torsion bounce?  Combines Phase 3C.2 (axisymmetric Bianchi I, spin fluid rho_s = p_s = -B n^2, critical spin strength B_c set by the shear)
with Phase 3F.0 (Israel-Stewart-type shear stress, eta = c^2 tau n, tau = K / max(|Theta|, Theta_floor)).  kappa = 1.
Model (stated approximations): dust number density n = exp(-alpha-2beta); dust energy density rho_d evolves with viscous heating, rho_d' = -Theta rho_d - (alpha'-beta') S, heat treated as pressureless;
shear stress S = pi_xx (pi_yy = pi_zz = -S/2): tau S' + S = -(4/3) eta (alpha'-beta'); total rho = rho_d - B n^2, p = -B n^2, p_x = p + S, p_y = p - S/2.
Outputs: B_c (smallest B giving a bounce, Theta crossing zero upward before V < 1e-3) for (K, c^2) with the causal bound c^2 <= 1 and the hydrodynamic bound K <= 1, the isotropic threshold, and the Knudsen number actually reached.
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
THF = 1e-2


def run(B, K, c2, s0=0.05, tmax=6000.0, Vmin=1e-9):
    rho0 = 1.0 - B
    h = -np.sqrt(rho0 / (3 * (1 - s0 ** 2)))
    y0 = [0.0, 0.0, h * (1 + 2 * s0), h * (1 - s0), 0.0, 1.0]

    def rhs(t, y):
        al, be, ad, bd, S, rd = y
        n = np.exp(-al - 2 * be)
        Th = ad + 2 * bd
        tau = K / max(abs(Th), THF) if K > 0 else 1.0
        eta = c2 * tau * n if K > 0 else 0.0
        p = -B * n * n
        rho = rd - B * n * n
        bdd = 0.5 * (-(p + S) - 3 * bd * bd)
        add = -(p - 0.5 * S) - bdd - ad * ad - bd * bd - ad * bd
        Sd = (-S - (4.0 / 3.0) * eta * (ad - bd)) / tau if K > 0 else 0.0
        rdd = -Th * rd - (ad - bd) * S
        return [ad, bd, add, bdd, Sd, rdd]

    bounce = lambda t, y: (y[2] + 2 * y[3]); bounce.direction = 1; bounce.terminal = True
    small = lambda t, y: y[0] + 2 * y[1] - np.log(Vmin); small.terminal = True
    sol = solve_ivp(rhs, (0, tmax), y0, method="DOP853", rtol=1e-11, atol=1e-14, events=[bounce, small], max_step=0.5)
    al, be, ad, bd, S, rd = sol.y
    n = np.exp(-al - 2 * be)
    rho = rd - B * n * n
    resid = float(np.max(np.abs(2 * ad * bd + bd ** 2 - rho)))
    Th = ad + 2 * bd
    kn = float(np.max(np.minimum(np.abs(Th), 1e9) * (K / np.maximum(np.abs(Th), THF)))) if K > 0 else 0.0
    return dict(bounce=bool(len(sol.t_events[0])), V_end=float(np.exp(al[-1] + 2 * be[-1])), residual=resid, knudsen_max=kn,
                Sigma2_end=float(((ad[-1] - bd[-1]) ** 2) / max((ad[-1] + 2 * bd[-1]) ** 2, 1e-30)), S_over_rho_max=float(np.max(np.abs(S) / np.maximum(np.abs(rho), 1e-12))), S_over_rhod_max=float(np.max(np.abs(S) / np.maximum(rd, 1e-12))), S_over_n_max=float(np.max(np.abs(S) / n)))


def critical(K, c2, s0=0.05, lo=1e-8, hi=0.5):
    if not run(hi, K, c2, s0)["bounce"]:
        return None
    if run(lo, K, c2, s0)["bounce"]:
        return lo
    for _ in range(44):
        mid = np.sqrt(lo * hi)
        if run(mid, K, c2, s0)["bounce"]: hi = mid
        else: lo = mid
    return float(hi)


out = {}
rows = []
for s0 in (0.05, 0.2):
    base = critical(0.0, 0.0, s0)
    iso = critical(0.0, 0.0, 0.0)
    print("H: s0 = %.2f: inviscid B_c = %.5g (isotropic reference %.3g)" % (s0, base, iso))
    rows.append(dict(s0=s0, K=0.0, c2=0.0, B_c=base, ratio=1.0, isotropic_ref=iso))
    for K, c2 in ((0.1, 1.0), (0.3, 1.0), (1.0, 0.25), (1.0, 1.0), (3.0, 1.0)):
        Bc = critical(K, c2, s0)
        r_ = run(Bc * 1.5, K, c2, s0) if Bc else None
        rows.append(dict(s0=s0, K=K, c2=c2, B_c=Bc, ratio=(Bc / base if Bc else None), knudsen_max_at_1p5Bc=None if r_ is None else r_["knudsen_max"], residual=None if r_ is None else r_["residual"],
                         S_over_rho_max=None if r_ is None else r_["S_over_rho_max"], S_over_rhod_max=None if r_ is None else r_["S_over_rhod_max"], S_over_n_max=None if r_ is None else r_["S_over_n_max"]))
        print("H:   K = %.1f, c^2 = %.2f: B_c = %s (x%.3f of inviscid) ; Knudsen reached %s ; constraint residual %s ; max |S|/rho_dust %s, |S|/n %s" %
              (K, c2, ("%.5g" % Bc) if Bc else None, (Bc / base) if Bc else float("nan"), ("%.2f" % r_["knudsen_max"]) if r_ else "-", ("%.1e" % r_["residual"]) if r_ else "-", ("%.2f" % r_["S_over_rhod_max"]) if r_ else "-", ("%.2f" % r_["S_over_n_max"]) if r_ else "-"))
out["rows"] = rows
out["approximations"] = "dust heated by viscous work treated as pressureless; viscosity acts on the dust only; truncated Israel-Stewart; Theta_floor = 1e-2 caps tau near the turning point"
json.dump(out, open(os.path.join(HERE, "phase3h0_viscous_spin_bounce.json"), "w"), indent=1)
