# -*- coding: utf-8 -*-
"""Phase 3C (preliminary feasibility, homogeneous): effective fluid of a spin (Einstein-Cartan, Weyssenhoff) dust ball, and the matching to a vacuum exterior.  Units 8 pi G = 1 (kappa = 1).

Effective fluid (literature form for a Weyssenhoff fluid, taken as given here, NOT re-derived):   rho_eff = n - B n^2 ,   p_eff = -B n^2 ,   B = kappa sigma^2 / 4 >= 0   (s = sigma n),
with n = 3C/a^3 the conserved number (rest-mass) density.  Internal consistency (verified below): thermodynamic identity p = n d rho/dn - rho, continuity d(rho_eff)/dt + 3H(rho_eff + p_eff) = 0, Bianchi.
The loop-quantum-cosmology-inspired modified Friedmann equation H^2 = (kappa/3) rho (1 - rho/rho_c) for dust is, as a conserved effective fluid, the SAME fluid with B = 1/rho_c (shown below).

Closed ball (K = +1):   a_dot^2 = (a^2/3)(n - B n^2) - 1 = C/a - 3 B C^2/a^4 - 1.      Bounce / turnaround: f(a) = C a^3 - a^4 - 3 B C^2 = 0 (two positive roots iff B <= 9 C^2/256).
Ball of coordinate angle chi0 with vacuum Schwarzschild exterior of mass parameter m_ext = (rho_dust R^3 / 3) = C sin^3(chi0) (F = 1 + R_dot^2 - m/R, horizon R = m).

    python phase3c_ec_homogeneous.py          (writes phase3c_ec_homogeneous.json)
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}


def roots(C, B):
    f = lambda a: C * a ** 3 - a ** 4 - 3 * B * C * C
    amid = 0.75 * C
    if f(amid) <= 0:
        return None
    return brentq(f, 1e-9 * C, amid), brentq(f, amid, C)       # (a_bounce, a_turnaround)


# ---------------------------------------------------------------- (a) the fluid and its consistency
B_ = 3.0e-7; C_ = 0.01
nn = np.linspace(1.0, 2.0 / B_ * 0.9, 20001)
rho = nn - B_ * nn ** 2; p = -B_ * nn ** 2
thermo = float(np.max(np.abs(nn * np.gradient(rho, nn) - rho - p)[5:-5] / np.abs(p[5:-5])))
cs2 = np.gradient(p, nn) / np.gradient(rho, nn)
n_star = 1.0 / (2 * B_)
out["fluid"] = dict(B=B_, n_star=n_star, rho_max=float(1.0 / (4 * B_)), thermodynamic_identity_residual=thermo,
                    cs2_formula="c_s^2 = -2 B n / (1 - 2 B n)", cs2_below_nstar_max=float(np.max(cs2[nn < 0.98 * n_star])), cs2_below_nstar_min=float(np.min(cs2[nn < 0.98 * n_star])),
                    cs2_above_nstar_min=float(np.min(cs2[nn > 1.02 * n_star])), NEC_violated_above_nstar=bool(np.all((rho + p)[nn > 1.001 * n_star] < 0)), NEC_holds_below_nstar=bool(np.all((rho + p)[nn < 0.999 * n_star] > 0)))
print("fluid: thermodynamic identity residual %.1e; c_s^2 in [%.3g, %.3g] for n < n*, min c_s^2 above n* = %.3f (superluminal if > 1); NEC violated only above n* (%s)" %
      (thermo, out["fluid"]["cs2_below_nstar_min"], out["fluid"]["cs2_below_nstar_max"], out["fluid"]["cs2_above_nstar_min"], out["fluid"]["NEC_violated_above_nstar"]))

# LQC-inspired equivalence: rho_eff = rho (1 - rho/rho_c); rho_dot = -3 H rho  =>  p_eff = -rho_eff - rho_eff_dot/(3H)
rc = 1.0 / B_
rho_m = np.linspace(1.0, rc * 0.99, 4001)
rho_eff = rho_m * (1 - rho_m / rc)
p_eff = -rho_eff + rho_m * (1 - 2 * rho_m / rc)                    # -rho_eff_dot/(3H) with rho_dot = -3 H rho
out["lqc_equivalence_max_diff"] = float(np.max(np.abs(p_eff - (-rho_m ** 2 / rc))))
print("LQC-type effective fluid p_eff = -rho^2/rho_c: max deviation %.1e (identical to the spin fluid with B = 1/rho_c)" % out["lqc_equivalence_max_diff"])

# ---------------------------------------------------------------- (b) homogeneous closed ball: bounce and cycle
def cycle(C, B):
    ab, at = roots(C, B)
    f = lambda t, y: [y[1], -0.5 * C / y[0] ** 2 + 6 * B * C * C / y[0] ** 5]         # from a_dot^2 = C/a - 3BC^2/a^4 - 1 :  a'' = -C/(2a^2) + 6 B C^2 / a^5
    ev = lambda t, y: y[1]
    ev.direction = 1
    ev.terminal = True
    s1 = solve_ivp(f, (0, 1e-3 * C), [ab, 0.0], method="DOP853", rtol=1e-13, atol=1e-16)
    s2 = solve_ivp(f, (1e-3 * C, 20 * np.pi * C), s1.y[:, -1], method="DOP853", rtol=1e-13, atol=1e-16, events=[ev], dense_output=True)
    T = float(s2.t_events[0][0])
    tt = np.linspace(0, T, 4000)
    Y = np.concatenate([s1.y[:, -1:], s2.sol(np.linspace(1e-3 * C, T, 3999))], axis=1) if False else None
    sol = solve_ivp(f, (0, T), [ab, 0.0], method="DOP853", rtol=1e-13, atol=1e-16, dense_output=True)
    Y = sol.sol(tt)
    return ab, at, T, tt, Y


C, B = 0.01, 2.0e-7
ab, at, T, tt, Y = cycle(C, B)
a, ad = Y[0], Y[1]
fried = np.abs(ad ** 2 - (C / a - 3 * B * C * C / a ** 4 - 1))
n = 3 * C / a ** 3; rho_e = n - B * n ** 2; p_e = -B * n ** 2
H = ad / a
rho_dot = np.gradient(rho_e, tt)
cont = np.abs(rho_dot[5:-5] + 3 * H[5:-5] * (rho_e[5:-5] + p_e[5:-5])) / np.max(np.abs(rho_dot[5:-5]))
out["cycle"] = dict(C=C, B=B, a_bounce=float(ab), a_turnaround=float(at), period=T, return_error=float(abs(Y[0, -1] - ab) / ab), friedmann_residual=float(fried.max()), continuity_residual=float(cont.max()),
                    n_at_bounce=float(3 * C / ab ** 3), n_star=float(1 / (2 * B)), bounce_above_nstar=bool(3 * C / ab ** 3 > 1 / (2 * B)), max_density_ratio=float(rho_e.max() / (1 / (4 * B))))
print("cycle (C=%g, B=%g): bounce a=%.5f, turnaround a=%.5f, period %.5f, return error %.1e, Friedmann residual %.1e, continuity residual %.1e; n at bounce %.3g vs n* = %.3g" %
      (C, B, ab, at, T, out["cycle"]["return_error"], fried.max(), cont.max(), 3 * C / ab ** 3, 1 / (2 * B)))

# ---------------------------------------------------------------- (c) matching to a vacuum exterior
rows = []
for chi0 in (0.4, 0.7, 1.0, 1.3, 1.5):
    sn = np.sin(chi0); m_ext = C * sn ** 3
    for Bv in (1e-9, 1e-8, 1e-7, 2e-7, 5e-7, 1e-6, 3e-6):
        r = roots(C, Bv)
        if r is None:
            continue
        ab_, at_, T_, tt_, Y_ = cycle(C, Bv)
        Rs = Y_[0] * sn
        mfloor = float(Rs.min() / m_ext)                              # minimum boundary radius over the cycle in units of the horizon radius m
        peak = float((1 + Y_[1] ** 2).max() * sn ** 2)                  # boundary trapped (rho_eff a^2 sin^2 chi0 / 3 >= 1) iff (1 + a_dot^2) sin^2 chi0 >= 1
        m_eff = (Y_[0] * sn) ** 3 * ((3 * C / Y_[0] ** 3) - Bv * (3 * C / Y_[0] ** 3) ** 2) / 3
        rows.append(dict(chi0=chi0, B=Bv, m_ext=float(m_ext), R_bounce_over_m=mfloor, trapped_at_boundary=bool(peak >= 1.0), m_eff_min_over_m=float(m_eff.min() / m_ext), m_eff_max_over_m=float(m_eff.max() / m_ext),
                         p_over_rho_at_bounce=float(-Bv * (3 * C / ab_ ** 3) ** 2 / ((3 * C / ab_ ** 3) - Bv * (3 * C / ab_ ** 3) ** 2))))
cls = lambda r: ("no_horizon_bounce_outside" if (r["R_bounce_over_m"] > 1.0) else ("horizon_forms_bounce_inside_exterior_horizon" if r["trapped_at_boundary"] else "other"))
out["matching"] = dict(rows=rows, counts={k: sum(1 for r in rows if cls(r) == k) for k in ("no_horizon_bounce_outside", "horizon_forms_bounce_inside_exterior_horizon", "other")})
print("matching classification over %d (chi0, B) cases: %s" % (len(rows), out["matching"]["counts"]))
for r in rows[::4]:
    print("  chi0=%.1f B=%.0e: R_bounce/m = %.3g, trapped %s, m_eff/m in [%.3g, %.3g], p/rho at bounce %.3g -> %s" % (r["chi0"], r["B"], r["R_bounce_over_m"], r["trapped_at_boundary"], r["m_eff_min_over_m"], r["m_eff_max_over_m"], r["p_over_rho_at_bounce"], cls(r)))
# the theorem: a shell cannot bounce inside the horizon of a static vacuum exterior:  Gamma_out^2 = F = 1 + R_dot^2 - m/R >= 0  at R_dot = 0  =>  R >= m
out["theorem_check"] = dict(any_case_with_bounce_inside_horizon_and_F_nonneg=bool(any(r["R_bounce_over_m"] < 1.0 for r in rows)), n_cases_bounce_inside=int(sum(1 for r in rows if r["R_bounce_over_m"] < 1.0)))
print("cases whose boundary bounces at R < m (F < 0 at R_dot = 0, impossible for a vacuum exterior): %d of %d" % (out["theorem_check"]["n_cases_bounce_inside"], len(rows)))
json.dump(out, open(os.path.join(HERE, "phase3c_ec_homogeneous.json"), "w"), indent=1)
