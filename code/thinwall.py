# -*- coding: utf-8 -*-
"""Exact thin-wall (Israel) analysis of a false-vacuum region in a true-vacuum exterior, and the cosmological mass and entropy bookkeeping.

Units: 8 pi G = 1 (M_Pl = 1). Interior: de Sitter with density rho_in (R_dS = sqrt(3 / rho_in)); exterior: Schwarzschild-(A)dS of mass parameter m = 2 G M and density
rho_out = -q rho_in; wall of tension sigma at areal radius R(tau), U = dR/dtau, Delta = rho_in - rho_out = (1 + q) rho_in.
Israel junction condition in generalized-Lorentz-factor form (Gamma_in - Gamma_out = sigma R / 2):
    Gamma_in  = sigma R / 4 + m / (sigma R^2) - Delta R / (3 sigma)                ,   Gamma_out = Gamma_in - sigma R / 2 ,
    Gamma_in^2 = 1 + U^2 - R^2 / R_dS^2 ,      Gamma_out^2 = 1 + U^2 - m / R - rho_out R^2 / 3 .
Allowed motion needs U^2 >= 0, i.e. Gamma_in^2 >= 1 - R^2 / R_dS^2, which for R < R_dS reads  m >= m_+(R)  or  m <= m_-(R),
    m_pm(R) = (Delta / 3 - sigma^2 / 4) R^3 pm sigma R^2 sqrt(1 - R^2 / R_dS^2).
The two branches meet at R = R_dS, where m = (Delta / 3 - sigma^2 / 4) R_dS^3. For m < m_c a window of forbidden radii exists around R_M ~ (3 m / Delta)^(1/3);
for m > m_c there is none.  (Exact; the check below scans R numerically.)

    python thinwall.py          (writes thinwall.json)
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- physical constants (SI) and Planck 2018 cosmology (matches the benchmark of the dilution-triggered model)
C = 2.99792458e8
G = 6.67430e-11
HBAR = 1.054571817e-34
MPC = 3.0856775814913673e22
MSUN = 1.98847e30
H0 = 67.66e3 / MPC                       # 1/s
OM_L, OM_M = 0.6889, 0.3110
LP = np.sqrt(HBAR * G / C ** 3)
M_PL_GEV = 2.435e18
HBAR_GEV_S = 6.582119569e-25


def m_pm(R, sigma, q, rho_in=3.0, sign=+1):
    Delta = (1.0 + q) * rho_in
    RdS = np.sqrt(3.0 / rho_in)
    return (Delta / 3.0 - sigma ** 2 / 4.0) * R ** 3 + sign * sigma * R ** 2 * np.sqrt(np.clip(1.0 - (R / RdS) ** 2, 0.0, None))


def m_crit(sigma, q, rho_in=3.0):
    """Largest mass at which a window of forbidden radii exists: max over R <= R_dS of m_+(R); equals (Delta/3 - sigma^2/4) R_dS^3 + O(sigma^2) -> (1 + q) R_dS for small tension."""
    RdS = np.sqrt(3.0 / rho_in)
    R = np.linspace(0.5 * RdS, RdS, 2000001)
    return float(m_pm(R, sigma, q, rho_in, +1).max())


def forbidden_window(m, sigma, q, rho_in=3.0, n=400001):
    """Numerical check: radii R < R_dS at which no real motion exists (m_-(R) < m < m_+(R))."""
    RdS = np.sqrt(3.0 / rho_in)
    R = np.linspace(1e-6, RdS, n)
    bad = (m > m_pm(R, sigma, q, rho_in, -1)) & (m < m_pm(R, sigma, q, rho_in, +1))
    return (float(R[bad].min()), float(R[bad].max())) if bad.any() else None


def horizon_radius(m, q, rho_in=3.0):
    """Largest root of f(r) = 1 - m / r + r^2 / L^2 (L^2 = 3 / (q rho_in)); for q = 0 the Schwarzschild value m."""
    if q == 0:
        return float(m)
    L2 = 3.0 / (q * rho_in)
    return float(np.roots([1.0 / L2, 0.0, 1.0, -m]).real.max())


def exterior_roots(m, q):
    """Positive real roots (R_dS = 1) of the horizon condition q r^3 + r - m = 0 of f(r) = 1 - m/r + q r^2, sorted (black-hole horizon first)."""
    if q == 0:
        return [float(m)]
    return sorted(float(x.real) for x in np.roots([q, 0.0, 1.0, -m]) if abs(x.imag) < 1e-12 and x.real > 0)


def m_crit_series(sigma, q):
    """Small-tension expansion of the exact critical mass: (1 + q) - sigma^2 [1/4 - 1/(6 (1 + q))]  (R_dS = 1)."""
    return (1.0 + q) - sigma ** 2 * (0.25 - 1.0 / (6.0 * (1.0 + q)))


def main():
    out = {}
    # ------------------------------------------------------------ critical mass check (moderate tension, then the model's)
    chk = []
    for q in (0.5, 1.0, 2.0):
        for sigma in (0.01, 0.3, 0.9):
            mc = m_crit(sigma, q)
            below = forbidden_window(mc * (1 - 1e-3), sigma, q)
            above = forbidden_window(mc * (1 + 1e-3), sigma, q)
            chk.append(dict(q=q, sigma=sigma, m_c=float(mc), analytic_leading=float((1 + q)), window_below=below, window_above=above,
                            r_h_at_mc=horizon_radius(mc, q)))
    out["critical_mass_check"] = chk
    # ------------------------------------------------------------ finite-tension critical mass and the horizon-area coincidence
    fin = []
    for q in (1.0,):
        for sigma in (0.0, 0.01, 0.1, 0.3, 0.9):
            mc = m_crit(sigma, q) if sigma > 0 else (1.0 + q)
            rh = horizon_radius(mc, q)
            fin.append(dict(q=q, sigma=sigma, m_c_exact=float(mc), m_c_series=float(m_crit_series(sigma, q)), m_c_limit=float(1 + q), rel_dev_from_limit=float(mc / (1 + q) - 1),
                            r_h=float(rh), entropy_ratio_dS_over_BH=float(1.0 / rh ** 2), series_entropy_ratio=float(1.0 + 2.0 * (1 + q - m_crit_series(sigma, q)) / (1 + 3 * q))))
    out["finite_tension"] = fin
    ext = []
    for q in (-0.5, 0.0, 0.5, 1.0, 2.0):
        m = 1.0 + q
        ext.append(dict(q=q, m=m, roots=exterior_roots(m, q), nariai_mass=(float(2.0 / (3.0 * np.sqrt(3.0 * -q))) if q < 0 else None)))
    out["exterior_horizons_at_mc"] = ext
    # ------------------------------------------------------------ the model's wall (physical)
    RdS_gev = 1.0                                                        # set below
    H_L = np.sqrt(OM_L) * H0                                              # 1/s (rho_Lambda = V_f = Omega_L rho_crit)
    RdS_m = C / H_L
    out["R_dS_Gpc"] = RdS_m / (MPC * 1e3)
    out["R_dS_m"] = float(RdS_m)
    H_L_gev = H_L * HBAR_GEV_S
    sigma_model_gev3 = 19.242e-36                                         # tension of the dilution-triggered model at its trigger (meV^3 -> GeV^3)
    out["sigma_model_hat"] = float(sigma_model_gev3 / (M_PL_GEV ** 2 * H_L_gev))       # sigma R_dS / M_Pl^2
    out["wall_thickness_mm"] = 0.239
    out["thickness_over_RdS"] = float(0.239e-3 / RdS_m)
    out["gravity_parameter_model"] = float(1.5 * out["sigma_model_hat"] ** 2 / (2.0))   # 3 sigma^2 / (2 Delta) with Delta = 2 (R_dS = 1), q = 1: sigma^2 * 3/4
    out["gravity_parameter_model"] = float(0.75 * out["sigma_model_hat"] ** 2)
    # ------------------------------------------------------------ cosmological mass and entropy
    rows = []
    rho_crit = 3 * H0 ** 2 / (8 * np.pi * G)                              # kg/m^3 (mass density)
    M_H = 4.0 * np.pi / 3.0 * OM_M * rho_crit * (C / H0) ** 3             # matter inside one Hubble radius today
    M_dS_sphere = 4.0 * np.pi / 3.0 * OM_M * rho_crit * RdS_m ** 3
    for q in (0.25, 0.5, 1.0, 2.0, 4.0):
        mc = m_crit(0.0, q)                                               # in units of R_dS (R_dS = 1)
        Mc = mc / 2.0 * C ** 2 * RdS_m / G                                # m = 2 G M / c^2  in length: M = m c^2 / (2 G) with m in meters (m_hat * R_dS)
        dens_ratio = (1.0 + q) * OM_L / OM_M                              # required mean density over the lump, in units of today's mean matter density
        rows.append(dict(q=q, m_c_over_RdS=float(mc), M_c_kg=float(Mc), M_c_Msun=float(Mc / MSUN), M_c_over_matter_in_Hubble_sphere=float(Mc / M_H),
                         deficit_over_matter_in_RdS_sphere=float(dens_ratio), matter_fraction=float(1.0 / dens_ratio), r_h_over_RdS_at_mc=horizon_radius(mc, q),
                         entropy_ratio_dS_over_BH_at_mc=float((1.0 / horizon_radius(mc, q)) ** 2)))
    out["mass_table"] = rows
    out["M_Hubble_matter_Msun"] = float(M_H / MSUN)
    out["M_matter_in_RdS_sphere_Msun"] = float(M_dS_sphere / MSUN)
    out["S_dS"] = float(np.pi * RdS_m ** 2 / LP ** 2)
    out["Planck_length_m"] = float(LP)
    out["H_Lambda_km_s_Mpc"] = float(np.sqrt(OM_L) * 67.66)
    # matter density of the parent at scale factor a, in units of the vacuum-energy deficit Delta = (1 + q) V_f of the critical region (q = 1)
    a = np.linspace(0.4, 4.0, 19)
    out["matter_over_deficit_vs_a"] = [dict(a=float(x), ratio=float(OM_M * x ** -3 / (2.0 * OM_L))) for x in a]
    out["a_matter_equals_deficit"] = float((OM_M / (2.0 * OM_L)) ** (1.0 / 3.0))
    out["matter_fraction_of_Mc_today"] = float(OM_M / (2.0 * OM_L))
    json.dump(out, open(os.path.join(HERE, "thinwall.json"), "w"), indent=1)
    print("R_dS = %.3f Gpc; sigma_hat(model) = %.2e; gravity parameter %.1e; thickness/R_dS %.1e" % (out["R_dS_Gpc"], out["sigma_model_hat"], out["gravity_parameter_model"], out["thickness_over_RdS"]))
    for r in rows:
        print("q=%.2f: m_c = %.3f R_dS, M_c = %.2e Msun = %.2f x matter in a Hubble sphere, r_h(m_c)/R_dS = %.4f, density ratio %.2f" %
              (r["q"], r["m_c_over_RdS"], r["M_c_Msun"], r["M_c_over_matter_in_Hubble_sphere"], r["r_h_over_RdS_at_mc"], r["deficit_over_matter_in_RdS_sphere"]))
    for f in fin:
        print("finite tension q=%.1f sigma=%.2f: m_c=%.5f (series %.5f), r_h/R_dS=%.5f, S_dS/S_BH=%.5f" % (f["q"], f["sigma"], f["m_c_exact"], f["m_c_series"], f["r_h"], f["entropy_ratio_dS_over_BH"]))
    for e in ext:
        print("exterior roots at m=(1+q)R_dS:", e)
    for c in chk:
        print("check q=%.1f sigma=%.2f m_c=%.5f window below %s above %s r_h/R_dS %.5f" % (c["q"], c["sigma"], c["m_c"], c["window_below"], c["window_above"], c["r_h_at_mc"]))


if __name__ == "__main__":
    main()
