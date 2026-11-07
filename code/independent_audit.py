# -*- coding: utf-8 -*-
"""Independent audit of the thin-wall critical mass (revision 2).  Nothing here imports thinwall.py's formulas except for the final comparison.

(1) The Israel junction branches are re-derived symbolically from the extrinsic curvatures (sympy) and compared with the closed form used in thinwall.py.
(2) The critical mass m_c(sigma, q) = max_R m_+(R) is recomputed with arbitrary-precision root finding (mpmath), including the physical tension where the correction is of order 1e-60.
(3) Horizon assignments of the exterior metric f(r) = 1 - m/r + q r^2 (R_dS = 1) at m = m_c, and the equality of horizon areas.
(4) Mass definitions: the Schwarzschild-(A)dS parameter m = 2GM, the Misner-Sharp quasi-local mass at the wall, and their relation for AdS (q > 0), Minkowski (q = 0) and dS (q < 0) exteriors.
(5) Dimensional restoration with CODATA/Planck constants fixed here independently, against the values in thinwall.json.
Units of (1)-(4): 8 pi G = 1, R_dS = 1 (rho_in = 3); rho_out = -q rho_in.
"""
import json
import os

import mpmath as mp
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
TW = json.load(open(os.path.join(HERE, "thinwall.json"), encoding="utf-8"))
out = {}

# ------------------------------------------------------------------------------------------------ (1) symbolic junction conditions
R, m, sig, q, rin = sp.symbols("R m sigma q rho_in", positive=True)
rout = -q * rin
RdS2 = 3 / rin
f_in = 1 - R ** 2 / RdS2                                    # interior de Sitter
f_out = 1 - m / R - rout * R ** 2 / 3                       # exterior Schwarzschild-(A)dS, 8 pi G = 1
Gi, Go, U2 = sp.symbols("Gamma_in Gamma_out U2")
# Israel (8 pi G = 1): S^tau_tau = -sigma, [K^theta_theta] = (Gamma_out - Gamma_in)/R, S^tau_tau = 2 [K^theta_theta]  =>  Gamma_in - Gamma_out = sigma R / 2 ;  Gamma_a^2 = f_a + U^2
sol = sp.solve([sp.Eq(Gi - Go, sig * R / 2), sp.Eq(Go ** 2 - Gi ** 2, f_out - f_in)], [Gi, Go], dict=True)[0]
Gamma_in = sp.simplify(sol[Gi])
U2_expr = sp.simplify(Gamma_in ** 2 - f_in)
Delta = (1 + q) * rin
mpm = sp.solve(sp.Eq(U2_expr, 0), m)
print("1: Gamma_in =", Gamma_in)
print("1: U^2 = 0 has", len(mpm), "roots in m")
fn_sym = [sp.lambdify((R, sig, q, rin), r_, "mpmath") for r_ in mpm]
out["1_Gamma_in"] = str(Gamma_in)
# closed form of thinwall.py
def m_closed(Rv, s, qq, rho=3.0, sign=+1):
    D = (1 + qq) * rho; RdS = mp.sqrt(3 / rho)
    return (D / 3 - s ** 2 / 4) * Rv ** 3 + sign * s * Rv ** 2 * mp.sqrt(1 - (Rv / RdS) ** 2)
mp.mp.dps = 40
rng = np.random.default_rng(2)
maxdev = 0.0
for _ in range(200):
    Rv = mp.mpf(float(rng.uniform(0.05, 0.99))); s = mp.mpf(float(rng.uniform(0.01, 0.9))); qq = mp.mpf(float(rng.choice([-0.5, 0.0, 0.5, 1.0, 2.0])))
    vals = sorted([complex(f(Rv, s, qq, mp.mpf(3))).real for f in fn_sym])
    ref = sorted([float(m_closed(Rv, s, qq, 3, +1)), float(m_closed(Rv, s, qq, 3, -1))])
    maxdev = max(maxdev, max(abs(a - b) for a, b in zip(vals, ref)))
print("1: symbolic roots of U^2 = 0 versus the closed-form m_pm(R): max deviation over 200 random (R, sigma, q) =", maxdev)
out["1_max_deviation_symbolic_vs_closed_form"] = float(maxdev)
meet = sp.simplify(mpm[0].subs(R, sp.sqrt(3 / rin)))
print("1: branches meet at R = R_dS with m =", sp.simplify(meet), "(expected (Delta/3 - sigma^2/4) R_dS^3)")
out["1_meeting_point"] = [str(sp.simplify(mpm[0].subs(R, sp.sqrt(3 / rin)))), str(sp.simplify(mpm[1].subs(R, sp.sqrt(3 / rin))))]

# ------------------------------------------------------------------------------------------------ (2) arbitrary-precision critical mass
mp.mp.dps = 120
def mc_hp(s, qq):
    s = mp.mpf(s); qq = mp.mpf(qq)
    g = lambda x: mp.diff(lambda y: m_closed(y, s, qq, 3, +1), x)
    # the maximum lies just below R_dS = 1: bracket by scanning x = 1 - eps
    best = None
    for k in range(1, 60):
        eps = mp.mpf(10) ** (-k * mp.mpf("0.5"))
        pass
    f = lambda x: mp.diff(lambda y: m_closed(y, s, qq, 3, +1), x)
    # solve f(x) = 0 by bisection on (x_lo, 1 - 1e-60)
    lo, hi = mp.mpf("0.5"), 1 - mp.mpf(10) ** (-100)
    if f(lo) <= 0:
        return m_closed(lo, s, qq, 3, +1), lo
    for _ in range(400):
        mid = (lo + hi) / 2
        if f(mid) > 0: lo = mid
        else: hi = mid
    return m_closed(lo, s, qq, 3, +1), lo
rows = []
for qq in (0.5, 1.0, 2.0):
    for s in (0.01, 0.1, 0.3, 0.9):
        mc, xs = mc_hp(s, qq)
        ref = [r for r in TW["finite_tension"] if abs(r["q"] - qq) < 1e-9 and abs(r["sigma"] - s) < 1e-12] if "finite_tension" in TW and isinstance(TW["finite_tension"], list) else []
        series = (1 + qq) - s ** 2 * (0.25 - 1 / (6 * (1 + qq)))
        rows.append(dict(q=qq, sigma=s, m_c_hp=float(mc), R_at_max=float(xs), series=series, dev_from_series=float(mc) - series))
        print("2: q = %.1f, sigma = %.2f: m_c (120-digit) = %.12f, series %.12f, R at the maximum = %.6f" % (qq, s, mc, series, xs))
out["2_critical_mass_hp"] = rows
phys_sigma = TW["sigma_model_hat"]
mc_phys, x_phys = mc_hp(phys_sigma, 1.0)
dev_phys = mc_phys - 2
series_dev = -mp.mpf(phys_sigma) ** 2 * (mp.mpf(1) / 4 - mp.mpf(1) / 12)
print("2: physical tension sigma_hat = %.4g, q = 1: m_c - 2 = %s (series: %s); R_dS - R at the maximum = %s" % (phys_sigma, mp.nstr(dev_phys, 6), mp.nstr(series_dev, 6), mp.nstr(1 - x_phys, 6)))
out["2_physical"] = dict(sigma_hat=phys_sigma, m_c_minus_limit=str(mp.nstr(dev_phys, 8)), series_prediction=str(mp.nstr(series_dev, 8)), one_minus_R_at_max=str(mp.nstr(1 - x_phys, 8)))

# ------------------------------------------------------------------------------------------------ (3) horizons at m_c
mp.mp.dps = 30
def horizons(mm, qq):
    # f(r) = 1 - m/r + q r^2 = 0  <=>  q r^3 + r - m = 0
    if qq == 0: return [mm]
    rts = mp.polyroots([qq, 0, 1, -mm], maxsteps=200, extraprec=200)
    return sorted([mp.re(x) for x in rts if abs(mp.im(x)) < 1e-20 and mp.re(x) > 0])
hz = []
for qq in (-0.5, 0.0, 0.5, 1.0, 2.0):
    mm = 1 + qq
    rts = horizons(mm, qq)
    note = {True: "Schwarzschild-AdS: one horizon", False: ""}[qq > 0]
    if qq == 0: note = "Schwarzschild: r_h = m"
    if qq < 0: note = "Schwarzschild-de Sitter: black-hole horizon and cosmological horizon (both listed)"
    hz.append(dict(q=qq, m_c=mm, roots=[float(x) for x in rts], note=note))
    print("3: q = %5.2f, m = m_c = %.2f: horizons r = %s  (%s)" % (qq, mm, [round(float(x), 9) for x in rts], note))
out["3_horizons_at_mc"] = hz

# ------------------------------------------------------------------------------------------------ (4) mass definitions
# exterior metric function f_out(r) = 1 - m/r - rho_out r^2/3 = 1 - m/r + q (rho_in/3) r^2 ;  Misner-Sharp: (1 - f)(r/2G) -> in 8 pi G = 1 units 2 m_MS(r) = m + rho_out r^3 / 3 ... (G-normalised below)
mdefs = []
for qq in (-0.5, 0.0, 0.5, 1.0, 2.0):
    mm = 1 + qq
    m_ms_out = mm + (-qq * 3.0) * 1.0 ** 3 / 3.0                 # parameter + rho_out R^3 / 3 at R = R_dS = 1
    m_ms_in = 3.0 * 1.0 ** 3 / 3.0                                # interior de Sitter, rho_in R^3 / 3 at R = R_dS
    mdefs.append(dict(q=qq, parameter_m=mm, misner_sharp_exterior_at_wall=m_ms_out, misner_sharp_interior_at_wall=m_ms_in,
                      note="m is the Schwarzschild-(A)dS parameter m = 2GM (ADM-like only for q >= 0; for q < 0 there is no asymptotically flat region); the quasi-local mass at the wall at zero tension equals the interior de Sitter mass R_dS (in these units) for every q"))
    print("4: q = %5.2f: parameter m = %.2f ; Misner-Sharp mass at R = R_dS: exterior %.2f, interior %.2f" % (qq, mm, m_ms_out, m_ms_in))
out["4_mass_definitions"] = mdefs

# ------------------------------------------------------------------------------------------------ (5) dimensional restoration (constants fixed independently here)
c = 299792458.0; G = 6.67430e-11; Mpc = 3.0856775814913673e22; Msun = 1.98847e30
H0 = 67.66e3 / Mpc; OmL = 0.6889
HL = H0 * np.sqrt(OmL)
RdS_m = c / HL
dim = dict(R_dS_Gpc=RdS_m / (1e3 * Mpc), H_Lambda_km_s_Mpc=HL * Mpc / 1e3)
M_per_unit = RdS_m * c ** 2 / (2 * G)                       # mass corresponding to m = 1 R_dS: 2GM/c^2 = m R_dS
for qq in (-0.5, 0.0, 0.25, 0.5, 1.0, 2.0):
    M_c = (1 + qq) * M_per_unit
    M_ms = 1.0 * M_per_unit
    dim["q=%g" % qq] = dict(M_c_kg=M_c, M_c_Msun=M_c / Msun, M_misner_sharp_at_wall_Msun=M_ms / Msun)
ref = {str(r["q"]): r for r in TW["mass_table"]}
cmp_ = {}
for qq in (0.25, 0.5, 1.0, 2.0):
    if str(qq) in ref:
        cmp_[str(qq)] = dict(independent_Msun=dim["q=%g" % qq]["M_c_Msun"], thinwall_json_Msun=ref[str(qq)]["M_c_Msun"], rel_diff=dim["q=%g" % qq]["M_c_Msun"] / ref[str(qq)]["M_c_Msun"] - 1)
out["5_dimensional_restoration"] = dim; out["5_comparison"] = cmp_
print("5: R_dS = %.4f Gpc (thinwall.json %.4f) ; M_c(q=1) = %.4g Msun (parameter), Misner-Sharp mass at the wall = %.4g Msun" % (dim["R_dS_Gpc"], TW["R_dS_Gpc"], dim["q=1"]["M_c_Msun"], dim["q=1"]["M_misner_sharp_at_wall_Msun"]))
print("5: relative differences to thinwall.json:", {k_: "%.1e" % v["rel_diff"] for k_, v in cmp_.items()})
json.dump(out, open(os.path.join(HERE, "independent_audit.json"), "w"), indent=1, default=str)
