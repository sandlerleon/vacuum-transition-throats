# -*- coding: utf-8 -*-
"""Phase 3D.4: proper-time and instability-window diagnostic on the prescribed Simpson-Visser benchmark ds^2 = -f dt^2 + dr^2/f + R^2 dOmega^2, R^2 = r^2 + a^2, f = 1 - 2m/R.  G = c = 1, m = 1 unless stated.
This is a diagnostic on a background; it does not show formation and does not touch the source problem (3D.3).

Part A  invariant curvature profile (orthonormal Riemann components; validated against Schwarzschild and 3D.1) and the invariantly defined high-curvature region K >= eps K(0).
Part B  proper-time transit through that region: Delta tau(E, L, a/m, eps); horizon-regular (proper time), plus Killing-time behaviour.
Part C  redshift for freely falling observers from photon four-momenta and observer four-velocities (exact, regular across the horizon); static-observer ratio only where static observers exist.
Part D  critical growth rate gamma_c = A_max / Delta tau (no growth rate invented); one COMPUTED instability model: a non-minimally coupled test scalar (xi R coupling) on the wormhole branch, tachyonic bound state => gamma(xi); amplification A = gamma * (residence time).
Part E  verdict.
"""
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}

# ------------------------------------------------------------------------------------------------ Part A: curvature
def geom(r, a, m):
    R = np.sqrt(r ** 2 + a ** 2); f = 1 - 2 * m / R
    Rp = r / R; Rpp = a ** 2 / R ** 3
    fp = 2 * m * r / R ** 3
    fpp = 2 * m * (R ** 2 - 3 * r ** 2) / R ** 5 * 1.0                        # d/dr (2 m r R^-3) = 2m (R^-3 - 3 r^2 R^-5)
    fpp = 2 * m * (1.0 / R ** 3 - 3 * r ** 2 / R ** 5)
    return R, f, Rp, Rpp, fp, fpp
def kretschmann(r, a, m):
    R, f, Rp, Rpp, fp, fpp = geom(r, a, m)
    A = fpp / 2; B = fp * Rp / (2 * R); C = -(fp * Rp / (2 * R) + f * Rpp / R); D = (1 - f * Rp ** 2) / R ** 2
    return 4 * A ** 2 + 8 * B ** 2 + 8 * C ** 2 + 4 * D ** 2
chk = {"3D.1 K(0; a=3,m=1)": (kretschmann(0.0, 3.0, 1.0), 0.06584362139917696), "3D.1 K(0; a=1,m=1)": (kretschmann(0.0, 1.0, 1.0), 16.0), "3D.1 K(0; a=0.5,m=1)": (kretschmann(0.0, 0.5, 1.0), 1472.0),
       "Schwarzschild limit 48 m^2/r^6 at r=2": (kretschmann(2.0, 1e-9, 1.0), 48.0 / 2.0 ** 6)}
for k_, (x, y) in chk.items():
    print("A: %-40s %.10g (reference %.10g)" % (k_, x, y)); assert abs(x / y - 1) < 1e-9
out["A_validation"] = {k_: [float(x), float(y)] for k_, (x, y) in chk.items()}
# Kmax at the throat? and the half-width r_eps of the region K >= eps K(0)
prof = {}
for aa in (0.5, 1.0, 1.9, 3.0, 6.0):
    rr = np.linspace(0, 20 * aa + 20, 200001)
    Kp = kretschmann(rr, aa, 1.0)
    peak_at_throat = bool(np.argmax(Kp) == 0 or Kp[0] >= 0.999999 * Kp.max())
    prof[str(aa)] = dict(K0=float(Kp[0]), Kmax=float(Kp.max()), peak_at_throat=peak_at_throat, curvature_length_K0=float(Kp[0] ** -0.25))
    print("A: a/m = %.1f: K(0) = %.4g, max over r>=0 = %.4g, curvature length K^-1/4 at the throat = %.3f m" % (aa, Kp[0], Kp.max(), Kp[0] ** -0.25))
out["A_profile"] = prof
def r_eps(aa, eps, mm=1.0):
    K0 = kretschmann(0.0, aa, mm); Kmax = max(K0, kretschmann(np.linspace(0, 50, 50001), aa, mm).max())
    g = lambda x: kretschmann(x, aa, mm) - eps * Kmax
    xs = np.linspace(0, 1e3, 200001); vals = g(xs)
    idx = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    return float(brentq(g, xs[idx[-1]], xs[idx[-1] + 1])) if len(idx) else None

# ------------------------------------------------------------------------------------------------ Part B: transit proper time
def U(r, aa, mm, L):
    R = np.sqrt(r ** 2 + aa ** 2); return (1 - 2 * mm / R) * (1 + L ** 2 / R ** 2)
def Lcrit(aa, E, mm=1.0, R0=100.0):
    Rg = np.concatenate([np.linspace(aa, 40 + 4 * aa, 40001), np.geomspace(40 + 4 * aa, R0, 2000)])
    Um = lambda LL: float(np.max((1 - 2 * mm / Rg) * (1 + LL ** 2 / Rg ** 2)))
    return brentq(lambda LL: Um(LL) - E ** 2, 0.0, 80.0, xtol=1e-10)
def transit(aa, E, L, eps, mm=1.0):
    re = r_eps(aa, eps, mm)
    if re is None: return None
    d = lambda x: E ** 2 - U(x, aa, mm, L)
    if min(d(np.linspace(0, re, 4001))) <= 0: return None
    dtau = 2 * quad(lambda x: 1 / np.sqrt(d(x)), 0, re, limit=400)[0]
    # Killing-time elapsed (finite only without a horizon inside the region): dt = E dtau / f
    f_ = lambda x: 1 - 2 * mm / np.sqrt(x ** 2 + aa ** 2)
    dt = None
    if aa > 2 * mm:
        dt = 2 * quad(lambda x: E / (f_(x) * np.sqrt(d(x))), 0, re, limit=400)[0]
    return dict(r_eps=re, dtau=float(dtau), dt_killing=None if dt is None else float(dt))
maps = {}
print("B: transit proper time through the invariant region K >= eps K(0), in units of m (c = G = 1)")
for aa in (0.5, 1.0, 1.9, 3.0, 6.0):
    for eps in (0.5, 0.1, 0.01):
        for E in (1.0, 2.0, 5.0):
            Lc = Lcrit(aa, E)
            for Lfrac in (0.0, 0.5, 0.9):
                res = transit(aa, E, Lfrac * Lc, eps)
                maps["a/m=%.1f eps=%g E=%g L/Lc=%g" % (aa, eps, E, Lfrac)] = None if res is None else dict(res, L=float(Lfrac * Lc), Lc=float(Lc))
    sel = [(eps, maps["a/m=%.1f eps=%g E=1 L/Lc=0" % (aa, eps)]["dtau"], maps["a/m=%.1f eps=%g E=5 L/Lc=0" % (aa, eps)]["dtau"]) for eps in (0.5, 0.1, 0.01)]
    print("B: a/m = %.1f radial: (eps, dtau[E=1], dtau[E=5]) = %s" % (aa, [(e_, round(x, 3), round(y, 3)) for e_, x, y in sel]))
out["B_transit_maps"] = maps
# angular momentum lengthens the transit
for aa in (1.0, 3.0):
    d0 = maps["a/m=%.1f eps=0.1 E=1 L/Lc=0" % aa]["dtau"]; d9 = maps["a/m=%.1f eps=0.1 E=1 L/Lc=0.9" % aa]["dtau"]
    print("B: a/m = %.1f, eps = 0.1, E = 1: dtau(L=0) = %.3f, dtau(L=0.9 Lc) = %.3f (factor %.2f)" % (aa, d0, d9, d9 / d0))
    out.setdefault("B_L_effect", {})[str(aa)] = [d0, d9]
# scaling of the proper transit time with a (E = 1, L = 0, eps = 0.1) and time-to-horizon from R0 = 20
sc = {str(aa): maps["a/m=%.1f eps=0.1 E=1 L/Lc=0" % aa]["dtau"] for aa in (0.5, 1.0, 1.9, 3.0, 6.0)}
out["B_scaling_dtau_vs_a"] = sc
print("B: dtau(E=1, eps=0.1) vs a/m:", {k_: round(v, 3) for k_, v in sc.items()})
# proper time from release at R0 = 20m to the throat and to the horizon; Killing time to approach the horizon
def tau_to(rtarget, aa=1.0, mm=1.0, R0=20.0, E=None):
    r0 = np.sqrt(R0 ** 2 - aa ** 2); E = np.sqrt(1 - 2 * mm / R0) if E is None else E
    return quad(lambda x: 1 / np.sqrt(E ** 2 - (1 - 2 * mm / np.sqrt(x ** 2 + aa ** 2))), rtarget, r0, limit=400)[0]
rh = np.sqrt(3.0)
horizon = dict(tau_R0_to_horizon=tau_to(rh), tau_R0_to_throat=tau_to(0.0))
kt = {}
for delta in (1e-2, 1e-4, 1e-6):
    E = np.sqrt(0.9)
    t_ = quad(lambda x: E / ((1 - 2 / np.sqrt(x ** 2 + 1)) * np.sqrt(E ** 2 - (1 - 2 / np.sqrt(x ** 2 + 1)))), rh + delta, np.sqrt(400 - 1), limit=800, points=[rh + 10 * delta])[0]
    kt[str(delta)] = float(t_)
horizon["killing_time_to_r_h_plus_delta"] = kt
print("B: a/m = 1, rest at 20m: proper time to horizon %.2f m, to the throat %.2f m; Killing time to reach r_h + delta: %s (log divergence, coordinate artifact)" % (horizon["tau_R0_to_horizon"], horizon["tau_R0_to_throat"], {k_: round(v, 1) for k_, v in kt.items()}))
out["B_horizon_clock"] = horizon

# ------------------------------------------------------------------------------------------------ Part C: redshifts from four-momenta
def z_infall(f, E):
    """frequency of a radially ingoing photon (Killing energy w) measured by a radial infaller with energy E: w_obs/w = 1/(E + sqrt(E^2 - f)) ; regular across the horizon."""
    return 1.0 / (E + np.sqrt(E ** 2 - f))
def z_signal_out(f, E):
    """ratio w_inf / w_emitted of an outgoing radial photon emitted by the infaller (only for f > 0): f / (E + sqrt(E^2 - f))."""
    return f / (E + np.sqrt(E ** 2 - f))
# consistency check of the formulas against explicit four-vectors (null condition and normalisation)
fv = 0.37; Ev = 0.9
u_up = np.array([Ev / fv, -np.sqrt(Ev ** 2 - fv)]); g2 = np.diag([-fv, 1 / fv])
k_dn = np.array([-1.0, -1.0 / fv])                                           # ingoing: k_t = -w, k_r = -w/f (w = 1)
nullc = 1 / g2[0, 0] * k_dn[0] ** 2 + 1 / g2[1, 1] * k_dn[1] ** 2
w_obs = -(k_dn @ u_up)
print("C: null condition residual %.1e ; observed/Killing frequency %.6f vs formula %.6f" % (abs(nullc), w_obs, z_infall(fv, Ev)))
assert abs(w_obs - z_infall(fv, Ev)) < 1e-12 and abs(nullc) < 1e-12
red = {}
for aa in (1.0, 3.0):
    rows = []
    for rr in (30.0, 5.0, 2.0, 1.0, 0.0):
        f_ = 1 - 2 / np.sqrt(rr ** 2 + aa ** 2)
        row = dict(r=rr, f=float(f_), infaller_sees_ingoing_photon=float(z_infall(f_, 1.0)))
        if f_ > 0:
            row["static_observer_clock_ratio"] = float(np.sqrt(f_))
            row["signal_from_infaller_received_at_infinity"] = float(z_signal_out(f_, 1.0))
        rows.append(row)
    red["a/m=%.1f" % aa] = rows
    print("C: a/m = %.1f:" % aa, [(rw["r"], round(rw["f"], 3), round(rw["infaller_sees_ingoing_photon"], 3), round(rw.get("signal_from_infaller_received_at_infinity", float("nan")), 3)) for rw in rows])
out["C_redshift"] = red
out["C_statement"] = ("The infalling clock is regular across the horizon (proper time 98 m to the horizon, a few m inside); the 'freezing' is only in signals received outside: outgoing signals from near r_h are redshifted by f/(E+sqrt(E^2-f)) -> 0. "
                      "Static observers do not exist inside |r| < r_h; nothing relates a clock in the second exterior to a clock in the first.")

# ------------------------------------------------------------------------------------------------ Part D: critical growth rates and a computed instability model
dm = 4.925490947e-6                                                          # G M_sun / c^3 in seconds
Amax = {"delta0=1e-3 -> O(1)": np.log(1e3), "delta0=1e-6 -> O(1)": np.log(1e6)}
crit = {}
for aa in (0.5, 1.0, 1.9, 3.0, 6.0):
    for E in (1.0, 5.0):
        d = maps["a/m=%.1f eps=0.1 E=%g L/Lc=0" % (aa, E)]["dtau"]
        crit["a/m=%.1f E=%g" % (aa, E)] = dict(dtau_m=d, gamma_c_per_m=[A / d for A in Amax.values()], gamma_c_per_s_Msun=[A / d / dm for A in Amax.values()], gamma_c_per_s_4e6Msun=[A / d / (dm * 4e6) for A in Amax.values()])
out["D_critical_growth_rates"] = crit
print("D: critical growth rate gamma_c = A_max/dtau (eps = 0.1, radial): a/m=1, E=1: %.3f /m (A=6.9), %.3g /s for 1 M_sun; a/m=1, E=5: %.3f /m" %
      (crit["a/m=1.0 E=1"]["gamma_c_per_m"][0], crit["a/m=1.0 E=1"]["gamma_c_per_s_Msun"][0], crit["a/m=1.0 E=5"]["gamma_c_per_m"][0]))
# computed model: non-minimally coupled test scalar, ell = 0, wormhole branch; V = f[ell(ell+1)/R^2 + xi Ricci] + (f/R)(f a^2/R^3 + 2 m r^2/R^4),  Ricci = -2 a^2 (R - 3m)/R^5
def Vpot(r, aa, mm, ell, xi):
    R = np.sqrt(r ** 2 + aa ** 2); f = 1 - 2 * mm / R
    Ric = -2 * aa ** 2 * (R - 3 * mm) / R ** 5
    return f * (ell * (ell + 1) / R ** 2 + xi * Ric) + (f / R) * (f * aa ** 2 / R ** 3 + 2 * mm * r ** 2 / R ** 4)
# validate against the symbolic 3D.1 potential at a point
rs_, as_, ms_, xs_ = sp.symbols("rs as ms xs", positive=True)
R_ = sp.sqrt(rs_ ** 2 + as_ ** 2); f_ = 1 - 2 * ms_ / R_
Ric_s = -2 * as_ ** 2 * (R_ - 3 * ms_) / R_ ** 5
V_s = f_ * (xs_ * Ric_s) + f_ * sp.diff(f_ * sp.diff(R_, rs_), rs_) / R_
v1 = float(V_s.subs({rs_: 1.3, as_: 3.0, ms_: 1.0, xs_: 0.7})); v2 = float(Vpot(1.3, 3.0, 1.0, 0, 0.7))
print("D: potential check (symbolic vs closed form): %.12g vs %.12g" % (v1, v2)); assert abs(v1 - v2) < 1e-12
def ground(aa, mm, ell, xi, Lb=150.0, h=0.02):
    ss = np.arange(-Lb, Lb + h / 2, h)
    sol = solve_ivp(lambda s, y: [1 - 2 * mm / np.sqrt(y[0] ** 2 + aa ** 2)], (0, Lb), [0.0], dense_output=True, rtol=1e-12, atol=1e-12)
    rpos = sol.sol(np.abs(ss))[0] * np.sign(ss)
    V = Vpot(rpos, aa, mm, ell, xi)
    w = eigh_tridiagonal(2 / h ** 2 + V, -np.ones(len(ss) - 1) / h ** 2, select="i", select_range=(0, 0), eigvals_only=True)
    return float(w[0])
gm = {}
for name, aa, mm in (("wormhole a/m=3", 3.0, 1.0), ("wormhole a/m=6", 6.0, 1.0), ("Ellis m=0, a=1", 1.0, 0.0)):
    rows = []
    for xi in (0.0, 1 / 6, 1.0, 3.0, 6.0, 12.0):
        lam = ground(aa, mm, 0, xi)
        rows.append(dict(xi=xi, lowest_eigenvalue=lam, gamma_killing=float(np.sqrt(-lam)) if lam < -1e-6 else 0.0))
    # critical xi by bisection on the existence of a bound state (lowest eigenvalue < -1e-5)
    lo, hi = 0.0, 12.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if ground(aa, mm, 0, mid) < -1e-5: hi = mid
        else: lo = mid
    gm[name] = dict(rows=rows, xi_critical=hi)
    print("D: %-18s xi_crit(l=0) = %.3f ; gamma_t(xi) = %s" % (name, hi, [(round(r_["xi"], 3), round(r_["gamma_killing"], 4)) for r_ in rows]))
out["D_tachyonic_test_field"] = gm
# amplification A = gamma * residence time (Killing time of the packet in the region K >= eps K(0)); delta0 = 1e-3 -> A_max = 6.9
Atab = []
for name, aa, mm in (("wormhole a/m=3", 3.0, 1.0), ("wormhole a/m=6", 6.0, 1.0)):
    xi = 12.0
    gam = [r_["gamma_killing"] for r_ in gm[name]["rows"] if r_["xi"] == xi][0]
    for E in (1.0, 1.5, 2.0, 5.0, 20.0):
        res = transit(aa, E, 0.0, 0.1, mm)
        if res is None or res["dt_killing"] is None: continue
        Atab.append(dict(case=name, xi=xi, E=E, gamma_killing=gam, dt_killing=res["dt_killing"], dtau=res["dtau"], A=gam * res["dt_killing"]))
        print("D: %s, xi = 12: E = %5.1f: Killing residence time %.3f m, gamma = %.4f /m -> A = %.3f (A_max = 6.9, gamma_c = %.3f)" % (name, E, res["dt_killing"], gam, gam * res["dt_killing"], 6.9 / res["dt_killing"]))
floor = {}
for name, aa, mm in (("wormhole a/m=3", 3.0, 1.0), ("wormhole a/m=6", 6.0, 1.0)):
    re = r_eps(aa, 0.1, mm)
    tl = 2 * quad(lambda x: 1.0 / (1 - 2 * mm / np.sqrt(x ** 2 + aa ** 2)), 0, re, limit=400)[0]
    gam = [r_["gamma_killing"] for r_ in gm[name]["rows"] if r_["xi"] == 12.0][0]
    floor[name] = dict(light_crossing_killing_time=float(tl), gamma=gam, A_floor=float(gam * tl), A_max=float(np.log(1e3)))
    print("D: %s: light-crossing Killing time of the region (floor for ANY E) = %.3f m ; A_floor = gamma*t_light = %.3f (A_max = 6.9)" % (name, tl, gam * tl))
out["D_killing_time_floor"] = floor
out["D_amplification"] = Atab
out["D_interpretation"] = ("A = gamma * residence time measures the growth, during the packet's stay in the high-curvature region, of a perturbation seeded at entry in an unstable Killing-time mode. The mode is a global bound state of the TEST FIELD: a toy "
                           "instability (tachyonic xi R coupling), NOT the instability of the source of the benchmark. For the forced phantom source no growth rate is derived here (literature: unstable radial mode, not recomputed).")
json.dump(out, open(os.path.join(HERE, "phase3d4_propertime.json"), "w"), indent=1, default=str)
