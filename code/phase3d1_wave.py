# -*- coding: utf-8 -*-
"""Phase 3D.1: quantum-wave (test scalar field) propagation through a hypothetical recycling throat.  Benchmark geometry: the Simpson-Visser "black-bounce" family
    ds^2 = -f dt^2 + dr^2/f + R^2 dOmega^2 ,  R^2 = r^2 + a^2 ,  f = 1 - 2m/R ,  r in (-inf, +inf),
which contains all three classes with one parameter:  a > 2m  traversable two-way wormhole (no horizon);  a = 2m one-way extremal transition;  a < 2m regular black hole whose interior bounces into an expanding (white-hole) region
of a second asymptotic universe (a = 0: Schwarzschild).  G = c = hbar = 1, kappa = 8 pi.  The geometry is a mathematical benchmark: it is NOT shown to solve physically admissible field equations (Part C).

Part A  curvature (symbolic), stress-energy required, energy conditions, horizons, causal passability (radial geodesics).
Part B  Klein-Gordon  (Box - mu^2 - xi Ricci) Phi = 0  ->  u'' + (omega^2 - V_l) u = 0 in the tortoise coordinate r*, dr*/dr = 1/f;  V_l = f [l(l+1)/R^2 + mu^2 + xi Ricci] + (1/R) d^2R/dr*^2.
        Transmission T (wormhole) or greybody Gamma (black bounce), reflection, flux conservation, convergence in the box size.
Part C  backreaction requirements: energy conditions, volume integral, classical wave stress, order-of-magnitude quantum-inequality bound.
"""
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
kappa = 8 * np.pi

# ------------------------------------------------------------------------------------------------ Part A: curvature and stress-energy (sympy)
t, r, th, ph = sp.symbols("t r theta phi", real=True)
a, m = sp.symbols("a m", positive=True)
Rs = sp.sqrt(r ** 2 + a ** 2)
fs = 1 - 2 * m / Rs
X = [t, r, th, ph]
g = sp.diag(-fs, 1 / fs, Rs ** 2, Rs ** 2 * sp.sin(th) ** 2)
gi = g.inv()
Gam = [[[sum(gi[i, l] * (sp.diff(g[l, j], X[k]) + sp.diff(g[l, k], X[j]) - sp.diff(g[j, k], X[l])) for l in range(4)) / 2 for k in range(4)] for j in range(4)] for i in range(4)]
def Riem(i, j, k, l):
    e = sp.diff(Gam[i][j][l], X[k]) - sp.diff(Gam[i][j][k], X[l])
    e += sum(Gam[i][k][s] * Gam[s][j][l] - Gam[i][l][s] * Gam[s][j][k] for s in range(4))
    return e
Ric = sp.Matrix(4, 4, lambda j, l: sp.simplify(sum(Riem(i, j, i, l) for i in range(4))))
Rsc = sp.simplify(sum(gi[i, i] * Ric[i, i] for i in range(4)))
Gmix = sp.Matrix(4, 4, lambda i, j: sp.simplify(sum(gi[i, k] * Ric[k, j] for k in range(4)) - (Rsc / 2 if i == j else 0)))
rho_s = sp.simplify(-Gmix[0, 0] / kappa)
pr_s = sp.simplify(Gmix[1, 1] / kappa)
pt_s = sp.simplify(Gmix[2, 2] / kappa)
print("A: Ricci scalar =", Rsc)
print("A: rho =", rho_s, "\n   p_r =", pr_s, "\n   p_t =", pt_s)
out["A_ricci_scalar"] = str(Rsc); out["A_rho"] = str(rho_s); out["A_pr"] = str(pr_s); out["A_pt"] = str(pt_s)
# Kretschmann via the standard static spherical formula is long; regularity is checked through the curvature components at r = 0 and as r -> infinity numerically
Kret = sp.simplify(sum(Riem(i, j, k, l) * Riem(i, j, k, l) * gi[j, j] * gi[k, k] * gi[l, l] * g[i, i] for i in range(4) for j in range(4) for k in range(4) for l in range(4) if Riem(i, j, k, l) != 0))
Kfun = sp.lambdify((r, a, m), Kret, "numpy")
print("A: Kretschmann K(r=0) for (a,m) = (3,1), (1,1), (0.5,1):", [float(Kfun(0.0, aa, 1.0)) for aa in (3.0, 1.0, 0.5)], " ; max over r in [-50,50] finite:", max(float(Kfun(x, 1.0, 1.0)) for x in np.linspace(-50, 50, 2001)))
out["A_kretschmann_r0"] = {str(aa): float(Kfun(0.0, aa, 1.0)) for aa in (3.0, 1.0, 0.5)}
rho_f = sp.lambdify((r, a, m), rho_s, "numpy"); pr_f = sp.lambdify((r, a, m), pr_s, "numpy"); pt_f = sp.lambdify((r, a, m), pt_s, "numpy")
# energy conditions.  Closed form: rho + p_r = -2 f a^2 /(kappa R^4) = -(2/kappa) f R''/R  ;  the radial null combination T_kk = (T^r_r - T^t_t)/f = -2 a^2/(kappa R^4) < 0 at EVERY r (all m, all a > 0)
chk = sp.simplify(rho_s + pr_s + 2 * fs * a ** 2 / (kappa * Rs ** 4))
print("A: rho + p_r + 2 f a^2/(kappa R^4) =", chk, "(0 => radial NEC T_kk = -2 a^2/(kappa R^4) < 0 everywhere)")
out["A_nec_closed_form_residual"] = str(chk)
ec = {}
for aa in (0.5, 1.0, 2.5, 3.0, 6.0):
    rr = np.linspace(-30, 30, 6000)               # even point count: r = 0 is not a grid point (a = 2m is degenerate there)
    with np.errstate(all="ignore"):
        rho = rho_f(rr, aa, 1.0); pr = pr_f(rr, aa, 1.0); pt = pt_f(rr, aa, 1.0)
    Tkk = -2 * aa ** 2 / (kappa * (rr ** 2 + aa ** 2) ** 2)
    ec[str(aa)] = dict(Tkk_radial_max=float(Tkk.max()), Tkk_radial_min=float(Tkk.min()), Tkk_at_throat=float(-2 / (kappa * aa ** 2)),
                       rho_pt_min=float(np.nanmin(rho + pt)), rho_at_throat=float(rho_f(1e-9, aa, 1.0)), pr_at_throat=float(pr_f(1e-9, aa, 1.0)), pt_at_throat=float(pt_f(1e-9, aa, 1.0)))
    print("A: a/m = %.1f: radial null T_kk in [%.3g, %.3g] (always negative); throat: rho = %.3g, p_r = %.3g, p_t = %.3g" % (aa, Tkk.min(), Tkk.max(), ec[str(aa)]["rho_at_throat"], ec[str(aa)]["pr_at_throat"], ec[str(aa)]["pt_at_throat"]))
out["A_energy_conditions"] = ec
# horizons and causal passability (radial timelike geodesic from rest at large r)
def fnum(rr, aa, mm=1.0): return 1 - 2 * mm / np.sqrt(rr ** 2 + aa ** 2)
caus = {}
for aa in (3.0, 2.0, 1.0):
    if aa > 2.0:
        horizon = None; fmin = float(fnum(0.0, aa))
    else:
        horizon = float(np.sqrt(4.0 - aa ** 2)); fmin = float(1 - 2.0 / aa)
    # radial geodesic dr/dtau = -sqrt(E^2 - f), E = 1 (released from rest at infinity): crosses r = 0 and reaches r = -R_end if E^2 - f > 0 everywhere
    rr = np.linspace(-200, 200, 400001)
    E2mf = 1.0 - fnum(rr, aa)
    ok = bool(np.all(E2mf > 0))
    # proper time from r = +50 to r = -50 (not through coordinate singularities: dtau = dr / sqrt(E^2 - f), regular at horizons)
    tau = quad(lambda x: 1.0 / np.sqrt(1.0 - fnum(x, aa)), -50, 50, limit=400)[0]
    caus[str(aa)] = dict(horizon_radius_r=horizon, f_min=fmin, radial_geodesic_crosses=ok, proper_time_50_to_minus50=float(tau))
    print("A: a/m = %.1f: horizon r_h = %s, f_min = %.3f, E=1 radial geodesic crosses the throat: %s, proper time (+50 -> -50) = %.3f" % (aa, horizon, fmin, ok, tau))
out["A_causal"] = caus
# null expansions: outgoing/ingoing radial null congruences have theta_+- = +-(2/R) (dR/dr) sqrt?  -> theta_+- ~ (2/R) * (+-f) * (r/R) in r-ingoing coordinates; trapped iff f < 0
out["A_trapped_iff"] = "f < 0, i.e. R < 2m, only for a < 2m (|r| < sqrt(4 m^2 - a^2))"

# ------------------------------------------------------------------------------------------------ Part B: wave scattering
Rn = lambda x, aa: np.sqrt(x ** 2 + aa ** 2)
def make_V(aa, mm, ell, mu, xi):
    rs_, as_, ms_ = sp.symbols("rs as ms", positive=True)
    R_ = sp.sqrt(rs_ ** 2 + as_ ** 2); f_ = 1 - 2 * ms_ / R_
    d2R = f_ * sp.diff(f_ * sp.diff(R_, rs_), rs_)                     # d^2 R / dr*^2
    Ricci_ = Rsc.subs({r: rs_, a: as_, m: ms_})
    V = f_ * (ell * (ell + 1) / R_ ** 2 + mu ** 2 + xi * Ricci_) + d2R / R_
    return sp.lambdify(rs_, V.subs({as_: aa, ms_: mm}), "numpy")
def scatter(aa, mm, ell, omega, mu=0.0, xi=0.0, L=None, rtol=1e-11):
    k = np.sqrt(omega ** 2 - mu ** 2)
    V = make_V(aa, mm, ell, mu, xi)
    fN = lambda x: 1 - 2 * mm / np.sqrt(x ** 2 + aa ** 2)
    L = L if L is not None else 2000.0 / k
    two_way = aa > 2 * mm
    # r(r*) from the reference point r*=0 (wormhole: r=0; black bounce: choose r = r_h + 1 as a convenient origin: only the lateral shift of r* changes)
    if two_way:
        r0 = 0.0
        rm = solve_ivp(lambda s, y: [fN(y[0])], (0, -L), [r0], rtol=1e-12, atol=1e-12).y[0, -1]
        rstart = rm
    else:
        rh = np.sqrt(4 * mm ** 2 - aa ** 2); r0 = rh + 1.0
        rm = solve_ivp(lambda s, y: [fN(y[0])], (0, -L), [r0], rtol=1e-12, atol=1e-30).y[0, -1]
        rstart = rm
    kh = omega if not two_way else k                                   # near-horizon wavenumber (V -> 0) or asymptotic wavenumber at the other end
    Vs = float(V(rstart)) if two_way else 0.0
    kloc = np.sqrt(omega ** 2 - Vs) if two_way else omega
    def rhs(s, y):
        x, ur, ui, vr, vi = y
        Vv = float(V(x))
        return [fN(x), vr, vi, -(omega ** 2 - Vv) * ur, -(omega ** 2 - Vv) * ui]
    y0 = [rstart, np.cos(-kloc * (-L)), np.sin(-kloc * (-L)), 0.0, 0.0]          # u = exp(-i kloc r*) at r* = -L
    u0 = np.exp(-1j * kloc * (-L)); v0 = -1j * kloc * u0
    y0 = [rstart, u0.real, u0.imag, v0.real, v0.imag]
    sol = solve_ivp(rhs, (-L, L), y0, method="DOP853", rtol=rtol, atol=1e-13)
    x, ur, ui, vr, vi = sol.y[:, -1]
    u = ur + 1j * ui; v = vr + 1j * vi
    kend = np.sqrt(max(omega ** 2 - float(V(x)), 1e-30))
    al = 0.5 * (u + v / (-1j * kend)) * np.exp(1j * kend * L)                     # coefficient of exp(-i k r*)
    be = 0.5 * (u - v / (-1j * kend)) * np.exp(-1j * kend * L)                    # coefficient of exp(+i k r*)
    # flux: k_inf (|al|^2 - |be|^2) = kloc  (Wronskian)
    Tprob = kloc / (k * abs(al) ** 2)
    Rprob = abs(be) ** 2 / abs(al) ** 2
    return dict(T=float(Tprob), R=float(Rprob), flux_residual=float(Tprob + Rprob - 1.0))
rows = []
cases = [("wormhole a/m=3", 3.0, 1.0), ("wormhole a/m=6", 6.0, 1.0), ("Ellis m=0 a=3", 3.0, 0.0), ("black bounce a/m=1", 1.0, 1.0), ("black bounce a/m=1.9", 1.9, 1.0), ("Schwarzschild a=0", 1e-8, 1.0)]
for name, aa, mm in cases:
    for ell in (0, 1, 2):
        for om in (0.05, 0.1, 0.2, 0.4, 0.8, 1.6):
            if mm == 0.0 and aa > 0:
                pass
            try:
                s_ = scatter(aa, mm, ell, om)
            except Exception as e:
                s_ = dict(T=float("nan"), R=float("nan"), flux_residual=float("nan"))
            rows.append(dict(case=name, a=aa, m=mm, ell=ell, omega=om, **s_))
print("B: scattering table done (%d rows); worst |T + R - 1| = %.2e" % (len(rows), np.nanmax([abs(rw["flux_residual"]) for rw in rows])))
out["B_scattering"] = rows
for name in ("wormhole a/m=3", "wormhole a/m=6", "Ellis m=0 a=3", "black bounce a/m=1", "black bounce a/m=1.9", "Schwarzschild a=0"):
    sel = [rw for rw in rows if rw["case"] == name and rw["ell"] == 0]
    print("B: %-22s l=0: T(omega) = %s" % (name, ", ".join("%.3g@%.2f" % (rw["T"], rw["omega"]) for rw in sel)))
# convergence in the box size and tolerance for one case
conv = []
for L_fac in (1000.0, 2000.0, 4000.0):
    s1 = scatter(3.0, 1.0, 1, 0.4, L=L_fac / 0.4)
    conv.append(dict(L_times_k=L_fac, T=s1["T"], R=s1["R"], flux_residual=s1["flux_residual"]))
s_t = scatter(3.0, 1.0, 1, 0.4, rtol=1e-13)
print("B: convergence (wormhole a/m=3, l=1, omega=0.4):", [(c["L_times_k"], round(c["T"], 7)) for c in conv], " rtol 1e-13:", round(s_t["T"], 7))
out["B_convergence"] = dict(box=conv, tight_tolerance_T=s_t["T"])
# low-frequency scaling T_l ~ omega^p
lowfreq = {}
for name, aa, mm in (("wormhole a/m=3", 3.0, 1.0), ("Ellis m=0 a=3", 3.0, 0.0)):
    for ell in (0, 1):
        om = np.array([0.02, 0.04, 0.08])
        Tt = np.array([scatter(aa, mm, ell, o_)["T"] for o_ in om])
        p_ = np.polyfit(np.log(om), np.log(Tt), 1)[0]
        lowfreq["%s l=%d" % (name, ell)] = dict(T=Tt.tolist(), exponent=float(p_))
        print("B: %s l=%d: T(omega) at omega = 0.02, 0.04, 0.08: %s, fitted exponent %.2f" % (name, ell, np.round(Tt, 6), p_))
out["B_lowfreq"] = lowfreq
# massive field: Compton wavelength vs throat radius (wormhole a/m=3): mu a = 0.3, 1.0
massive = []
for mu in (0.1, 0.3):
    for om_fac in (1.05, 1.5, 3.0):
        om = mu * om_fac
        s_ = scatter(3.0, 1.0, 0, om, mu=mu)
        massive.append(dict(mu=mu, mu_a=mu * 3.0, omega=om, **s_))
        print("B: massive wormhole a/m=3: mu a = %.2f, omega/mu = %.2f: T = %.4g, flux residual %.1e" % (mu * 3.0, om_fac, s_["T"], s_["flux_residual"]))
out["B_massive"] = massive
# bound-state / well check of the potential (wormhole a/m=3, l=0, mu=0.3): does V dip below mu^2?
V_ = make_V(3.0, 1.0, 0, 0.3, 0.0)
xs = np.linspace(-40, 40, 8001)
out["B_potential_minus_mu2_min"] = float((V_(xs) - 0.09).min())

# ------------------------------------------------------------------------------------------------ Part C: backreaction requirements
volume = {}
for aa in (3.0, 6.0):
    f_ = lambda x: 1 - 2.0 / np.sqrt(x ** 2 + aa ** 2)
    integrand = lambda x: (-2 * f_(x) * aa ** 2 / (kappa * (x ** 2 + aa ** 2) ** 2)) * 4 * np.pi * (x ** 2 + aa ** 2) / np.sqrt(f_(x))
    I = quad(integrand, -400, 400, limit=800)[0]
    # scale: max |rho + p_r| a^2 kappa
    xx = np.linspace(-20 * aa, 20 * aa, 40001)
    c_ = float(np.max(np.abs(-2 * f_(xx) * aa ** 2 / (kappa * (xx ** 2 + aa ** 2) ** 2))) * aa ** 2 * kappa)
    volume[str(aa)] = dict(volume_integral_rho_plus_pr=float(I), nec_density_scale_c=c_)
    print("C: a/m = %.0f: integral (rho + p_r) dV over the wormhole = %.4g (Planck units, tends to a negative constant), max |rho + p_r| kappa a^2 = %.3g" % (aa, I, c_))
out["C_volume"] = volume
# Ford-Roman-type quantum inequality (flat-space massless minimally coupled scalar, Lorentzian sampling): |rho| <= 3 N / (32 pi^2 tau0^4) with tau0 <~ a ; requirement rho ~ c/(kappa a^2)
c_req = max(v["nec_density_scale_c"] for v in volume.values())
Nreq = (32 * np.pi ** 2 / 3) * c_req / kappa                                  # N >= (32 pi^2/3) c a^2/kappa  = Nreq * (a/l_P)^2
print("C: QI estimate: required number of massless scalar species N >= %.3g * (a / l_P)^2" % Nreq)
out["C_QI_species_coefficient"] = float(Nreq)
sizes = {"1 m": 1.0 / 1.616e-35, "1 km": 1e3 / 1.616e-35, "R_sun": 6.96e8 / 1.616e-35}
out["C_QI_species_needed"] = {k_: float(Nreq * v ** 2) for k_, v in sizes.items()}
print("C: species needed for throat scale a:", {k_: "%.2e" % (Nreq * v ** 2) for k_, v in sizes.items()})
# classical minimally coupled wave: T_kk = (k.dPhi)^2 >= 0  -> NEC satisfied; non-minimal xi: required amplitude Phi^2 ~ c/(kappa xi) -> energy density of the wave ~ omega^2 Phi^2 ~ the geometry's own scale for omega ~ 1/a
out["C_classical_wave"] = "minimal coupling: T_kk = (k.nabla Phi)^2 >= 0 (NEC satisfied: cannot hold a throat open); non-minimal xi: needs Phi^2 ~ c/(kappa |xi|) i.e. Planckian amplitude, test-field approximation inconsistent"
json.dump(out, open(os.path.join(HERE, "phase3d1_wave.json"), "w"), indent=1, default=str)
