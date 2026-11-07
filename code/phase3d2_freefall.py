# -*- coding: utf-8 -*-
"""Phase 3D.2: free-fall transport through the Simpson-Visser black-bounce benchmark  ds^2 = -f dt^2 + dr^2/f + R^2 dOmega^2, R^2 = r^2 + a^2, f = 1 - 2m/R.   G = c = 1, kappa = 8 pi.
The geometry is PRESCRIBED (Phase 3D.1); nothing here shows that it forms dynamically.

Part A  geodesics: (dr/dtau)^2 = E^2 - U(r),  U = f (1 + L^2/R^2),  rddot = -U'/2 (regular across horizons and turning points); crossing criterion, direct integration, conserved quantities.
Part B  ISS-style comparison: circular geodesics, throat orbit stability, critical angular momentum, whether L helps.
Part C  tidal tensor in the freely falling frame (symbolic), Jacobi fields (transverse solution must equal R(tau)), maximum tidal strength.
Part D  causal structure: region classification along the worldline, null expansions, anti-trapped region, congruence expansion Theta, one-way property.
Part E  stress-energy required, Penrose singularity theorem consistency, what free fall does not change.
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

# ------------------------------------------------------------------------------------------------ symbolic setup
t, r, th, ph = sp.symbols("t r theta phi", real=True)
a, m, E, L = sp.symbols("a m E L", positive=True)
Rs = sp.sqrt(r ** 2 + a ** 2)
fs = 1 - 2 * m / Rs
X = [t, r, th, ph]
g = sp.diag(-fs, 1 / fs, Rs ** 2, Rs ** 2 * sp.sin(th) ** 2)
gi = g.inv()
Gam = [[[sum(gi[i, l] * (sp.diff(g[l, j], X[k]) + sp.diff(g[l, k], X[j]) - sp.diff(g[j, k], X[l])) for l in range(4)) / 2 for k in range(4)] for j in range(4)] for i in range(4)]
def Riem_up(i, j, k, l):
    e = sp.diff(Gam[i][j][l], X[k]) - sp.diff(Gam[i][j][k], X[l])
    return e + sum(Gam[i][k][s] * Gam[s][j][l] - Gam[i][l][s] * Gam[s][j][k] for s in range(4))
Rl = lambda i, j, k, l: sp.simplify(g[i, i] * Riem_up(i, j, k, l))                 # R_{ijkl} (diagonal metric)
R_trtr = Rl(0, 1, 0, 1); R_thtth = Rl(2, 0, 2, 0); R_thrth = Rl(2, 1, 2, 1)
rd2 = E ** 2 - fs                                                                  # rdot^2 on a radial geodesic
K11 = sp.simplify(R_trtr * ((rd2 - E ** 2) / fs) ** 2)                            # R_{e1 u e1 u} = R_trtr (e1^t u^r - e1^r u^t)^2 ,  (rd^2 - E^2)/f = -1
Kth = sp.simplify((R_thtth * E ** 2 / fs ** 2 + R_thrth * rd2) / Rs ** 2)          # R_{e_th u e_th u}
print("C: radial tidal K_11 =", K11, "\n   transverse tidal K_th =", Kth)
out["C_K11"] = str(K11); out["C_Kth"] = str(Kth)
K11f = sp.lambdify((r, a, m, E), K11, "numpy"); Kthf = sp.lambdify((r, a, m, E), Kth, "numpy")
# regularity at the horizon: limit of K as f -> 0 (E fixed)
rh = sp.sqrt(4 * m ** 2 - a ** 2)
lim_th = sp.simplify(sp.limit(Kth.subs({a: 1, m: 1}), r, sp.sqrt(3)))
print("C: transverse tidal at the horizon r_h (a=m=1) finite:", lim_th)
out["C_Kth_at_horizon_a1m1"] = str(lim_th)

# ------------------------------------------------------------------------------------------------ geodesic machinery
Ufun = lambda rr, aa, mm, LL: (1 - 2 * mm / np.sqrt(rr ** 2 + aa ** 2)) * (1 + LL ** 2 / (rr ** 2 + aa ** 2))
Us = fs * (1 + L ** 2 / Rs ** 2)
dUf = sp.lambdify((r, a, m, L), sp.diff(Us, r), "numpy")
def geodesic(aa, mm, LL, r0, rd0, tau_max, stop=None, events=(), dense=True):
    f_ = lambda y: dUf(y[0], aa, mm, LL)
    def rhs(tau, y):
        return [y[1], -0.5 * dUf(y[0], aa, mm, LL), LL / (y[0] ** 2 + aa ** 2)]
    return solve_ivp(rhs, (0, tau_max), [r0, rd0, 0.0], method="DOP853", rtol=1e-12, atol=1e-13, events=list(events), dense_output=dense, max_step=0.5)
def Rnorm(rr, aa): return np.sqrt(rr ** 2 + aa ** 2)

# ---- Part A: radial infall from rest at R0 (bound, E^2 = f(R0)): a/m = 1 (black bounce) and a/m = 3 (wormhole)
Apart = {}
for name, aa, mm, R0 in (("black bounce a/m=1, rest at R0=20m", 1.0, 1.0, 20.0), ("wormhole a/m=3, rest at R0=20m", 3.0, 1.0, 20.0)):
    r0 = np.sqrt(R0 ** 2 - aa ** 2)
    E2 = 1 - 2 * mm / R0
    ev_zero = lambda tau, y: y[0]; ev_zero.terminal = False
    rh_ = np.sqrt(max(4 * mm ** 2 - aa ** 2, 0.0))
    evs = [ev_zero]
    if rh_ > 0:
        e_h1 = lambda tau, y: y[0] - rh_; e_h2 = lambda tau, y: y[0] + rh_
        evs += [e_h1, e_h2]
    tmax = 3000.0
    sol = geodesic(aa, mm, 0.0, r0, 0.0, tmax, events=evs)
    tau = np.linspace(0, tmax, 200001)
    y = sol.sol(tau)
    cons = np.max(np.abs(y[1] ** 2 + Ufun(y[0], aa, mm, 0.0) - E2))
    turn = [(float(tau[i]), float(y[0][i])) for i in range(1, len(tau) - 1) if (y[1][i - 1] * y[1][i + 1] < 0 and abs(y[1][i]) < 5e-3 and abs(y[0][i]) > 1)]
    cross0 = [float(x) for x in sol.t_events[0]]
    # collapse repeated turning detections
    turns = []
    for tt, rr in turn:
        if not turns or abs(tt - turns[-1][0]) > 20:
            turns.append((tt, rr))
    info = dict(E2=float(E2), conserved_E_drift=float(cons), r_zero_crossings=cross0[:6], turning_points_tau_r=turns[:6],
                R_min_along=float(Rnorm(y[0], aa).min()), tau_first_zero=cross0[0] if cross0 else None)
    if rh_ > 0:
        info["tau_cross_horizon_plus"] = float(sol.t_events[1][0]) if len(sol.t_events[1]) else None
        info["tau_cross_horizon_minus"] = [float(x) for x in sol.t_events[2][:3]]
    Apart[name] = info
    print("A: %s: E^2 = %.4f, drift of E^2 = %.2e, first r=0 at tau = %.2f, turning points (tau, r) = %s" % (name, E2, cons, info["tau_first_zero"], [(round(tt, 1), round(rr, 2)) for tt, rr in turns[:4]]))
out["A_radial"] = Apart
# E = 1 (released from rest at infinity): reaches r = -infinity in finite proper time? proper time from r=+50 to r=-50 (3D.1) ~ 335; no turning point
# ---- crossing criterion: crossing iff E^2 > max_{R in [a, R0]} U(R); critical L; compare with Schwarzschild
def Umax(aa, mm, LL, R0=100.0):
    Rg = np.concatenate([np.linspace(aa, 12 * mm + 4 * aa, 40001), np.geomspace(12 * mm + 4 * aa, R0, 4000)])
    return float(np.max((1 - 2 * mm / Rg) * (1 + LL ** 2 / Rg ** 2)))
crit = {}
for aa in (0.5, 1.0, 1.9, 3.0, 6.0):
    for E_ in (1.0, 1.2):
        Lc = brentq(lambda LL: Umax(aa, 1.0, LL) - E_ ** 2, 0.0, 60.0, xtol=1e-10) if Umax(aa, 1.0, 0.0) < E_ ** 2 else 0.0
        # direct check just below / above
        res = {}
        for fac in (0.98, 1.02):
            LL = Lc * fac
            R0 = 100.0
            r0 = np.sqrt(R0 ** 2 - aa ** 2)
            ev_turn = lambda tau, y: y[1]; ev_turn.terminal = True; ev_turn.direction = 1
            ev_gone = lambda tau, y: y[0] + 3.0 * aa - 0.0; ev_gone.terminal = True; ev_gone.direction = -1
            sol = geodesic(aa, 1.0, LL, r0, -np.sqrt(max(E_ ** 2 - Ufun(r0, aa, 1.0, LL), 1e-30)), 4000.0, events=[ev_turn, ev_gone], dense=False)
            res["L=%.2f Lc" % fac] = "crosses throat" if (len(sol.t_events[1]) > 0 and len(sol.t_events[0]) == 0) else "scattered back"
        crit["a/m=%.1f E=%.2f" % (aa, E_)] = dict(L_critical=float(Lc), direct_integration=res)
        print("A/B: a/m = %.1f, E = %.2f: critical L for crossing = %.4f (Schwarzschild plunge at E=1: 4.0000), direct: %s" % (aa, E_, Lc, res))
out["B_critical_L"] = crit

# ---- Part B: circular geodesics and throat orbit
Rg = np.linspace(0.01, 60, 600001)
circ = {}
for aa in (0.5, 1.0, 1.9, 3.0, 6.0):
    Ls = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0])
    rows = []
    for LL in Ls:
        Rr = np.linspace(max(aa, 1e-3), 80, 400001)
        U = (1 - 2 / Rr) * (1 + LL ** 2 / Rr ** 2)
        dU = np.gradient(U, Rr)
        idx = np.where(np.sign(dU[:-1]) != np.sign(dU[1:]))[0]
        for i in idx:
            rr = float(Rr[i])
            if rr > aa + 1e-3 and rr > 2.0:                                          # exterior orbits (r>0 side with f>0)
                d2 = (dU[i + 1] - dU[i - 1]) / (Rr[i + 1] - Rr[i - 1])
                rows.append(dict(L=float(LL), R=rr, stable=bool(d2 > 0)))
    # throat orbit at r = 0 (exists as a static region only for a > 2m): d^2U/dr^2 = U_R(a)/a, U_R = 2m/R^2 + (L^2/R^3)(6m/R - 2)
    throat = {}
    if aa > 2.0:
        for LL in (0.0, 2.0, 6.0, 12.0):
            UR = 2 / aa ** 2 + LL ** 2 / aa ** 3 * (6 / aa - 2)
            throat["L=%g" % LL] = "stable (potential minimum)" if UR > 0 else "unstable (barrier at throat)"
        throat["L_threshold"] = float(np.sqrt(aa ** 2 / (aa - 3))) if aa > 3 else None
    circ["a/m=%.1f" % aa] = dict(exterior_circular=rows, throat=throat)
    nst = sum(1 for rw in rows if rw["stable"]); nun = sum(1 for rw in rows if not rw["stable"])
    print("B: a/m = %.1f: exterior circular geodesics among L in {1..8}: %d stable, %d unstable; throat orbit: %s" % (aa, nst, nun, throat))
out["B_circular"] = circ
# monotonicity: the barrier U increases with L^2 wherever f>0 ; so more angular momentum can only remove crossing orbits
out["B_L_never_helps"] = "dU/d(L^2) = f / R^2 > 0 outside the horizon, so the allowed-crossing set is {L <= L_c(E)}; inside the horizon f < 0 there is no turning point"

# ------------------------------------------------------------------------------------------------ Part C: Jacobi fields along the radial geodesic (black bounce a/m=1, E^2 = f(R0=20))
aa, mm, R0 = 1.0, 1.0, 20.0
E_ = np.sqrt(1 - 2 * mm / R0)
r0 = np.sqrt(R0 ** 2 - aa ** 2)
def rhs_full(tau, y):
    rr, rd, xr, xrd, xt, xtd = y
    return [rd, -0.5 * dUf(rr, aa, mm, 0.0), xrd, -K11f(rr, aa, mm, E_) * xr, xtd, -Kthf(rr, aa, mm, E_) * xt]
# the geodesic must be integrated in tau with rd = -sqrt(E^2 - f); the tidal expressions depend on r and E only
y0 = [r0, 0.0, 1e-3, 0.0, 1e-3 * Rnorm(r0, aa) / Rnorm(r0, aa), 0.0]
# transverse comparison solution: xi_th = R(tau) * c  (neighbouring radial geodesics of the same shell); initial xi_th(0) = R0 * eps, xi_th'(0) = R'(0) rdot * eps = 0
y0 = [r0, 0.0, 1e-3, 0.0, 1e-3 * R0, 0.0]
sol = solve_ivp(rhs_full, (0, 300.0), y0, method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True, max_step=0.2)
tt = np.linspace(0, 300.0, 60001); yy = sol.sol(tt)
Rtau = Rnorm(yy[0], aa)
ratio_th = yy[4] / (1e-3 * Rtau)
dev = float(np.max(np.abs(ratio_th - 1.0)))
print("C: transverse Jacobi field vs R(tau): max |xi_th/(eps R) - 1| = %.2e (a shell of radial geodesics has proper size R(tau): finite, never zero, min R = %.3f)" % (dev, Rtau.min()))
Kr = K11f(yy[0], aa, mm, E_); Kt = Kthf(yy[0], aa, mm, E_)
xr_growth = float(np.max(np.abs(yy[2])) / 1e-3)
print("C: tidal eigenvalues along the fall (m=a=1): K_11 in [%.3f, %.3f], K_th in [%.3f, %.3f]; radial separation max growth = %.3g, min |xi_r| = %.3g" % (Kr.min(), Kr.max(), Kt.min(), Kt.max(), xr_growth, np.abs(yy[2]).min() / 1e-3))
zc_r = int(np.sum(np.diff(np.sign(yy[2])) != 0)); zc_t = int(np.sum(np.diff(np.sign(yy[4])) != 0))
print("C: sign changes of the radial / transverse deviation (conjugate points):", zc_r, zc_t)
out["C_jacobi"] = dict(zero_crossings_radial=zc_r, zero_crossings_transverse=zc_t, transverse_vs_R_max_dev=dev, K11_range=[float(Kr.min()), float(Kr.max())], Kth_range=[float(Kt.min()), float(Kt.max())], radial_max_growth=xr_growth,
                       radial_xi_min=float(np.abs(yy[2]).min() / 1e-3), R_min=float(Rtau.min()))
# scaling of the maximum tidal strength with a (E^2 = f(R0 = 20) ), in units of m
sc = {}
for aa_ in (0.25, 0.5, 1.0, 1.5, 1.9):
    rr = np.linspace(-20, 20, 40001)
    sc[str(aa_)] = dict(Kth_abs_max=float(np.max(np.abs(Kthf(rr, aa_, 1.0, E_)))), K11_abs_max=float(np.max(np.abs(K11f(rr, aa_, 1.0, E_)))), kretschmann_throat=None)
    print("C: a/m = %.2f: max |K_11| = %.3f, max |K_th| = %.3f (units 1/m^2) ; Schwarzschild at R = a would give 2m/a^3 = %.3f" % (aa_, sc[str(aa_)]["K11_abs_max"], sc[str(aa_)]["Kth_abs_max"], 2 / aa_ ** 3))
out["C_tidal_scaling"] = sc

# ------------------------------------------------------------------------------------------------ Part D: causal structure along the worldline (black bounce a/m = 1)
aa, mm = 1.0, 1.0
rh_ = np.sqrt(3.0)
fnum = lambda rr: 1 - 2 * mm / Rnorm(rr, aa)
sol = geodesic(aa, mm, 0.0, np.sqrt(R0 ** 2 - aa ** 2), 0.0, 400.0, dense=True)
tau = np.linspace(0, 400, 400001); yy = sol.sol(tau)
rr, rd = yy[0], yy[1]
fr = fnum(rr); Rdot = rr * rd / Rnorm(rr, aa)
labels = np.empty(len(tau), dtype=object)
for i in range(len(tau)):
    if fr[i] > 0:
        labels[i] = "exterior (r>0)" if rr[i] > 0 else "exterior (r<0)"
    elif fr[i] < 0:
        # both radial null expansions have the sign of Rdot (r timelike) inside the horizon
        labels[i] = "trapped (R decreasing)" if Rdot[i] < 0 else "anti-trapped (R increasing)"
    else:
        labels[i] = "horizon"
seq = []
for i, lb in enumerate(labels):
    if not seq or seq[-1][0] != lb:
        seq.append((lb, float(tau[i])))
print("D: regions visited along the radial fall (tau):", [(lb, round(t_, 1)) for lb, t_ in seq])
out["D_region_sequence"] = seq
# congruence expansion Theta = (1/R^2) d(R^2 rd)/dr = 2 r rd/R^2 + d rd/dr ,  rd = -sqrt(E^2 - f)
def Theta(rr, Ee):
    f_ = fnum(rr); rd_ = -np.sqrt(np.maximum(Ee ** 2 - f_, 0)); fp = 2 * mm * rr / Rnorm(rr, aa) ** 3
    return 2 * rr * rd_ / Rnorm(rr, aa) ** 2 + fp / (2 * rd_)
xs = np.linspace(-rh_ + 1e-3, -1e-3, 2000)                                            # anti-trapped part of the interior
Th_a = Theta(xs, E_); Th_b = Theta(xs, 0.0 + 1e-9); Th_c = Theta(xs, 1.0)
print("D: volume expansion Theta of the radial congruence in the anti-trapped region (-r_h < r < 0): E^2=f(20): [%.3f, %.3f]; E=0: [%.3f, %.3f]; E=1: [%.3f, %.3f]" %
      (Th_a.min(), Th_a.max(), Th_b.min(), Th_b.max(), Th_c.min(), Th_c.max()))
out["D_theta_antitrapped"] = dict(E2_f20=[float(Th_a.min()), float(Th_a.max())], E0=[float(Th_b.min()), float(Th_b.max())], E1=[float(Th_c.min()), float(Th_c.max())])
# bound geodesic beyond the white-hole exit: turning point in the second exterior and the return
sol2 = geodesic(aa, mm, 0.0, np.sqrt(R0 ** 2 - aa ** 2), 0.0, 400.0, dense=True)
yy2 = sol2.sol(np.linspace(0, 400, 400001))
out["D_second_exterior"] = dict(r_min_reached=float(yy2[0].min()), R_at_turn=float(Rnorm(yy2[0].min(), aa)), returns_to_horizon_branch=bool(yy2[0][-1] > -rh_))
print("D: bound fall reaches r_min = %.3f (R = %.3f) in the second exterior, then returns toward the horizon branch (next trapped region): %s" % (yy2[0].min(), Rnorm(yy2[0].min(), aa), out["D_second_exterior"]["returns_to_horizon_branch"]))
out["D_one_way"] = "inside |r| < r_h the coordinate r is timelike and dr/dtau has a fixed sign for future-directed curves (U < 0 there), so no future-directed curve returns to the first exterior"

# ------------------------------------------------------------------------------------------------ Part E: stress-energy (3D.1 closed form) and Penrose theorem
Tkk = -2 * a ** 2 / (kappa * Rs ** 4)
integral = sp.integrate(Tkk * Rs ** 2, (r, -sp.oo, sp.oo))
print("E: T_kk = -2 a^2/(kappa R^4) < 0 at every r (independent of E, L, m); integral T_kk R^2 dr =", sp.simplify(integral), "  (= -a/4 in G=1 units)")
out["E_Tkk_integral"] = str(sp.simplify(integral))
out["E_penrose"] = ("The maximal extension contains closed trapped surfaces (r > 0 inside the horizon) and is geodesically complete. The Penrose theorem (NEC + closed trapped surface + non-compact Cauchy surface => null geodesic incompleteness) "
                    "then requires NEC violation, consistent with T_kk < 0. Formation from NEC-satisfying (classical ordinary) collapsing matter is excluded by the same theorem; any formation channel must violate the NEC (quantum/torsion/unknown).")
json.dump(out, open(os.path.join(HERE, "phase3d2_freefall.json"), "w"), indent=1, default=str)
