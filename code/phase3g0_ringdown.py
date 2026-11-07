# -*- coding: utf-8 -*-
"""Phase 3G.0: damped oscillations after a high-energy transition ("nine diminishing bounces"): cosmological ringdown (Direction B) versus repeated universe-wide bounces (Direction A).  kappa = 8 pi G = 1.
The pasted proposal left section 6 ("Recommended Phase 3G.0") empty; this is the natural reading: (A) arithmetic of a damped sequence; (B) ringdown of a scalar/modulus in an expanding flat universe: is there any bounce, what damps it, what does the
amplitude sequence look like; (C) repeated closed-universe bounces using the Phase 1B model: do the bounces diminish, where does the energy go, how does it end; (D) scale constraints.
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp

import closed_cycle as CC
import closed_dissipation as CD

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}

# ------------------------------------------------------------------------------------------------ Part A
q = 0.7
gam_over_wd = -np.log(q) / (2 * np.pi)
A9 = q ** 8
print("A: peak ratio q = %.2f -> gamma/omega_d = %.5f ; amplitude of the 9th peak relative to the 1st = q^8 = %.4f (%.1f%%)" % (q, gam_over_wd, A9, 100 * A9))
out["A_gamma_over_omegad"] = float(gam_over_wd); out["A_ninth_peak"] = float(A9)
nthr = {str(th): float(1 + np.log(1 / th) / (-np.log(q))) for th in (0.1, 0.05, 0.01, 0.001)}
print("A: number of peaks before the amplitude falls below a threshold (q = 0.7):", {k: round(v, 1) for k, v in nthr.items()}, "-> nine is not selected; it follows from the threshold")
out["A_peaks_to_threshold"] = nthr

# ------------------------------------------------------------------------------------------------ Part B: ringdown in a flat expanding universe
m = 1.0
def ringdown(H0_over_m, Gam, ncyc=14):
    # flat universe, massive scalar (potential m^2 phi^2/2) + radiation from decay Gamma: phi'' + (3H + Gamma) phi' + m^2 phi = 0 ; rho_r' + 4H rho_r = Gamma phi'^2 ; H^2 = rho/3
    phi0 = np.sqrt(6.0) * H0_over_m * m / m                                         # rho = m^2 phi0^2/2 = 3 H0^2  ->  phi0 = sqrt(6) H0 / m
    def rhs(t, y):
        la, ph, pd, rr = y
        rho = 0.5 * pd ** 2 + 0.5 * m ** 2 * ph ** 2 + rr
        H = np.sqrt(rho / 3.0)
        return [H, pd, -(3 * H + Gam) * pd - m ** 2 * ph, -4 * H * rr + Gam * pd ** 2]
    ev = lambda t, y: y[2]; ev.direction = 0
    tmax = ncyc * 2 * np.pi / m
    sol = solve_ivp(rhs, (0, tmax), [0.0, phi0, 0.0, 0.0], method="DOP853", rtol=1e-11, atol=1e-14, events=[ev], max_step=0.02, dense_output=True)
    te = sol.t_events[0]; ye = sol.y_events[0]
    peaks = [(float(t_), float(abs(y_[1]))) for t_, y_ in zip(te, ye) if t_ > 1e-9]
    ts = np.linspace(0, tmax, 200001); Y = sol.sol(ts)
    ph, pd, rr = Y[1], Y[2], Y[3]
    rho = 0.5 * pd ** 2 + 0.5 * m ** 2 * ph ** 2 + rr
    H = np.sqrt(rho / 3.0)
    Hdot = np.gradient(H, ts)
    p = 0.5 * pd ** 2 - 0.5 * m ** 2 * ph ** 2 + rr / 3.0
    addot_over_a = -(rho + 3 * p) / 6.0                                              # kappa = 1
    acc = addot_over_a > 0
    episodes = int(np.sum(np.diff(acc.astype(int)) == 1))
    # relative H ripple per cycle
    T = 2 * np.pi / m
    ripple = []
    for c in range(1, ncyc):
        sel = (ts >= (c - 0.5) * T) & (ts < (c + 0.5) * T)
        ripple.append(float((H[sel].max() - H[sel].min()) / H[sel].mean())) if sel.any() else None
    amps = np.array([a_ for _, a_ in peaks])
    return dict(H0_over_m=H0_over_m, Gamma_over_m=Gam, peaks=peaks, ratio=[float(amps[i + 1] / amps[i]) for i in range(len(amps) - 1)], A9_over_A1=float(amps[8] / amps[0]) if len(amps) > 8 else None,
                H_min=float(H.min()), Hdot_max=float(Hdot.max()), H_zero_crossings=0 if H.min() > 0 else 1, accel_episodes=episodes, acceleration_sign_changes=int(np.sum(np.diff(acc.astype(int)) != 0)),
                H_fractional_change_per_cycle=ripple[:12], radiation_fraction_end=float(rr[-1] / rho[-1]))
H0m = -np.log(0.7) / (1.5 * np.pi)                                                    # Hubble friction: ratio between successive extrema of |phi| (half period pi/m) = exp(-1.5 pi H/m) = 0.7
cases = {}
for name, Gam in (("Hubble friction only", 0.0), ("Gamma/m = 0.2267 (extremum ratio 0.7 from decay alone)", -2 * np.log(0.7) / np.pi), ("Gamma/m = 0.03", 0.03)):
    r_ = ringdown(H0m, Gam)
    cases[name] = r_
    print("B: %-52s: peak ratios %s ; A9/A1 = %s ; H(t) zero crossings %d, max Hdot = %.2e ; acceleration episodes %d ; fractional change of H per cycle (first 5) %s" %
          (name, [round(x, 3) for x in r_["ratio"][:6]] + ["..."] + [round(x, 3) for x in r_["ratio"][-2:]], ("%.4f" % r_["A9_over_A1"]) if r_["A9_over_A1"] else "-", r_["H_zero_crossings"], r_["Hdot_max"], r_["accel_episodes"], [round(x, 4) for x in r_["H_fractional_change_per_cycle"][:5]]))
out["B_ringdown"] = cases
out["B_statement"] = ("Flat GR + scalar: Hdot = -(1/2) phidot^2 - ... <= 0 and rho > 0, so H never changes sign: the oscillation is a ripple in addot/a, not a bounce. Hubble friction alone gives amplitudes ~ a^(-3/2) ~ 1/t: peak ratios exp(-3 pi H/m) -> 1 as H/m falls (power law, not geometric 0.7^n); "
                      "a geometric sequence needs a constant damping rate (decay width): ratio exp(-pi Gamma/m).")

# ------------------------------------------------------------------------------------------------ Part C: repeated universe-wide bounces (Phase 1B closed model)
fam = dict(q=1.0, rc=0.2, mu=20.0, x=0.3); Gm = 0.3
mod = CC.Model(q=fam["q"], rc=fam["rc"], mu=fam["mu"], C=fam["x"] * CC.C_CRIT, ups=Gm)
best = None
for ls in np.linspace(-2.6, -0.3, 300):
    y0 = CD.start_state(mod, 10 ** ls * mod.v)
    if y0 is None: continue
    B, T, end, res, sol = CD.necks(mod, y0, 150.0, rtol=1e-11)
    if best is None or len(B) > best[0]:
        best = (len(B), float(ls), B, T, end, float(sol.t[-1]), res)
nB, ls, B, T, end, tend, res = best
aneck = [float(b[1][0]) for b in B]; aturn = [float(t_[1][0]) for t_ in T]
S = [float(b[1][4] * b[1][0] ** 4) for b in B]
radV = [float(b[1][4] / CC.VF) for b in B]
print("C: Phase 1B closed model (Gamma = %.1f): best seed 10^%.3f gives %d consecutive bounces, ends %s at t = %.1f ; turnaround radii a_max: %s ; neck radii a_min: %s ; entropy S per bounce: %s" %
      (Gm, ls, nB, end, tend, [round(x, 2) for x in aturn[:12]], [round(x, 3) for x in aneck[:12]], [round(x, 4) for x in S[:12]]))
out["C_closed_sequence"] = dict(n_bounces=nB, log10_seed=ls, end=end, t_end=tend, a_turnaround=aturn, a_neck=aneck, S=S, radiation_over_Vf=radV, max_constraint=res)
nn_ = min(len(aneck), len(aturn))
delta = [aturn[i] - aneck[i] for i in range(nn_)]
out["C_amplitude_a_turn_minus_a_neck"] = delta
print("C: oscillation amplitude a_turn - a_neck per bounce:", [round(x, 3) for x in delta[:12]])
dec = lambda xs: float(np.polyfit(np.arange(len(xs)), np.log(xs), 1)[0]) if len(xs) > 3 else None
out["C_trends"] = dict(amplitude_log_slope=(float(np.polyfit(np.arange(len(delta)), np.log(np.maximum(delta, 1e-12)), 1)[0]) if len(delta) > 3 else None), a_turn_log_slope=dec(aturn), a_neck_log_slope=dec(aneck), S_monotone_increasing=bool(all(S[i + 1] >= S[i] for i in range(len(S) - 1))) if len(S) > 1 else None)
print("C: trends: d ln(a_turn)/dn = %s ; d ln(a_neck)/dn = %s ; entropy S monotone increasing: %s" % (out["C_trends"]["a_turn_log_slope"], out["C_trends"]["a_neck_log_slope"], out["C_trends"]["S_monotone_increasing"]))
lt = json.load(open(os.path.join(HERE, "closed_dissipation_lifetime.json")))
out["C_lifetime_stats"] = [dict(gamma=r["gamma"], family=r["family"], max_bounces=r["max_bounces"], frac_ge10=r["frac_ge10"], end=r["best_end"]) for r in lt if r["family"]["mu"] == 20.0]
print("C: stored lifetime scan (family mu=20): max bounces by Gamma:", [(r["gamma"], r["max_bounces"], r["best_end"]) for r in lt if r["family"]["mu"] == 20.0])

# ------------------------------------------------------------------------------------------------ Part D: scale constraints
MPl = 1.221e19                                                                           # GeV
H_BBN = 1.66 * np.sqrt(10.75) * (1e-3) ** 2 / MPl                                         # at T = 1 MeV
hbar = 6.582e-25                                                                         # GeV s
out["D_scale"] = dict(H_at_T_1MeV_GeV=float(H_BBN), t_BBN_s=float(hbar / (2 * H_BBN)), Gamma_decay_must_exceed_GeV=float(H_BBN),
                      note="the ringing field must have decayed (Gamma >~ H) before BBN; for any mass m the number of remaining oscillations before decay is ~ m / (pi Gamma) ln(A/A_th)")
for mm in (1e-6, 1e3, 1e13):                                                              # GeV
    out["D_scale"]["oscillations_to_threshold_m_%g_GeV" % mm] = float(mm / (np.pi * H_BBN) * np.log(1 / 0.05))
print("D: H(T = 1 MeV) = %.3g GeV (t ~ %.2f s): a ringing field decaying only at Gamma ~ H_BBN completes N ~ m/(pi Gamma) ln(1/0.05) = %.2g oscillations for m = 1 TeV: the number of oscillations is set by m/Gamma, not by any fixed nine" %
      (H_BBN, hbar / (2 * H_BBN), out["D_scale"]["oscillations_to_threshold_m_1000_GeV"]))
json.dump(out, open(os.path.join(HERE, "phase3g0_ringdown.json"), "w"), indent=1, default=str)
