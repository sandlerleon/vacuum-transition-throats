# -*- coding: utf-8 -*-
"""Phase 3F.0: causal relativistic dissipative fluids (Israel-Stewart type) and matter-recycling feasibility.  kappa = 8 pi G = 1 in the numerical parts.

Part 0  Independent verification of the Phase 3E.1 flat-bounce statement: (i) the derived Einstein-scalar-Gauss-Bonnet equations are compared EXACTLY with the published Guo-Schwarz form
        (6H^2 = phid^2 + 2V + 24 xid H^3 ;  2 Hdot = -phid^2 + 4 xidd H^2 + 4 xid H (2 Hdot - H^2), action R/2 - (dphi)^2/2 - V - xi R_GB^2/2), (ii) the Noether (reparametrisation) identity holds, (iii) nondegeneracy at H = 0.
Part 1  Bulk viscosity in closed-universe collapse (GR + Lambda + dust + causal bulk pressure Pi, tau Pi' + Pi = -zeta Theta): Hdot at a turning point; critical Lambda for a bounce with and without viscosity.
Part 2  Shear viscosity in anisotropic (axisymmetric Bianchi I) collapse with causal shear stress: can viscosity damp the spacetime shear while staying causal and hydrodynamic?  Derived inequality (1 + w) c^2 K >= 3; numerical test.
Part 3  Null energy condition of the dissipative stress in the causal regime.
"""
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from sympy.calculus.euler import euler_equations

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}

# ------------------------------------------------------------------------------------------------ Part 0
t = sp.Symbol("t")
N = sp.Function("N")(t); a = sp.Function("a")(t); ph = sp.Function("phi")(t)
V = sp.Function("V"); xig = sp.Function("xig")                                       # Guo-Schwarz xi, coupling in the action:  - xig R_GB^2 / 2
ad = sp.diff(a, t); Nd = sp.diff(N, t)
Rsc = 6 * (sp.diff(a, t, 2) / (a * N ** 2) - ad * Nd / (a * N ** 3) + ad ** 2 / (a ** 2 * N ** 2))
sqrtgG = 24 * (ad ** 2 / N ** 2) * sp.diff(ad / N, t)
Lag = a ** 3 * N * (Rsc / 2 + sp.diff(ph, t) ** 2 / (2 * N ** 2) - V(ph)) - sp.Rational(1, 2) * xig(ph) * sqrtgG
EN, Ea, Ep = [e.lhs for e in euler_equations(Lag, [N, a, ph], t)]
# Noether identity with a general lapse: a' E_a + phi' E_phi - N d(E_N)/dt = 0 ?
noether = sp.simplify(sp.expand(sp.diff(a, t) * Ea + sp.diff(ph, t) * Ep - N * sp.diff(EN, t) - 0))
noether2 = sp.simplify(sp.expand(sp.diff(a, t) * Ea + sp.diff(ph, t) * Ep - sp.diff(N * EN, t)))
print("0: Noether identity residuals: a'E_a + phi'E_phi - N E_N' =", noether, " ; with d(N E_N)/dt:", noether2)
out["0_noether_residuals"] = [str(noether), str(noether2)]
def at1(e):
    for n_ in (3, 2, 1):
        e = e.subs(sp.diff(N, t, n_), 0)
    return e.subs(N, 1).doit()
ENs, Eas, Eps = [sp.simplify(at1(e)) for e in (EN, Ea, Ep)]
A, H, Hd, P, Pd, Vs, x1, x2 = sp.symbols("A H Hd Pd Pdd Vv xi1 xi2")
Pdot = sp.Symbol("Pdot")
def sym(e):
    e = e.subs(sp.diff(ph, t, 2), Pd).subs(sp.diff(ph, t), Pdot)
    e = e.subs(sp.diff(a, t, 2), A * (Hd + H ** 2)).subs(sp.diff(a, t), A * H).subs(a, A)
    return e
ENx, Eax, Epx = [sym(e) for e in (ENs, Eas, Eps)]
# replace the function values / derivatives by symbols
phis = sp.Symbol("ph0")
def funcs(e):
    e = e.replace(lambda x: isinstance(x, sp.Derivative) and x.expr.func in (V, xig), lambda x: x)
    return e
Vp = sp.Symbol("Vp"); Vsym = sp.Symbol("Vv")
Ens = ENx.subs({sp.Derivative(xig(ph), ph): x1}) if False else ENx
# explicit substitution of xi(phi) derivatives: xi(phi(t)) -> xig, d xi/dt = x1*phidot, d2 xi/dt2 = x2*phidot^2 + x1*phiddot
def repl(e):
    ph_sym = sp.Symbol("phi_val")
    e = e.subs(ph, ph_sym)
    e = e.doit()
    e = e.subs(sp.Derivative(xig(ph_sym), (ph_sym, 2)), x2).subs(sp.Derivative(xig(ph_sym), ph_sym), x1).subs(sp.Derivative(V(ph_sym), ph_sym), Vp).subs(V(ph_sym), Vs)
    return sp.simplify(e)
ENr, Ear, Epr = [repl(e) for e in (ENx, Eax, Epx)]
# Guo-Schwarz (kappa = 1):
xid = x1 * Pdot; xidd = x2 * Pdot ** 2 + x1 * Pd
GS_F = 6 * H ** 2 - Pdot ** 2 - 2 * Vs - 24 * xid * H ** 3
GS_R = 2 * Hd + Pdot ** 2 - 4 * xidd * H ** 2 - 4 * xid * H * (2 * Hd - H ** 2)
GS_S = Pd + 3 * H * Pdot + Vp + 12 * x1 * H ** 2 * (Hd + H ** 2)
c1, c2, c3 = sp.symbols("c1 c2 c3")
# my constraint E_N should be proportional to GS_F ; my a-equation to a combination of GS_R and GS_F ; my phi-equation to GS_S
rN = sp.simplify(ENr / GS_F)
print("0: E_N / (Guo-Schwarz Friedmann) =", rN)
sol = sp.solve(sp.Poly(sp.expand(Ear - c1 * GS_R - c2 * GS_F), Hd, Pd, Pdot, H, Vs, x1, x2, Vp).coeffs(), [c1, c2], dict=True)
print("0: a-equation = c1 * GS_R + c2 * GS_F with", sol)
rP = sp.simplify(Epr / GS_S)
print("0: E_phi / (Guo-Schwarz scalar equation) =", rP)
out["0_EN_over_GSF"] = str(rN); out["0_Ea_decomposition"] = [str(s_) for s_ in sol]; out["0_Ephi_over_GSS"] = str(rP)
# nondegeneracy: coefficient of Hdot in GS_R at H = 0
coef = sp.simplify(sp.diff(GS_R, Hd).subs(H, 0))
print("0: coefficient of Hdot at H = 0 in the Raychaudhuri equation:", coef, " -> Hdot(H=0) = -phid^2/2 exactly")
out["0_coeff_Hd_at_H0"] = str(coef)
Hd_solved = sp.solve(GS_R.subs(H, 0), Hd)[0]
out["0_Hdot_at_H0"] = str(Hd_solved)

# ------------------------------------------------------------------------------------------------ Part 1: bulk viscosity, closed universe + Lambda + dust (Lemaitre-type bounce)
KAP = 1.0; LAM = 0.02; A0 = 30.0
def run_bulk(M, zeta_frac, K, k=1, tmax=3000.0):
    kap = KAP; Lam = LAM
    rho0 = M / A0 ** 3
    HL = np.sqrt(kap * Lam / 3.0)
    tau = K / HL
    h2 = kap * (rho0 + Lam) / 3 - k / A0 ** 2
    if h2 <= 0: return None
    H0 = -np.sqrt(h2)
    def rhs(tt, y):
        A_, H_, rho, Pi = y
        zeta = zeta_frac * tau * rho                                                  # c_b^2 = zeta / (tau (rho + p)), dust c_s^2 = 0
        Hd_ = -0.5 * kap * (rho + Pi) + k / A_ ** 2
        return [H_ * A_, Hd_, -3 * H_ * (rho + Pi), (-Pi - 3 * zeta * H_) / tau]
    ev = lambda tt, y: y[1]; ev.direction = 1; ev.terminal = True
    evs = lambda tt, y: y[0] - 0.05; evs.terminal = True; evs.direction = -1
    sol = solve_ivp(rhs, (0, tmax), [A0, H0, rho0, 0.0], method="DOP853", rtol=1e-11, atol=1e-14, events=[ev, evs], max_step=0.5)
    y = sol.y
    res = np.max(np.abs(y[1] ** 2 - (kap * (y[2] + Lam) / 3 - k / y[0] ** 2)) * y[0] ** 2)
    return dict(bounce=bool(len(sol.t_events[0]) > 0), a_min=float(y[0].min()), Pi_over_rho_max=float(np.max(np.abs(y[3]) / np.maximum(y[2], 1e-30))), Pi_min=float(y[3].min()), residual_a2=float(res))
def critical_M(zeta_frac, K):
    lo, hi = 1e-4, 400.0
    r_lo = run_bulk(lo, zeta_frac, K)
    if r_lo is None or not r_lo["bounce"]: return None
    for _ in range(45):
        mid = np.sqrt(lo * hi)
        r_ = run_bulk(mid, zeta_frac, K)
        if r_ is not None and r_["bounce"]: lo = mid
        else: hi = mid
    return float(lo)
bulk = {}
print("1: closed universe with Lambda = %.2f, dust, a0 = %.0f (R_dS = %.2f): critical dust mass M_c below which the contracting universe bounces (Lemaitre); larger M_c = easier bounce" % (LAM, A0, np.sqrt(3 / (KAP * LAM))))
for K in (0.1, 0.3, 1.0):
    for cb in (0.0, 0.25, 1.0):
        Mc = critical_M(cb, K)
        bulk["K=%g c_b2=%g" % (K, cb)] = Mc
        print("1:   K = %.1f, c_b^2 = %.2f: M_c = %s" % (K, cb, ("%.5g" % Mc) if Mc else None))
ex = run_bulk(1.0, 1.0, 0.3)
out["1_critical_dust_mass"] = bulk; out["1_example_run"] = ex
print("1: example (M = 1, K = 0.3, c_b^2 = 1):", ex)
out["1_statement"] = ("Pi(t) = -zeta int Theta(t') e^{-(t-t')/tau} dt'/tau >= 0 after a contracting phase, so at a turning point Hdot = -(kappa/2)(rho + p + Pi) + k/a^2 is smaller than the inviscid value: "
                      "bulk viscosity can only raise the bounce threshold (lower M_c) and never produces the bounce.")

# ------------------------------------------------------------------------------------------------ Part 2: shear viscosity in axisymmetric Bianchi I collapse
def run_shear(w, K, c2, s0=0.05, Vend=1e-3, causal_only=True):
    kap = 1.0
    rho0 = 1.0
    h = -np.sqrt(rho0 / (3 * (1 - s0 ** 2)))
    y0 = [h * (1 + 2 * s0), h * (1 - s0), 0.0, rho0, 0.0]                                 # alpha_dot, beta_dot, S, rho, ln V
    def rhs(tt, y):
        ad, bd, S, rho, lnV = y
        Th = ad + 2 * bd
        p = w * rho
        tau = K / max(abs(Th), 1e-8)
        eta = c2 * tau * (rho + p)                                                       # shear-mode speed^2 = eta / (tau (rho + p)) = c2
        px = p + S; py = p - 0.5 * S
        bdd = 0.5 * (-px - 3 * bd ** 2)
        add = -py - bdd - ad ** 2 - bd ** 2 - ad * bd
        Sd = (-S - (4.0 / 3.0) * eta * (ad - bd)) / tau
        rd = -Th * (rho + p) - (ad - bd) * S
        return [add, bdd, Sd, rd, Th]
    ev = lambda tt, y: y[4] - np.log(Vend); ev.terminal = True; ev.direction = -1
    sol = solve_ivp(rhs, (0, 400), y0, method="DOP853", rtol=1e-11, atol=1e-14, events=[ev], max_step=0.05)
    ad, bd, S, rho, lnV = sol.y
    Th = ad + 2 * bd
    Sig2 = (ad - bd) ** 2 / np.maximum(Th ** 2, 1e-30)
    resid = np.max(np.abs(2 * ad * bd + bd ** 2 - rho))
    p = w * rho
    # NEC of the stress along the principal directions: rho + p + min eigenvalue of pi_ab (S, -S/2, -S/2)
    nec_ratio = (rho + p + np.minimum(S, -S / 2)) / (rho + p)
    return dict(Sigma2_start=float(Sig2[0]), Sigma2_end=float(Sig2[-1]), Sigma2_max=float(Sig2.max()), constraint_residual=float(resid), reached_Vend=bool(len(sol.t_events[0]) > 0),
                nec_ratio_min=float(nec_ratio.min()), S_over_rho_max=float(np.max(np.abs(S) / rho)))
shear = []
for w in (1 / 3, 1.0):
    for K in (0.1, 0.3, 1.0, 3.0, 10.0):
        for c2 in (0.25, 1.0):
            r_ = run_shear(w, K, c2)
            thr = 3 * (1 - w) / (4 * (1 + w))
            r_.update(w=w, K=K, c2=c2, threshold_c2K=thr, predicted_damping=bool(c2 * K > thr), damped=bool(r_["Sigma2_end"] < r_["Sigma2_start"]))
            shear.append(r_)
out["2_shear_scan"] = shear
base = {str(w): run_shear(w, 1.0, 0.0) for w in (1 / 3, 1.0)}
out["2_inviscid_baseline"] = base
print("2: inviscid baseline Sigma^2 start -> end:", {k_: (round(v["Sigma2_start"], 4), round(v["Sigma2_end"], 4)) for k_, v in base.items()})
print("2: shear viscosity in contraction: Sigma^2 = sigma^2/(Theta^2/3), start -> end ; predicted net damping of Sigma^2 (quasi-static, K <~ 1) iff c^2 K > 3(1-w)/(4(1+w))")
for r_ in shear:
    print("2:   w=%.2f K=%5.1f c^2=%.2f : Sigma^2 %.4f -> %.4f ; damped=%s ; predicted=%s ; min NEC ratio %.3f ; constraint resid %.1e" % (r_["w"], r_["K"], r_["c2"], r_["Sigma2_start"], r_["Sigma2_end"], r_["damped"], r_["predicted_damping"], r_["nec_ratio_min"], r_["constraint_residual"]))
adm = [r_ for r_ in shear if r_["K"] <= 1.0 and r_["c2"] <= 1.0]
out["2_summary_admissible_agreement"] = sum(r_["damped"] == r_["predicted_damping"] for r_ in adm)
out["2_summary"] = dict(admissible_runs=len(adm), admissible_damped=sum(r_["damped"] for r_ in adm), all_runs=len(shear), damped_all=sum(r_["damped"] for r_ in shear),
                        damped_outside_regime=[(r_["w"], r_["K"], r_["c2"]) for r_ in shear if r_["damped"] and r_["K"] > 1.0],
                        agree_with_prediction=sum(r_["damped"] == r_["predicted_damping"] for r_ in shear))
print("2: admissible (K <= 1, c^2 <= 1) runs: %d, damped %d ; agreement of the refined inequality in the admissible domain: %d / %d (all %d runs: %d agree; K>1 outside the quasi-static derivation)" %
      (len(adm), out["2_summary"]["admissible_damped"], out["2_summary_admissible_agreement"], len(adm), len(shear), out["2_summary"]["agree_with_prediction"]))
out["2_inequality"] = ("Exact shear equation in these variables: d(alpha'-beta')/dt = -Theta (alpha'-beta') + (3/2) S, S = -(4/3) eta (alpha'-beta') quasi-statically => sigma'/sigma = |Theta| - 2 kappa eta. With Theta'/Theta = (1+w)|Theta|/2 in a power-law contraction, "
                       "d ln Sigma^2/dt = (1-w)|Theta| - 4 kappa eta. Causality eta = c^2 tau (rho+p), Knudsen tau|Theta| = K, Theta^2 = 3 kappa rho: net damping of Sigma^2 iff c^2 K > 3(1-w)/(4(1+w)). "
                       "(A first estimate used kappa eta >= |Theta| and a threshold (1+w)c^2K >= 3; the numerical test contradicted it and the derivation was corrected.)")
json.dump(out, open(os.path.join(HERE, "phase3f0_fluid.json"), "w"), indent=1, default=str)
