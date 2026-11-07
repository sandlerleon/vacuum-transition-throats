# -*- coding: utf-8 -*-
"""Phase 3C.2: symbolic verification of the kinematic identities used, built on the Einstein-Cartan 1+3 equations of Pasmatsiou, Tsagas & Barrow (arXiv:1611.07878, Eqs. 62-76, as transcribed in the
session; conventions: signature (-,+,+,+), kappa = 8 pi G, sigma^2 = sigma_ab sigma^ab/2, omega^2 = omega_ab omega^ab/2, s^2 = s_ab s^ab/2).

(1) effective fluid:  rho_eff = rho - kappa s^2/4, p_eff = p - kappa s^2/4  =>  -(kappa/2)(rho_eff + 3 p_eff) = -(kappa/2)(rho + 3p) + kappa^2 s^2 / 2   (the spin term of Eq. 76);
(2) with aligned spin s_ab = s * eps_abc n^c and vorticity omega~_ab = omega * eps_abc n^c:  2 omega~^2 - kappa s_ab omega~^ab + kappa^2 s^2/2 = 2 (omega~ - kappa s/2)^2  (spin and vorticity enter as a perfect square, consistent with
    omega_ab = omega~_ab - kappa s_ab/2, Eq. 74-75); the cross term vanishes for random orientation (<s_ab> = 0) leaving 2 omega~^2 + kappa^2 <s^2>/2 >= 0;
(3) circulation law from the vorticity equation (Eq. 64, geodesic, homogeneous limit): omega_dot = -(2/3) theta omega + 2 sigma_c[a omega_b]^c with axisymmetric shear gives omega ~ a_perp^{-2};
(4) power-law scalings of every term entering the Hamiltonian constraint, theta^2/3 = kappa rho + sigma^2 - omega^2 - (1/2) R3 + Lambda, and the dominance order as the volume V -> 0;
(5) scaling fit of the critical spin strength found in phase3c2_bianchi.py.
"""
import json
import os

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
kap, s, om, rho, p = sp.symbols("kappa s omega rho p", real=True)

# (1) effective-fluid identity
lhs = -(kap / 2) * ((rho - kap * s ** 2 / 4) + 3 * (p - kap * s ** 2 / 4))
rhs = -(kap / 2) * (rho + 3 * p) + kap ** 2 * s ** 2 / 2
out["effective_fluid_identity"] = str(sp.simplify(lhs - rhs))
print("(1) -(k/2)(rho_eff + 3 p_eff) - [-(k/2)(rho + 3p) + k^2 s^2/2] =", out["effective_fluid_identity"])

# (2) perfect square.  aligned: s_ab = s eps_abz, omega_ab = om eps_abz  ->  s_ab s^ab / 2 = s^2, s_ab omega^ab = 2 s om
alpha = 2 * om ** 2 - kap * (2 * s * om) + kap ** 2 * s ** 2 / 2
sq = 2 * (om - kap * s / 2) ** 2
out["perfect_square_difference"] = str(sp.simplify(sp.expand(alpha - sq)))
print("(2) 2 om^2 - kappa s_ab om^ab + kappa^2 s^2/2 - 2 (om - kappa s/2)^2 =", out["perfect_square_difference"])
# vanishing combination: corotating spin: om = kappa s / 2  -> total kinematic vorticity zero; the spin/rotation repulsion cancels
out["cancellation_at"] = "omega_tilde = kappa s / 2"
# random orientation: cross term averages to zero
out["random_orientation"] = "2 omega~^2 + kappa^2 <s^2>/2 (both >= 0)"

# (3) circulation law: omega_xy dot = -(2/3) theta omega + (sigma_{xc} omega_y^c - sigma_{yc} omega_x^c);  axisymmetric shear diag(s_perp, s_perp, s_par), s_par = -2 s_perp
th, sp_, H = sp.symbols("theta s_perp H_perp")
w_xy = sp.Symbol("w")
Om = sp.Matrix([[0, w_xy, 0], [-w_xy, 0, 0], [0, 0, 0]])                # omega_ab (spatial, lowered = raised for flat spatial metric)
Sg = sp.diag(sp_, sp_, -2 * sp_)
term = (Sg * Om.T - Om * Sg)[0, 1]                                      # sigma_{xc} omega_{y}^{c} - sigma_{yc} omega_{x}^{c} : (Sg * Om^T)[x,y] = sigma_xc omega_yc ; (Om * Sg)... check sign by explicit sum
explicit = sum(Sg[0, c] * Om[1, c] for c in range(3)) - sum(Sg[1, c] * Om[0, c] for c in range(3))
wdot = -sp.Rational(2, 3) * th * w_xy + explicit
wdot = sp.simplify(wdot.subs(sp_, H - th / 3))
out["omega_dot_over_omega"] = str(sp.simplify(wdot / w_xy))
print("(3) omega_dot/omega =", out["omega_dot_over_omega"], " -> omega ~ a_perp^{-2} (H_perp = a_perp_dot / a_perp)")

# (4) scalings in the contraction limit: exponent p with term ~ V^{-p}, isotropic estimate a ~ V^{1/3}, a_perp ~ V^{1/3}
V = sp.Symbol("V", positive=True)
terms = {"dust rho": 1, "curvature R3/2": sp.Rational(2, 3), "vorticity omega^2 (omega ~ a^-2)": sp.Rational(4, 3), "magnetic rho_B (flux, B ~ a^-2)": sp.Rational(4, 3),
         "orbital kinetic rho_r ~ n v^2, v ~ 1/a": sp.Rational(5, 3), "shear sigma^2 ~ V^-2": 2, "spin -B n^2": 2}
out["scaling_exponents_V^-p"] = {k: str(v) for k, v in terms.items()}
print("(4) exponents p of V^{-p}:", {k: str(v) for k, v in terms.items()})
out["dominance"] = "As V -> 0 only the V^-2 terms (spin, negative; shear, positive) matter: bounce iff the net V^-2 coefficient is negative: B N0^2 > Sigma_asymptotic (dust: V^-1, vorticity/magnetic/orbital: V^-5/3..-4/3 are subdominant)"

# (5) fit of the critical spin strength  B_c ~ c * X^q
Bj = json.load(open(os.path.join(HERE, "phase3c2_bianchi.json")))
rows = {r["case"].strip(): r for r in Bj["critical"]}
def fit(names, xs):
    y = [rows[n]["B_critical"] for n in names]
    q, lc = np.polyfit(np.log(xs), np.log(y), 1)
    return dict(exponent=float(q), coeff=float(np.exp(lc)), points=list(zip(xs, y)))
rot = fit(["C  rotating EC v0^2 = 0.01", "C  rotating EC v0^2 = 0.05", "C  rotating EC v0^2 = 0.2"], [0.01, 0.05, 0.2])
mag = fit(["D  magnetized nonrotating, rho_B0 = 0.01", "D  magnetized nonrotating, rho_B0 = 0.05", "D  magnetized nonrotating, rho_B0 = 0.2"], [0.01, 0.05, 0.2])
out["fit_rotation"] = rot; out["fit_magnetic"] = mag
print("(5) B_c ~ %.3g (v0^2)^%.2f  (rotation proxy);  B_c ~ %.3g (rho_B0)^%.2f  (magnetic)" % (rot["coeff"], rot["exponent"], mag["coeff"], mag["exponent"]))
# perpendicular magnetosonic speed at n* (w -> 0): x = [b^2 + n p'] / (w + b^2) stays finite
x_perp = (sp.Symbol("b2") + sp.Symbol("n_") * (-2 * sp.Symbol("B_") * sp.Symbol("n_"))) / (sp.Symbol("w_") + sp.Symbol("b2"))
out["perp_fast_speed2_at_nstar"] = str(sp.simplify(x_perp.subs(sp.Symbol("w_"), 0)))
print("perpendicular fast speed^2 at n* (w = 0): (b^2 - 2 B n^2)/b^2 =", out["perp_fast_speed2_at_nstar"], "; real and subluminal iff b^2 > 2 B n^2 (= 2 |p|)")
json.dump(out, open(os.path.join(HERE, "phase3c2_kinematics.json"), "w"), indent=1)
