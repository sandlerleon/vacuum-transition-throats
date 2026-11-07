# -*- coding: utf-8 -*-
"""Phase 3C.4: polarized spin in the Einstein-Cartan (Weyssenhoff) fluid.  Source: Pasmatsiou, Tsagas & Barrow, arXiv:1611.07878 (full text read in this session).  Conventions as in the paper:
signature (-,+,+,+), kappa = 8 pi G, s^2 = s_ab s^ab/2, Frenkel condition s_ab u^b = 0, effective stress (their Eq. 77)
    T_ab = (rho - kappa s^2/4) u_a u_b + (p - kappa s^2/4) h_ab - (1/2) h^{cd} nabla_d ( s_ca u_b + s_cb u_a )
(the last term is the part that only exists for polarized spin; the equations were NOT given in Riemannian shear/vorticity form in the paper, so the cross terms are derived here from Eq. 77).

Part A (sympy): cross-term stress for s_ab, omega_ab, sigma_ab along a common axis; fix the sign convention of omega by reproducing the paper's Raychaudhuri cross term -kappa s_ab omega^ab (Eq. 76).
Part B: stress-energy conservation check of the cross term in the homogeneous axisymmetric (Bianchi I) class, and the vorticity equation (their Eq. 64) for a Weyssenhoff fluid.
Part C: what-if bound: if a fraction f of the rotation-induced anisotropic stress could be cancelled, how far does the critical spin strength fall (phase3c2_bianchi model).
"""
import json
import os

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
kap, s, om, sg, th = sp.symbols("kappa s omega sigma_perp Theta", real=True)
eta = sp.Symbol("eta")                         # sign convention: D_d u_b = (Theta/3) h_db + sigma_db + eta * omega_bd
h3 = sp.eye(3)
S_ab = sp.Matrix([[0, s, 0], [-s, 0, 0], [0, 0, 0]])                         # s_ab (spatial), spin along z
W_ab = sp.Matrix([[0, om, 0], [-om, 0, 0], [0, 0, 0]])                       # omega_ab, vorticity along z
Sig = sp.diag(sg, sg, -2 * sg)                                               # axisymmetric shear
Du = lambda e: th / 3 * h3 + Sig + e * W_ab.T                                # (Du)_{db} with omega_{bd} = (W^T)_{db}
Du_e = Du(eta)

# effective stress cross term: T^(c)_ab = -(1/2) h^{cd} [ s_ca D_d u_b + s_cb D_d u_a ]  (spatial part) ;  uu part from the Frenkel derivative
def spatial_cross(Du_):
    # -(1/2) sum_{c,d} s_{ca} (Du)_{db} delta^{cd}  -> -(1/2) (S^T Du)_{ab}, plus the a<->b partner
    M = -sp.Rational(1, 2) * (S_ab.T * Du_)
    return M + M.T
Sab = sp.simplify(spatial_cross(Du_e))
# uu component: T_ab u^a u^b = - s_ca D^c u^a  (uses (nabla_d s_cb) u^b = - s_cb nabla_d u^b from the Frenkel condition)
Tuu = sp.simplify(-sum(S_ab[c, a] * Du_e[c, a] for c in range(3) for a in range(3)))
trace_spatial = sp.simplify(Sab.trace())
T_trace = sp.simplify(-Tuu + trace_spatial)                                  # g^{ab} T_ab = -T_uu + h^{ab} T_ab
Ray_cross = sp.simplify(-kap * Tuu - kap * T_trace / 2)                      # -kappa T_(ab) u^a u^b - (kappa/2) T
paper = -kap * sum(S_ab[a, b] * W_ab[a, b] for a in range(3) for b in range(3))   # -kappa s_ab omega^ab
sol_eta = sp.solve(sp.Eq(Ray_cross, paper), eta)
print("Raychaudhuri cross term: -kappa T_uu - kappa T/2 =", Ray_cross, " ; paper: -kappa s_ab omega^ab =", paper, " -> eta =", sol_eta)
out["eta_fixed_by_paper"] = [str(x) for x in sol_eta]
eta_val = sol_eta[0] if sol_eta else -1
Sab_v = sp.simplify(Sab.subs(eta, eta_val)); Tuu_v = sp.simplify(Tuu.subs(eta, eta_val)); Ttr = sp.simplify(T_trace.subs(eta, eta_val))
print("with eta =", eta_val, ": spatial cross stress =", Sab_v.tolist(), "; T_uu =", Tuu_v, "; trace T =", Ttr)
out["cross_stress"] = dict(spatial=str(Sab_v.tolist()), T_uu=str(Tuu_v), trace=str(Ttr))
# shear-propagation source (standard 1+3, no torsion once the effective stress is used):  sigma_dot_<ab> = ... - omega_c<a omega_b>^c + (kappa/2) pi_ab ,  pi = traceless spatial stress
pi = sp.simplify(Sab_v - sp.Rational(1, 3) * Sab_v.trace() * h3)
rot_src = sp.Matrix(3, 3, lambda a, b: -sum(W_ab[c, a] * W_ab[b, c] for c in range(3)))   # - omega_{ca} omega_{bc} (spatial indices raised with delta)
rot_src_tl = sp.simplify(rot_src - sp.Rational(1, 3) * rot_src.trace() * h3)
shear_src = sp.simplify(rot_src_tl + kap / 2 * pi)
print("rotation shear source (traceless):", rot_src_tl.tolist())
print("cross-term anisotropic stress pi_ab:", pi.tolist())
print("combined shear source:", shear_src.tolist())
# factor
fac = sp.factor(shear_src[0, 0])
out["shear_source_xx"] = str(fac)
print("combined shear source (xx) factorised:", fac, " -> proportional to omega (omega + kappa s / 2) up to sign conventions")
out["raychaudhuri_perfect_square"] = str(sp.factor(2 * om ** 2 - kap * 2 * s * om + kap ** 2 * s ** 2 / 2))

# ------------------------------------------------------------------ Part B
ald, bed, X, SW = sp.symbols("alpha_dot beta_dot X SW")           # X = d(s omega)/dt, SW = s omega
rho_c, p_perp = 2 * SW, SW                                                  # derived cross term (paper convention), p_parallel = 0
cons = 2 * X + ald * rho_c + 2 * bed * (rho_c + p_perp)                     # d_t rho + H_x (rho + p_x) + 2 H_y (rho + p_perp) = 0, rho = 2 s omega
req = sp.simplify(sp.solve(sp.Eq(cons, 0), X)[0] / SW)
print("conservation of the cross term requires d ln(s omega)/dt =", req)
actual = -(ald + 2 * bed) - 2 * bed                                          # s ~ n ~ e^{-alpha-2beta} and omega ~ e^{-2 beta} (separately conserved)
print("with s ~ n and omega ~ a_perp^-2 (separately conserved): d ln(s omega)/dt =", sp.simplify(actual), " ; mismatch =", sp.simplify(req - actual))
out["conservation_required"] = str(req); out["conservation_actual"] = str(sp.simplify(actual)); out["conservation_mismatch"] = str(sp.simplify(req - actual))
# vorticity equation for a Weyssenhoff fluid (paper Eq. 64, geodesic, spatially homogeneous, Frenkel):
#   omega_dot_ab = -(2/3) Theta omega_ab + 2 sigma_c[a omega_b]^c + kappa (sigma + omega)_[a^c s_b]c ,  omega = omega~ - kappa s / 2
# (index placement and sign of the last term as read from the paper's Eq. 64; flagged as a reading, not independently verified)
sigma_p = sp.Symbol("sigma_p")
Wt, Ss = sp.symbols("Wt Ss")
w_tors = Wt - kap * Ss / 2
mat = lambda w_: sp.Matrix([[0, w_, 0], [-w_, 0, 0], [0, 0, 0]])
Om_t, S_t = mat(w_tors), mat(Ss)
SigM = sp.diag(sigma_p, sigma_p, -2 * sigma_p)
term_sigma = sum(SigM[0, c] * Om_t[1, c] for c in range(3)) - sum(SigM[1, c] * Om_t[0, c] for c in range(3))
cross = sp.Rational(1, 2) * kap * (sum(((SigM + Om_t)[0, c]) * S_t[1, c] for c in range(3)) - sum(((SigM + Om_t)[1, c]) * S_t[0, c] for c in range(3)))
Hperp = sigma_p + th / 3
w_dot_torsionful = -sp.Rational(2, 3) * th * w_tors + term_sigma + cross
w_dot_torsionful = sp.simplify(w_dot_torsionful.subs(sigma_p, Hperp - th / 3))
# with s_dot = -Theta s (s ~ n): omega~ = omega + kappa s / 2 -> omega~_dot = omega_dot - kappa s_dot / 2
Wt_dot = sp.simplify(w_dot_torsionful.subs(Wt, Wt) + kap * th * Ss / 2)
print("omega~_dot (reading of Eq. 64, s ~ n) =", sp.collect(sp.expand(Wt_dot), [Wt, Ss]))
out["omega_tilde_dot_from_eq64"] = str(sp.expand(Wt_dot))

# ------------------------------------------------------------------ Part C: what-if bound on the Bianchi-I critical spin strength
import phase3c2_bianchi as P
base_v2 = 0.2
rows = []
for f in (0.0, 0.5, 0.9, 0.99):
    Bc = P.critical_B(v0sq=base_v2 * (1 - f), rhoB0=0.0, shear=0.0)
    rows.append(dict(f_cancelled=f, B_critical=Bc))
    print("what-if: fraction %.2f of the rotation-induced stress cancelled (v0^2 = %.3g effective): B_c = %.4g" % (f, base_v2 * (1 - f), Bc))
out["what_if"] = rows
json.dump(out, open(os.path.join(HERE, "phase3c4_polarized.json"), "w"), indent=1)
