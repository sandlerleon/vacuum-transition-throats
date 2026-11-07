# -*- coding: utf-8 -*-
"""Phase 3C.5: covariant spin transport and constraint closure for the anti-aligned polarized-spin configuration s = -2 omega~/kappa of Phase 3C.4.
Conventions as in Pasmatsiou, Tsagas & Barrow (arXiv:1611.07878): signature (-,+,+,+), kappa = 8 pi G, s^2 = s_ab s^ab/2, Frenkel condition s_ab u^b = 0, Riemannian (tilde) vorticity omega~.

Part A  spin-tensor algebra: aligned polarized spin satisfies the Frenkel condition; s is NOT a free parameter (s = n sigma0).
Part B  spin transport: the angular-momentum balance + Frenkel condition forces Fermi-Walker transport (the coefficient is fixed by preserving the Frenkel condition);
        consequence d ln s / dtau = -Theta EXACTLY, for any acceleration.
Part C  kinematic identities checked on an explicit metric (flat space, rigidly rotating observers): [D_a, D_b] f = c * 2 omega_ab fdot, and the vector (momentum) constraint sign structure.
Part D  closure theorem in the homogeneous axisymmetric class: s ~ n, Helmholtz omega~ ~ a_perp^-2 and the stress-energy conservation of the cross term (3C.4) are three equations
        for two rates (alpha_dot, beta_dot): they force beta_dot = 0.  The derivative of s + 2 omega~/kappa on the cancellation surface is  -alpha_dot * s  (and s (hdot/h - alpha_dot) with enthalpy).
Part E  numerical illustration along an Einstein-Cartan Bianchi-I collapse; cancellation of the source term is not cancellation of the shear.
"""
import json
import os

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}

# ---------------------------------------------------------------------------------------------- Part A
s, th = sp.symbols("s Theta", real=True)
Ax, Ay, Az = sp.symbols("A_x A_y A_z", real=True)
eta = sp.diag(-1, 1, 1, 1)
u_up = sp.Matrix([1, 0, 0, 0])
u_dn = eta * u_up
# aligned spin: s^{ab} = s eps^{abc} e_c, e = z : only s^{12} = s, s^{21} = -s (spatial, orthogonal to u)
S_up = sp.zeros(4, 4)
S_up[1, 2], S_up[2, 1] = s, -s
frenkel = sp.simplify(S_up * u_dn)
S_dn = eta * S_up * eta
s2 = sp.simplify(sum(S_dn[a, b] * S_up[a, b] for a in range(4) for b in range(4)) / 2)
out["A_frenkel_residual"] = [str(x) for x in frenkel]
out["A_s_squared"] = str(s2)
print("A: Frenkel residual s^{ab} u_b =", list(frenkel), "; s_ab s^ab/2 =", s2)

# ---------------------------------------------------------------------------------------------- Part B
lam = sp.Symbol("lambda")
A_dn = sp.Matrix([0, Ax, Ay, Az])                                          # acceleration, A_a u^a = 0
# general spin-transport ansatz allowed by antisymmetry:  Sdot^{ab} = -Theta S^{ab} + lam (u^b S^{ac} A_c - u^a S^{bc} A_c)
Sdot = sp.zeros(4, 4)
for a in range(4):
    for b in range(4):
        Sdot[a, b] = -th * S_up[a, b] + lam * (u_up[b] * sum(S_up[a, c] * A_dn[c] for c in range(4)) - u_up[a] * sum(S_up[b, c] * A_dn[c] for c in range(4)))
# d/dtau (S^{ab} u_b) = Sdot^{ab} u_b + S^{ab} A_b   (u_b dot = A_b)
frenkel_rate = sp.simplify(Sdot * u_dn + S_up * A_dn)
sol_lam = sp.solve(list(frenkel_rate), lam, dict=True)
print("B: d/dtau(S^{ab}u_b) =", list(frenkel_rate), " -> lambda =", sol_lam)
out["B_frenkel_rate"] = [str(x) for x in frenkel_rate]
out["B_lambda_forced"] = [str(x) for x in sol_lam]
Sdot1 = Sdot.subs(lam, 1)
d_s2 = sp.simplify(sum((eta * S_up * eta)[a, b] * Sdot1[a, b] for a in range(4) for b in range(4)))   # d(S_ab S^ab/2)/dtau = S_ab Sdot^ab
d_ln_s = sp.simplify(d_s2 / (2 * s2))
print("B: d ln s / dtau =", d_ln_s, " (independent of the acceleration)")
out["B_dlns_over_dtau"] = str(d_ln_s)

# ---------------------------------------------------------------------------------------------- Part C
t, x, y, z = sp.symbols("t x y z", real=True)
W = sp.Rational(1, 2)
X = [t, x, y, z]
gam = 1 / sp.sqrt(1 - W ** 2 * (x ** 2 + y ** 2))
uU = sp.Matrix([gam, -gam * W * y, gam * W * x, 0])
uD = eta * uU
hUD = sp.Matrix(4, 4, lambda a, b: (1 if a == b else 0) + uU[a] * uD[b])      # h^a_b
du = sp.Matrix(4, 4, lambda b, a: sp.diff(uD[a], X[b]))                       # du[b,a] = d_b u_a
A_dn_flat = sp.Matrix([sum(uU[b] * du[b, a] for b in range(4)) for a in range(4)])
def hproj(T):                                                                 # project both indices of a (0,2) tensor T[a,b]
    return sp.Matrix(4, 4, lambda a, b: sum(hUD[c, a] * hUD[d, b] * T[c, d] for c in range(4) for d in range(4)))
omega_dn = hproj(sp.Matrix(4, 4, lambda a, b: (du[b, a] - du[a, b]) / 2))     # omega_ab = h h  d_[b u_a]
# generic polynomial scalar
rng = np.random.default_rng(3)
mons = [1, t, x, y, z, t * x, t * y, x * y, x * z, y * z, t ** 2, x ** 2, y ** 2, z ** 2, t * x * y, x * y * z, x ** 3, y ** 3]
f = sum(sp.Rational(int(rng.integers(-5, 6)), int(rng.integers(1, 4))) * m for m in mons)
df = sp.Matrix([sp.diff(f, X[a]) for a in range(4)])
Df = sp.Matrix([sum(hUD[b, a] * df[b] for b in range(4)) for a in range(4)])          # D_a f = h_a^b d_b f  (h^b_a symmetric here)
DDf = hproj(sp.Matrix(4, 4, lambda c, d: sp.diff(Df[d], X[c])))                       # D_a D_b f
fdot = sum(uU[a] * df[a] for a in range(4))
pt = {t: sp.Rational(1, 3), x: sp.Rational(1, 3), y: sp.Rational(1, 5), z: sp.Rational(2, 7)}
anti = ((DDf - DDf.T) / 2).subs(pt).applyfunc(sp.nsimplify)
omf = (omega_dn * fdot).subs(pt).applyfunc(sp.nsimplify)
ratio = [sp.N(anti[a, b] / omf[a, b]) for a in range(4) for b in range(a + 1, 4) if abs(omf[a, b]) > 1e-12]
print("C: D_[a D_b] f / (omega_ab fdot) over nonzero components:", ratio)
out["C_commutator_ratio"] = [float(r) for r in ratio]
# vector constraint: curl omega + k eps A omega, flat space (Theta = sigma = 0 for rigid rotation, Ricci = 0)
epsU = lambda a, b, c, d: sp.LeviCivita(a, b, c, d)                            # eta_{abcd} with eta_{0123} = 1 (sqrt(-g) = 1)
eps3_dn = lambda a, b, c: sum(epsU(a, b, c, d) * uU[d] for d in range(4))        # eps_abc = eta_abcd u^d
omU = (eta * omega_dn * eta)                                                    # omega^{ab}
omvec_dn = sp.Matrix([sum(eps3_dn(a, b, c) * omU[b, c] for b in range(4) for c in range(4)) / 2 for a in range(4)])      # omega_a = (1/2) eps_abc omega^bc
omvec_up = eta * omvec_dn
dom = sp.Matrix(4, 4, lambda c, a: sp.diff(omvec_dn[a], X[c]))                    # d_c omega_a
Dom = sp.Matrix(4, 4, lambda b, c: sum(hUD[d, b] * hUD[e, c] * dom[d, e] for d in range(4) for e in range(4)))      # D_b omega_c (rest-space projection)
Domup = eta * Dom * eta
curl_om = sp.Matrix([sum(eps3_dn(a, b, c) * Domup[b, c] for b in range(4) for c in range(4)) for a in range(4)])
Aup = eta * A_dn_flat
AxO = sp.Matrix([sum(eps3_dn(a, b, c) * Aup[b] * omvec_up[c] for b in range(4) for c in range(4)) for a in range(4)])
curl_p = curl_om.subs(pt).applyfunc(sp.nsimplify); AxO_p = AxO.subs(pt).applyfunc(sp.nsimplify)
ks = []
for a in range(4):
    if abs(sp.N(AxO_p[a])) > 1e-12:
        ks.append(float(sp.N(-curl_p[a] / AxO_p[a])))
print("C: curl omega + k eps_abc A^b omega^c = 0  with k =", ks, "(GR: k = 2 for the vector constraint at Theta = sigma = Ricci = 0)")
out["C_vector_constraint_k"] = ks

# ---------------------------------------------------------------------------------------------- Part D
ad, bd, kap = sp.symbols("alpha_dot beta_dot kappa", real=True)
Theta = ad + 2 * bd
# (R2) spin density:  d ln s/dt = -Theta ;  (R3) Helmholtz circulation law:  d ln omega~/dt = -2 beta_dot (a_perp^-2);  (E) cross-term energy balance (3C.4): d ln(s omega~)/dt = -(alpha_dot + 3 beta_dot)
dlns = -Theta
dlnw = -2 * bd
E_req = -(ad + 3 * bd)
mismatch = sp.simplify(E_req - (dlns + dlnw))
print("D: energy balance required - (spin law + Helmholtz) =", mismatch, " -> closure only if", sp.solve(sp.Eq(mismatch, 0), bd))
out["D_mismatch"] = str(mismatch)
# derivative of (s + 2 omega~/kappa) on the cancellation surface omega~ = -kappa s / 2
ss = sp.Symbol("s")
hh, hd = sp.symbols("h hdot")
omt = -kap * ss / 2
D0 = sp.simplify(dlns * ss + (2 / kap) * (dlnw * omt))
Dh = sp.simplify(dlns * ss + (2 / kap) * ((dlnw - hd / hh) * omt))
print("D: d/dtau (s + 2 omega~/kappa) on s = -2 omega~/kappa :", D0, "   with Kelvin enthalpy factor (h omega~ ~ a_perp^-2):", Dh)
out["D_preservation_derivative"] = str(D0); out["D_preservation_derivative_enthalpy"] = str(Dh)
# trajectory that would keep the surface invariant with effective-fluid enthalpy h = 1 - 2 B n :  alpha_dot / beta_dot
B, n = sp.symbols("B n", positive=True)
hdot = 2 * B * n * Theta                                                        # d(1-2Bn)/dt with ndot = -Theta n
locus = sp.solve(sp.Eq(Dh.subs({hd: hdot, hh: 1 - 2 * B * n}), 0), ad)[0]
print("D: invariant-surface locus (effective-fluid enthalpy): alpha_dot / beta_dot =", sp.simplify(locus / bd))
out["D_locus_ratio"] = str(sp.simplify(locus / bd))
# shear fraction at alpha_dot = 0
sig2 = (ad - bd) ** 2 / 3
Sigma2 = sp.simplify(sig2 / (Theta ** 2 / 3))
print("D: shear fraction sigma^2/(Theta^2/3) =", Sigma2, "; at alpha_dot = 0:", Sigma2.subs(ad, 0))
out["D_shear_fraction_at_alpha_dot_zero"] = str(Sigma2.subs(ad, 0))
# the three laws also imply (without the mismatch-fixing beta_dot = 0):  ratio r = s / omega~ evolves as  d ln r/dt = -alpha_dot
out["D_ratio_law"] = "d ln(s/omega~)/dt = -alpha_dot   (s ~ n, omega~ ~ a_perp^-2)"

# ---------------------------------------------------------------------------------------------- Part E
import phase3c2_bianchi as P
rows = []
for B_ in (0.01, 0.03, 0.1):
    m = P.Model(B=B_)
    r = P.run(m, shear=0.05)
    if r["outcome"] != "bounce":
        rows.append(dict(B=B_, outcome=r["outcome"]))
        continue
    sl = r["sol"]
    al = sl.y[0]; ad_ = sl.y[2]; be_ = sl.y[1]; bd_ = sl.y[3]
    ratio_r = np.exp(-(al - al[0]))                                            # s/omega~ relative to its initial value (cancellation surface at t = 0)
    sig2n = (ad_ - bd_) ** 2 / 3
    V = np.exp(al + 2 * be_)
    rows.append(dict(B=B_, outcome="bounce", V_bounce=r["V_min"], r_over_r0_at_bounce=float(ratio_r[-1]), max_abs_deviation=float(np.max(np.abs(ratio_r - 1.0))), alpha_total_change=float(al[-1] - al[0]),
                     sigma2V2_ratio=float(sig2n[-1] * V[-1] ** 2 / (sig2n[0] * V[0] ** 2))))
    print("E: B = %.3g: V_b = %.4g ; s/omega~ at bounce / initial = %.4g ; max deviation %.3g ; sigma^2 V^2 end/start = %.3g" % (B_, r["V_min"], ratio_r[-1], np.max(np.abs(ratio_r - 1.0)), rows[-1]["sigma2V2_ratio"]))
out["E_bianchi_ratio"] = rows
# GR dust (no source of anisotropic stress): the shear is not removed; sigma^2 V^2 stays constant
m0 = P.Model(B=0.0)
r0 = P.run(m0, shear=0.05)
sl = r0["sol"]
V0 = np.exp(sl.y[0] + 2 * sl.y[1]); s0 = (sl.y[2] - sl.y[3]) ** 2 / 3
k_ = s0 * V0 ** 2
out["E_dust_sigma2V2_spread"] = float((k_.max() - k_.min()) / k_.mean())
out["E_dust_outcome"] = r0["outcome"]
print("E: GR dust, no anisotropic source: sigma^2 V^2 relative spread over the whole collapse =", out["E_dust_sigma2V2_spread"], "(outcome %s): the shear is not removed when the source term vanishes" % r0["outcome"])
json.dump(out, open(os.path.join(HERE, "phase3c5_closure.json"), "w"), indent=1)
