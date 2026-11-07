# -*- coding: utf-8 -*-
"""Phase 3C.3: local (WKB) wave analysis of the magnetized effective spin fluid.  Relativistic ideal MHD, locally flat, background at rest with magnetic field b along z; isentropic EOS rho(n), p(n)
(effective spin fluid: rho = n - B n^2, p = -B n^2; enthalpy density w = rho + p = n (1 - 2 B n); c_s^2 = dp/drho = -2Bn/(1 - 2Bn)).

Equations (c = 1):  d_mu (n u^mu) = 0 ;  d_mu T^{mu nu} = 0 with T^{mu nu} = (w + b^2) u^mu u^nu + (p + b^2/2) eta^{mu nu} - b^mu b^nu ;  d_mu (b^mu u^nu - b^nu u^mu) = 0 .
Plane waves exp(i (k.x - omega t)), k at angle psi to b.  Derivation by sympy; the result is compared with the known magnetosonic dispersion relation
    (omega^2 - k_par^2 v_A^2) [ omega^4 - omega^2 k^2 (v_A^2 + c_s^2 (1 - v_A^2)) + k^2 k_par^2 v_A^2 c_s^2 ] = 0 ,   v_A^2 = b^2 / (w + b^2).
"""
import json
import os

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
w, b, wp, pp, n0, om, k, psi = sp.symbols("w b wp pp n0 omega k psi", real=True)       # wp = dw/dn, pp = dp/dn
dn, ux, uy, uz, bx, by, bz = sp.symbols("dn ux uy uz bx by bz")
kx, kz = k * sp.sin(psi), k * sp.cos(psi)
I = sp.I
# perturbation amplitudes (common factor exp(i(kx x + kz z - om t))): d/dt -> -i om, d/dx -> i kx, d/dz -> i kz, d/dy -> 0
dt_, dx_, dz_ = -I * om, I * kx, I * kz
# background: u^mu = (1,0,0,0), b^mu = (0,0,0,b); perturbed: u^mu = (1, ux, uy, uz) (to first order u^0 = 1), b^0 = b uz (u_mu b^mu = 0), b^i = (bx, by, b + bz)
db0 = b * uz
# (1) number conservation:  d_t dn + n0 (d_i u^i) = 0
eq1 = dt_ * dn + n0 * (dx_ * ux + dz_ * uz)
# stress tensor first-order parts.  W = w + b^2 (delta: wp dn + 2 b bz),  P = p + b^2/2 (delta: pp dn + b bz)
dW = wp * dn + 2 * b * bz
dP = pp * dn + b * bz
Wtot = w + b ** 2
# T^{0 0} = W u^0 u^0 - P (eta^{00} = -1) - b^0 b^0 -> first order: dW - (-dP)... careful: T^{00} = W - P - (b^0)^2 ; delta T^{00} = dW - dP
dT00 = dW - dP
# T^{0 i} = W u^i - b^0 b^i : first order W u^i - b^0_pert * b^i_bg
dT0x = Wtot * ux
dT0y = Wtot * uy
dT0z = Wtot * uz - db0 * b
# T^{i j} = P delta_ij - b^i b^j  (W u^i u^j second order)
dTxx = dP; dTyy = dP; dTzz = dP - 2 * b * bz
dTxz = -bx * b
dTyz = -by * b
dTxy = 0
# momentum / energy equations d_mu T^{mu nu} = 0
eq_E = dt_ * dT00 + dx_ * dT0x + dz_ * dT0z
eq_x = dt_ * dT0x + dx_ * dTxx + dz_ * dTxz
eq_y = dt_ * dT0y + dz_ * dTyz
eq_z = dt_ * dT0z + dx_ * dTxz + dz_ * dTzz
# induction  d_mu (b^mu u^nu - b^nu u^mu) = 0, nu = i : d_t (b^0 u^i - b^i u^0) + d_j (b^j u^i - b^i u^j) = 0 ; first order with b^0 = b uz, u^0 = 1
# nu = x:  d_t(-bx) + d_z (b ux - bx*0 ...)  -> -dt bx + d_z (b ux) - d_x (b^x u^z) [2nd order] ... first order: -dt bx + b d_z ux - d_j(bx_bg...) b^x_bg = 0 -> only d_z
ind_x = -dt_ * bx + b * dz_ * ux
ind_y = -dt_ * by + b * dz_ * uy
# nu = z: d_t(b^0 u^z - b^z u^0) + d_j(b^j u^z - b^z u^j): first order: -dt bz + b d_z uz - b d_j u^j = -dt bz - b (dx_ ux)  [b d_z uz cancels]
ind_z = -dt_ * bz - b * dx_ * ux
unk = [dn, ux, uy, uz, bx, by, bz]
eqs = [eq1, eq_E, eq_x, eq_y, eq_z, ind_x, ind_y, ind_z]
Mx = sp.Matrix([[sp.diff(e, u) for u in unk] for e in eqs])
# The system has 8 equations for 7 unknowns (the equation d_mu b^mu-type constraint is automatically satisfied for the divergence condition); use the 7 independent ones, drop the redundant energy equation after checking rank
sub_n = {pp: sp.Symbol("cs2") * 0}  # placeholder, replaced below
cs2 = sp.Symbol("cs2", real=True)
rows = [eq1, eq_x, eq_y, eq_z, ind_x, ind_y, ind_z]
Mx7 = sp.Matrix([[sp.diff(e, u) for u in unk] for e in rows])
det7 = sp.factor(sp.simplify(Mx7.det()))
print("7x7 determinant factors:", det7)
# relations: pp = dp/dn = cs2 * drho/dn = cs2 * w / n0 ;  wp = dw/dn  (appears only in the energy equation, which is implied by the others for an isentropic fluid)
det = sp.simplify(det7.subs(pp, cs2 * w / n0))
print("with p' = c_s^2 w / n:", sp.factor(det))
vA2 = sp.Symbol("vA2", positive=True)
det_v = sp.simplify(det.subs(b, sp.sqrt(vA2 * w / (1 - vA2))))
print("in terms of v_A^2:", sp.factor(det_v))
res = dict(determinant=str(sp.factor(det_v)))
# product of the fast/slow roots (x = omega^2/k^2): from the quartic factor
x = sp.Symbol("x")
quartic = x ** 2 - x * (vA2 + cs2 * (1 - vA2 * sp.sin(psi) ** 2)) + sp.cos(psi) ** 2 * vA2 * cs2
print("magnetosonic quartic (derived):", quartic)
factor_check = sp.simplify(det_v * (vA2 - 1) ** 2 / (-sp.I * om * w ** 3) / ((k ** 2 * vA2 * sp.cos(psi) ** 2 - om ** 2)) + k ** 4 * (quartic.subs(x, om ** 2 / k ** 2)))
print("factorization check (should be 0):", factor_check)
res["factorization_check"] = str(factor_check)
# numerical scan: for c_s^2 < 0, is omega^2 < 0 for some branch at each angle?
rows_out = []
for cs2v in (-0.005, -0.05, -0.5):
    for vA2v in (0.0, 0.1, 0.5, 0.9, 0.999):
        for psi_v in (0.0, 0.3, 0.8, 1.2, 1.5, np.pi / 2):
            q = sp.Poly(quartic.subs({cs2: cs2v, vA2: vA2v, psi: psi_v}), x)
            r = np.roots([float(c) for c in q.all_coeffs()])
            rows_out.append(dict(cs2=cs2v, vA2=vA2v, psi=float(psi_v), omega2_over_k2=[float(np.real(z)) for z in r], alfven=float(vA2v * np.cos(psi_v) ** 2)))
res["scan"] = rows_out
unstable = [r for r in rows_out if min(r["omega2_over_k2"]) < -1e-12]
res["fraction_unstable"] = len(unstable) / len(rows_out)
perp_stable = [r for r in rows_out if abs(r["psi"] - np.pi / 2) < 1e-9 and min(r["omega2_over_k2"]) >= -1e-12]
res["perp_stable"] = [(r["cs2"], r["vA2"]) for r in perp_stable]
print("scan: %d of %d (c_s^2<0) cases have an unstable branch (product of the fast and slow roots c_s^2 v_A^2 cos^2 psi < 0 for every psi != pi/2); exactly perpendicular propagation is stable only for v_A^2 > -c_s^2/(1-c_s^2): %s" % (len(unstable), len(rows_out), res["perp_stable"]))
# enthalpy regularisation: total inertia w + b^2 does not vanish at n* when b != 0
B_ = 2e-7; n_ = np.linspace(0.5 / B_, 1.5 / B_, 9)
wv = n_ * (1 - 2 * B_ * n_)
res["inertia_at_nstar"] = dict(w_min=float(np.min(np.abs(wv))), with_field_b2_1e5=float(np.min(wv + 1e5)))
json.dump(res, open(os.path.join(HERE, "phase3c3_mhd.json"), "w"), indent=1)
