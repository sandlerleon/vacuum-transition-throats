# -*- coding: utf-8 -*-
"""Phase 3D.0: Einstein-Cartan-Dirac (ECD) feasibility.  Units hbar = c = 1, kappa = 8 pi G, Planck units G = 1 where numbers are quoted.

Part A  Fundamental equations.  (i) Dirac sector (explicit matrices, signature (+,-,-,-)): the contortion couples ONLY through its totally antisymmetric part, with coefficient |c_K| = 1/4:
            L_K = (i/4) K_{abc} psibar gamma^{abc} psi = c_K K_{abc} eps^{abcd} A_d ,  A^mu = psibar gamma^mu gamma5 psi.
        (ii) Gravity sector (symbolic): for K_{abc} = eps_{abcd} a^d the Riemann-Cartan scalar is R = R~ + 6 a.a  (signature (-,+,+,+)).
        (iii) Eliminating the algebraic torsion:  L_eff = -(3 kappa/16) A_a A^a  [(-,+,+,+)]  =  +lam [ (A^0)^2 - |A|^2 ] ,  lam = 3 kappa/16 = 3 pi G /2.
Part B  Matter states.  Hartree-Fock energy density of  H_int = -lam sum_mu eta_mumu (A^mu)^2  in a degenerate gas with local density matrix D:
            <:(psi^dag M psi)(psi^dag M psi):> = Tr(M D)^2 - Tr(M D M D) ,   eps_int = -lam sum_mu eta_mumu [Tr(M_mu D)^2 - Tr(M_mu D M_mu D)],  M_mu = gamma^0 gamma^mu gamma5.
        Result: every state gives eps_int = -B n^2, p_int = -B n^2 (the Weyssenhoff form), with B/lam = 1/2 (unpolarized massless), 2 (single-helicity Weyl), 3/2 (unpolarized nonrelativistic), 0 (fully polarized nonrelativistic: Pauli).
Part C  Barotropic baseline rho(n) = C n^(4/3) - B n^2 (flat FLRW), exact bounce, characteristic densities, curvature at the bounce, and the no-go theorem for ANY barotropic fluid bounce.
Part D  Embedding facts (Misner-Sharp mass of the homogeneous ball).
"""
import itertools
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}

# ------------------------------------------------------------------------------------------------ Part A (i): Dirac sector
I2 = np.eye(2); Z2 = np.zeros((2, 2))
sig = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
g0 = np.block([[I2, Z2], [Z2, -I2]]).astype(complex)
gs = [np.block([[Z2, s], [-s, Z2]]) for s in sig]
gam = [g0] + gs                                                     # gamma^mu (upper index)
g5 = 1j * g0 @ gs[0] @ gs[1] @ gs[2]
eta_m = np.diag([1.0, -1, -1, -1])                                  # (+,-,-,-)
assert all(np.allclose(gam[a] @ gam[b] + gam[b] @ gam[a], 2 * eta_m[a, b] * np.eye(4)) for a in range(4) for b in range(4))
def asym_prod(idx):
    M = np.zeros((4, 4), complex)
    for perm in itertools.permutations(range(len(idx))):
        sgn = np.linalg.det(np.eye(len(idx))[list(perm)])
        P = np.eye(4, dtype=complex)
        for j in perm:
            P = P @ gam[idx[j]]
        M += sgn * P
    return M / np.math.factorial(len(idx)) if hasattr(np, "math") else M / float(np.prod(range(1, len(idx) + 1)))
levi = lambda a, b, c, d: int(sp.LeviCivita(a, b, c, d))
# gamma^{abc} = c0 eps^{abcd} gamma_d gamma5  (eps^{0123} = +1)
ratios = []
for a, b, c in itertools.permutations(range(4), 3):
    d = [x for x in range(4) if x not in (a, b, c)][0]
    lhs = asym_prod((a, b, c))
    rhs = levi(a, b, c, d) * eta_m[d, d] * gam[d] @ g5
    k = np.trace(lhs @ rhs.conj().T) / np.trace(rhs @ rhs.conj().T)
    assert np.allclose(lhs, k * rhs)
    ratios.append(k)
c0 = ratios[0]
assert np.allclose(ratios, c0)
print("A: gamma^{abc} = %s * eps^{abcd} gamma_d gamma5 for all 24 index triples" % np.round(c0, 12))
out["A_gamma3_coefficient"] = [float(c0.real), float(c0.imag)]
# contortion coupling operator:  L_K = psi^dag O psi ,  O = (i/8) K_{abc} gamma^0 { gamma^c , gamma^{ab} },  gamma^{ab} = (1/2)[gamma^a, gamma^b]
def gab(a, b): return 0.5 * (gam[a] @ gam[b] - gam[b] @ gam[a])
def op_K(K):
    O = np.zeros((4, 4), complex)
    for a, b, c in itertools.product(range(4), repeat=3):
        if K[a, b, c] != 0:
            O += 1j / 8 * K[a, b, c] * g0 @ (gam[c] @ gab(a, b) + gab(a, b) @ gam[c])
    return O
rng = np.random.default_rng(11)
K = rng.normal(size=(4, 4, 4)); K = K - K.transpose(1, 0, 2)               # antisymmetric in (a,b): metric-compatible contortion
Kta = np.zeros_like(K)
for a, b, c in itertools.product(range(4), repeat=3):
    Kta[a, b, c] = sum(np.linalg.det(np.eye(3)[list(p)]) * K[[a, b, c][p[0]], [a, b, c][p[1]], [a, b, c][p[2]]] for p in itertools.permutations(range(3))) / 6
dev = np.max(np.abs(op_K(K) - op_K(Kta)))
print("A: coupling of a generic contortion minus coupling of its totally antisymmetric part:", dev, " (vector and tensor torsion decouple)")
out["A_nonaxial_decoupling_residual"] = float(dev)
# coefficient: O(K_[abc]) = c_K K_{abc} eps^{abcd} gamma^0 gamma_d gamma5 ,  A_d = psibar gamma_d gamma5 psi = psi^dag gamma^0 gamma_d gamma5 psi
S_eps = np.zeros((4, 4), complex)
for a, b, c, d in itertools.product(range(4), repeat=4):
    if levi(a, b, c, d) != 0:
        S_eps += Kta[a, b, c] * levi(a, b, c, d) * g0 @ (eta_m[d, d] * gam[d]) @ g5
cK = np.trace(op_K(Kta) @ S_eps.conj().T) / np.trace(S_eps @ S_eps.conj().T)
assert np.allclose(op_K(Kta), cK * S_eps)
print("A: |c_K| =", abs(cK), " (expected 1/4)")
out["A_cK_abs"] = float(abs(cK))
# ------------------------------------------------------------------------------------------------ Part A (ii)/(iii): gravity sector and elimination of torsion
etap = sp.diag(-1, 1, 1, 1)                                         # (-,+,+,+)
a_up = sp.symbols("a0:4", real=True)
A_dn = sp.symbols("A0:4", real=True)
a_dn = [sum(etap[i, j] * a_up[j] for j in range(4)) for i in range(4)]
eps_dn = lambda i, j, k, l: sp.LeviCivita(i, j, k, l)
def eps_up(i, j, k, l): return sp.LeviCivita(i, j, k, l) * etap[i, i] * etap[j, j] * etap[k, k] * etap[l, l]
Kd = lambda i, j, k: sum(eps_dn(i, j, k, l) * a_up[l] for l in range(4))          # K_{abc} = eps_{abcd} a^d
Ku = lambda i, j, k: Kd(i, j, k) * etap[i, i] * etap[j, j] * etap[k, k]               # all-up (diagonal metric)
# quadratic part of the Riemann-Cartan scalar:  K^a_{ca} K^{cb}_b - K^a_{cb} K^{cb}_a   (third index = form index)
Kmix = lambda i, j, k: Kd(i, j, k) * etap[i, i]                                    # K^a_{bc}
Kmix2 = lambda i, j, k: Kd(i, j, k) * etap[i, i] * etap[j, j]                      # K^{ab}_c
Q = sum(Kmix(a, c, a) * Kmix2(c, b, b) for a in range(4) for b in range(4) for c in range(4)) \
    - sum(Kmix(a, c, b) * Kmix2(c, b, a) for a in range(4) for b in range(4) for c in range(4))
aa = sum(a_up[i] * a_dn[i] for i in range(4))
Q_over_aa = sp.simplify(sp.expand(Q) / sp.expand(aa))
print("A: quadratic part of R(Gamma) / (a.a) =", Q_over_aa)
out["A_R_quadratic_over_aa"] = str(Q_over_aa)
kap, bet = sp.symbols("kappa beta", positive=True)
pair = sum(Kd(i, j, k) * eps_up(i, j, k, l) * A_dn[l] for i in range(4) for j in range(4) for k in range(4) for l in range(4))
aA = sum(a_up[i] * A_dn[i] for i in range(4))
beta_coef = sp.simplify(sp.expand(pair) / sp.expand(aA))
print("A: K_{abc} eps^{abcd} A_d / (a.A) =", beta_coef)
out["A_pairing_over_aA"] = str(beta_coef)
# L(a) = (1/(2 kappa)) Q + cK * pair ;  with |cK| = 1/4 ; eliminate a (solve the gradient with respect to a^0 .. a^3)
cKs = sp.Rational(1, 4)
La = sp.expand(Q / (2 * kap) + cKs * pair)
sol = sp.solve([sp.diff(La, a_up[i]) for i in range(4)], a_up, dict=True)[0]
Leff = sp.simplify(sp.expand(La.subs(sol)))
AA = sum(A_dn[i] * A_dn[i] * etap[i, i] for i in range(4))
coef = sp.simplify(Leff / AA)
print("A: torsion eliminated:  L_eff = (%s) A_a A^a   [(-,+,+,+)] ,  expected -3 kappa/16" % coef)
out["A_Leff_coefficient_times_AaAa"] = str(coef)
lam_sym = sp.Rational(3, 16) * kap
assert sp.simplify(coef + lam_sym) == 0

# ------------------------------------------------------------------------------------------------ Part B: Hartree-Fock table
def eps_int_over_lam_n2(D, n):
    tot = 0.0
    for mu in range(4):
        M = g0 @ gam[mu] @ g5
        tot += eta_m[mu, mu] * (np.trace(M @ D) ** 2 - np.trace(M @ D @ M @ D))
    return float((-tot / n ** 2).real)
n = 1.0
PR = (np.eye(4) + g5) / 2
states = {
    "unpolarized massless (g=2, particles)": n / 4 * np.eye(4),
    "single-helicity Weyl (g=1)": n / 2 * PR,
    "unpolarized nonrelativistic": n / 2 * np.diag([1, 1, 0, 0]).astype(complex),
    "polarized NR, fraction 0.5": n * np.diag([0.75, 0.25, 0, 0]).astype(complex),
    "fully polarized NR": n * np.diag([1, 0, 0, 0]).astype(complex),
}
hf = {}
for k_, D in states.items():
    c_ = eps_int_over_lam_n2(D, n)
    hf[k_] = c_
    print("B: %-40s eps_int / (lam n^2) = %+.4f   (B/lam = %.4f)" % (k_, c_, -c_))
out["B_hf_eps_over_lam_n2"] = hf
# Weyssenhoff comparison: B_W = kappa sigma^2/4 with sigma = 1/2  ->  kappa/16 = lam/3
out["B_weyssenhoff_over_lam"] = 1.0 / 3.0

# ------------------------------------------------------------------------------------------------ Part C: homogeneous baseline (Planck units, kappa = 8 pi)
kappa = 8 * np.pi
lam = 3 * kappa / 16
res = {}
for label, g, Bfac in (("unpolarized massless (g=2)", 2, 0.5), ("single-helicity Weyl (g=1)", 1, 2.0)):
    Bc = Bfac * lam
    a_ = 6 * np.pi ** 2 / g
    C = 0.75 * a_ ** (1 / 3)
    rho = lambda n_, C=C, Bc=Bc: C * n_ ** (4 / 3) - Bc * n_ ** 2
    p_ = lambda n_, C=C, Bc=Bc: C / 3 * n_ ** (4 / 3) - Bc * n_ ** 2
    nb = (C / Bc) ** 1.5; nstar = (2 * C / (3 * Bc)) ** 1.5; nsp = (2 * C / (9 * Bc)) ** 1.5
    kFb = (a_ * nb) ** (1 / 3)
    # integrate the contracting flat FLRW universe through the bounce: variables (ln a, H), n = n0 exp(-3 ln a), Hdot = -(kappa/2)(rho + p)
    n0 = 0.02 * nb
    H0 = -np.sqrt(kappa * rho(n0) / 3)
    def f(t, y):
        la, H = y; n_ = n0 * np.exp(-3 * la)
        return [H, -0.5 * kappa * (rho(n_) + p_(n_))]
    ev = lambda t, y: y[1]; ev.direction = 1
    sol = solve_ivp(f, (0, 1e6), [0.0, H0], method="DOP853", rtol=1e-12, atol=1e-14, events=ev, dense_output=True, max_step=50.0)
    tb = sol.t_events[0][0]; lb = sol.y_events[0][0][0]
    nbb = n0 * np.exp(-3 * lb)
    ts = np.linspace(0, tb, 4000); ys = sol.sol(ts); nn = n0 * np.exp(-3 * ys[0])
    cons = np.max(np.abs(ys[1] ** 2 - kappa * rho(nn) / 3)) / (kappa * rho(nn).max() / 3)
    # fluid sign structure along the contraction
    mu = (4 / 3) * C * nn ** (1 / 3) - 2 * Bc * nn
    dpdn = (4 / 9) * C * nn ** (1 / 3) - 2 * Bc * nn
    t_star = ts[np.argmax(mu < 0)] if np.any(mu < 0) else None
    t_sp = ts[np.argmax(dpdn < 0)] if np.any(dpdn < 0) else None
    pb = p_(nbb)
    Rb = kappa * (rho(nbb) - 3 * pb)
    res[label] = dict(B_over_lam=Bfac, C=C, B=Bc, n_bounce_analytic=nb, n_bounce_numeric=float(nbb), n_star_mu_zero=nstar, n_spinodal=nsp,
                      nstar_over_nb=nstar / nb, nsp_over_nb=nsp / nb, kF_over_Planck_at_bounce=float(kFb), eps_int_over_eps_kin_at_spinodal=float(Bc * nsp ** (2 / 3) / C),
                      p_at_bounce=float(pb), ricci_scalar_at_bounce=float(Rb), constraint_residual=float(cons), t_to_bounce=float(tb), t_mu_zero=None if t_star is None else float(t_star), t_spinodal=None if t_sp is None else float(t_sp))
    print("C: %s: n_b = %.4g (numeric %.4g), n*(mu=0) = %.3f n_b, n_spinodal = %.3f n_b (eps_int/eps_kin there = %.3f), kF(bounce) = %.3g, R~(bounce) = %.4g, p(bounce) = %.4g, constraint residual %.1e" %
          (label, nb, nbb, nstar / nb, nsp / nb, Bc * nsp ** (2 / 3) / C, kFb, Rb, pb, cons))
out["C_baseline"] = res

# ---- no-go theorem for barotropic fluids: F(n) = rho(n) - c n^(2/3) (c = 0 flat, c > 0 closed), H^2 = kappa F/3.  Bounce at F(n_b) = 0 with F'(n_b) < 0 (a-ddot > 0 iff F' < 0).
#      If F' > 0 at some lower density (positive chemical potential in the dilute regime), the mean-value theorem gives n* with F'(n*) = 0 (rho + p = n rho' = 0 for c = 0) and F'' < 0 somewhere below n_b,
#      i.e. dp/dn = n rho'' < -(2c/9) n^(-1/3) <= 0 : negative compressibility.  Random-EOS verification:
rng = np.random.default_rng(5)
nn_grid = np.linspace(1e-4, 1.0, 20001)
cnt = dict(tested=0, bounces=0, mu_cross=0, neg_compress=0)
for trial in range(4000):
    kk = rng.uniform(1.0, 3.0, 3); cc = rng.normal(size=3)
    cc[0] = abs(cc[0]) + 0.1
    curv = rng.choice([0.0, rng.uniform(0.0, 0.5)])
    F = cc[0] * nn_grid ** kk[0] + cc[1] * nn_grid ** kk[1] + cc[2] * nn_grid ** kk[2] - curv * nn_grid ** (2 / 3)
    dF = np.gradient(F, nn_grid)
    cnt["tested"] += 1
    # positive-F region ending at a root with dF < 0
    idx = np.where((F[:-1] > 0) & (F[1:] <= 0))[0]
    if len(idx) == 0:
        continue
    j = idx[0]
    if dF[j] >= 0 or not np.any(dF[:j] > 0):
        continue
    cnt["bounces"] += 1
    if np.any(np.diff(np.sign(dF[:j])) < 0):
        cnt["mu_cross"] += 1
    d2F = np.gradient(dF, nn_grid)
    if np.any(d2F[:j] < 0):
        cnt["neg_compress"] += 1
print("C: random barotropic EOS (positive mu at low density, bounce with F'<0): ", cnt)
out["C_random_theorem_check"] = cnt

# ------------------------------------------------------------------------------------------------ Part D: embedding facts
# homogeneous ball glued to vacuum: Misner-Sharp mass m = (kappa/6) rho R^3 (flat), = (a sin^3 chi0 / 2)(1 + adot^2) (closed, boundary at chi0): m ~ rho/n (closed) ; trapped iff 2m/R = sin^2 chi0 (1 + adot^2) > 1.
chi = np.linspace(0.05, np.pi, 200)
out["D_bounce_not_trapped"] = bool(np.all(np.sin(chi) ** 2 * (1 + 0.0) <= 1.0 + 1e-15))
out["D_flat_MS_mass_at_bounce"] = 0.0
print("D: at a bounce (adot = 0) the closed ball has 2m/R = sin^2(chi0) <= 1 (never trapped); flat ball: m_MS = (kappa/6) rho R^3 = 0 at rho = 0.")
json.dump(out, open(os.path.join(HERE, "phase3d0_ecd.json"), "w"), indent=1, default=str)
