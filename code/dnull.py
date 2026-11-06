# -*- coding: utf-8 -*-
"""Spherically symmetric Einstein-scalar evolution in double-null gauge with a regular centre.

    ds^2 = -f(u,v) du dv + r(u,v)^2 dOmega^2 ,   f = exp(s) ,   G_ab = T_ab ,   8 pi G = 1 ,
    T_ab = d_a p d_b p - g_ab [ (1/2)(dp)^2 + U(p) ] ,        (equations derived symbolically in derive_equations.py)

Evolution equations (r_uv, p_uv, s_uv) and constraint (uu):
    r_uv = - r_u r_v / r - f / (4 r) + U f r / 4
    p_uv = - (p_u r_v + p_v r_u) / r - f U'(p) / 4
    s_uv = [ U f - 2 p_u p_v - 4 r_uv / r ] / 2
    C_uu = r_uu - s_u r_u + (r/2) p_u^2 = 0 ;      C_vv likewise.
Misner-Sharp mass  m = r (1 + 4 r_u r_v / f)  (includes the vacuum energy: m = rho r^3 / 3 for a vacuum of density rho).
u = const are outgoing null rays (tangent d_v), v = const ingoing (tangent d_u). Normal region: r_v > 0, r_u < 0; trapped: r_v < 0, r_u < 0;
anti-trapped (expanding): r_v > 0, r_u > 0.

Scheme (Garfinkle 1995; Gundlach 1997): second-order diamond predictor-corrector on a uniform lattice u_j = j h, v_i = i h (row j valid for i >= j).
Data on the outgoing cone u = 0 from the regular centre, in the gauge r(0, v) = v/2: p(0, v) free, s(0, v) from the vv constraint, s_v = r p_v^2.
The centre point of every new row follows from regularity: r_v = sqrt(f)/2, p_t = 2 p_v (trapezoid in t), iterated.
"""
import math

import numpy as np
from numba import njit

KCON = 6                                   # points next to the centre where sigma is integrated from the constraint


# ------------------------------------------------------------------------------------------------ potential
def make_potential(q=1.0, kappa=24.0, rho_in=3.0):
    """U(psi) = (kappa/4)(psi^2-1)^2 + s psi + d with minima energies U_false = rho_in (psi > 0) and U_true = -q rho_in (psi < 0); psi = p / A."""
    from scipy.optimize import fsolve

    def minima(s):
        # psi near +1 and -1: roots of U'(psi) = kappa psi (psi^2 - 1) + s
        rts = np.roots([kappa, 0.0, -kappa, s])
        rts = np.sort(rts.real)
        return rts[0], rts[-1]

    def eqs(x):
        s, d = x
        pt, pf = minima(s)
        U = lambda p: kappa / 4 * (p * p - 1) ** 2 + s * p + d
        return [U(pf) - rho_in, U(pt) + q * rho_in]

    s, d = fsolve(eqs, [(1 + q) * rho_in / 2, (1 - q) * rho_in / 2])
    pt, pf = minima(s)
    return dict(kappa=kappa, s=float(s), d=float(d), psi_t=float(pt), psi_f=float(pf), rho_in=rho_in, rho_out=-q * rho_in)


@njit(cache=True)
def Upot(p, A, kappa, s, d):
    x = p / A
    return kappa / 4.0 * (x * x - 1.0) ** 2 + s * x + d


@njit(cache=True)
def dUpot(p, A, kappa, s):
    x = p / A
    return (kappa * x * (x * x - 1.0) + s) / A


@njit(cache=True)
def cell(rW, rA, rB, pW, pA, pB, sW, sA, sB, sC, h, A, kappa, s_, d_):
    """Diamond update of r and p: values at C from W (lower), A (same row, previous), B (previous row, next); sigma at C is given (constraint-integrated)."""
    rC = rA + rB - rW
    pC = pA + pB - pW
    for _ in range(3):
        rm = 0.25 * (rW + rA + rB + rC)
        pm = 0.25 * (pW + pA + pB + pC)
        sm = 0.25 * (sW + sA + sB + sC)
        f = math.exp(sm)
        ru = ((rA - rW) + (rC - rB)) / (2.0 * h)
        rv = ((rB - rW) + (rC - rA)) / (2.0 * h)
        pu = ((pA - pW) + (pC - pB)) / (2.0 * h)
        pv = ((pB - pW) + (pC - pA)) / (2.0 * h)
        U = Upot(pm, A, kappa, s_, d_)
        rr = -ru * rv / rm - f / (4.0 * rm) + U * f * rm / 4.0
        pp = -(pu * rv + pv * ru) / rm - f * dUpot(pm, A, kappa, s_) / 4.0
        rC = rA + rB - rW + h * h * rr
        pC = pA + pB - pW + h * h * pp
    return rC, pC


@njit(cache=True)
def sigma_row(r, p, sg, rn, pn, sn_guess, x, j, N, h, A, kappa, s_, d_, K):
    """sigma on row j+1 given the centre value x. Within K points of the centre it is integrated from the vv constraint (r_v > 0 there;
    s_v = (r_vv + r p_v^2 / 2) / r_v, s_v = 0 at the centre by regularity); beyond, the evolution equation
    s_uv = [U f - 2 p_u p_v - 4 r_uv / r] / 2 is used (the constraint divides by r_v, which vanishes on an apparent horizon).
    Also returns G = ln(4 r_v^2) at the centre (elementary flatness), which must equal x."""
    sn = np.zeros(N + 1)
    rv0 = (4.0 * rn[j + 2] - rn[j + 3]) / (2.0 * h)
    G = math.log(4.0 * rv0 * rv0)
    sn[j + 1] = x
    imax = min(j + 1 + K, N - 1)
    g = np.zeros(N + 1)
    for i in range(j + 2, imax + 1):
        rv = (rn[i + 1] - rn[i - 1]) / (2.0 * h)
        rvv = (rn[i + 1] - 2.0 * rn[i] + rn[i - 1]) / (h * h)
        pv = (pn[i + 1] - pn[i - 1]) / (2.0 * h)
        g[i] = (rvv + 0.5 * rn[i] * pv * pv) / rv
    for i in range(j + 1, imax):
        sn[i + 1] = sn[i] + 0.5 * h * (g[i] + g[i + 1])
    for i in range(imax + 1, N + 1):
        rW = r[i - 1]; rA = rn[i - 1]; rB = r[i]; rC = rn[i]
        pW = p[i - 1]; pA = pn[i - 1]; pB = p[i]; pC = pn[i]
        sW = sg[i - 1]; sA = sn[i - 1]; sB = sg[i]
        rm = 0.25 * (rW + rA + rB + rC); pm = 0.25 * (pW + pA + pB + pC)
        sm = 0.25 * (sW + sA + sB + sn_guess[i])
        f = math.exp(sm)
        ru = ((rA - rW) + (rC - rB)) / (2.0 * h)
        rv = ((rB - rW) + (rC - rA)) / (2.0 * h)
        pu = ((pA - pW) + (pC - pB)) / (2.0 * h)
        pv = ((pB - pW) + (pC - pA)) / (2.0 * h)
        suv = (f + 4.0 * ru * rv) / (2.0 * rm * rm) - pu * pv                  # = [U f - 2 p_u p_v - 4 r_uv / r] / 2 with r_uv from its equation
        sn[i] = sA + sB - sW + h * h * suv
    return sn, G


@njit(cache=True)
def evolve(r0, p0, s0, h, A, kappa, s_, d_, nrows, snap_every, smax, stop_after_trap):
    """Evolve nrows rows. Returns snapshots and per-row diagnostics.
    diag columns: 0 u, 1 min r_v over interior, 2 r at that min, 3 max r_u, 4 max s, 5 number of trapped points (r_v<0 & r_u<0),
                  6 number of anti-trapped points (r_u>0 & r_v>0), 7 centre s, 8 centre p, 9 min r over valid points beyond centre, 10 max |C_uu|"""
    N = r0.shape[0] - 1
    nsnap = nrows // snap_every + 1
    Rs = np.zeros((nsnap, N + 1)); Ps = np.zeros((nsnap, N + 1)); Ss = np.zeros((nsnap, N + 1)); Us = np.zeros(nsnap)
    Ru = np.zeros((nsnap, N + 1))
    diag = np.zeros((nrows + 1, 11))
    r = r0.copy(); p = p0.copy(); sg = s0.copy()
    pv_center = (-3.0 * p[0] + 4.0 * p[1] - p[2]) / (2.0 * h)
    k_snap = 0
    last = 0
    trapped_seen = False
    rows_after = 0
    for j in range(nrows):
        # snapshot of the current row (also r_u by backward difference is added from the next row)
        rn = np.zeros(N + 1); pn = np.zeros(N + 1); sn = np.zeros(N + 1)
        pc = p[j]
        pvc = pv_center
        x_old = sg[j]; F_old = 0.0
        x = sg[j]
        for i in range(j + 1, N + 1):
            sn[i] = sg[i] + (x - sg[j])
        for it in range(10):
            rn[j + 1] = 0.0; pn[j + 1] = pc
            for i in range(j + 2, N + 1):
                # W = (j, i-1), A = (j+1, i-1), B = (j, i)
                rc, pcc = cell(r[i - 1], rn[i - 1], r[i], p[i - 1], pn[i - 1], p[i], sg[i - 1], sn[i - 1], sg[i], sn[i], h, A, kappa, s_, d_)
                rn[i] = rc; pn[i] = pcc
            sn_new, Gx = sigma_row(r, p, sg, rn, pn, sn, x, j, N, h, A, kappa, s_, d_, KCON)
            Fx = x - Gx                                                       # residual of sigma_c = G(sigma_c); solved by the secant method (gain of the plain iteration is not < 1)
            pc = (p[j] + h * pv_center + 0.5 * (4.0 * pn[j + 2] - pn[j + 3])) / 2.5          # implicit trapezoid for p_t = 2 p_v with the one-sided p_v (plain iteration has gain -1.5)
            pvc = (-3.0 * pc + 4.0 * pn[j + 2] - pn[j + 3]) / (2.0 * h)
            for i in range(j + 1, N + 1):
                sn[i] = sn_new[i]
            if abs(Fx) < 1e-12:
                break
            if it > 0 and abs(Fx - F_old) > 1e-10 * (abs(Fx) + abs(F_old)) + 1e-300:
                x_new = x - Fx * (x - x_old) / (Fx - F_old)
            else:
                x_new = x - 0.3 * Fx
            x_old = x; F_old = Fx
            dx = x_new - x
            x = x_new
            for i in range(j + 1, N + 1):                                     # shift the whole row guess by the change of the centre value
                sn[i] += dx
        pn[j + 1] = pc
        # diagnostics for row j (need r_u: forward difference between rows j and j+1 at the same i)
        mrv = 1e30; rmin_at = 0.0; mru = -1e30; smax_row = -1e30; ntrap = 0; nanti = 0; cmax = 0.0
        for i in range(j + 2, N):
            rv = (r[i + 1] - r[i - 1]) / (2.0 * h)
            ru = (rn[i] - r[i]) / h
            if rv < mrv:
                mrv = rv; rmin_at = r[i]
            if ru > mru:
                mru = ru
            if sg[i] > smax_row:
                smax_row = sg[i]
            if rv < 0 and ru < 0:
                ntrap += 1
            if ru > 0 and rv > 0:
                nanti += 1
        if j % snap_every == 0 and k_snap < nsnap:
            Rs[k_snap, :] = r; Ps[k_snap, :] = p; Ss[k_snap, :] = sg; Us[k_snap] = j * h
            for i in range(j + 1, N + 1):
                Ru[k_snap, i] = (rn[i] - r[i]) / h
            k_snap += 1
        diag[j, 0] = j * h; diag[j, 1] = mrv; diag[j, 2] = rmin_at; diag[j, 3] = mru; diag[j, 4] = smax_row
        diag[j, 5] = ntrap; diag[j, 6] = nanti; diag[j, 7] = sg[j]; diag[j, 8] = p[j]
        diag[j, 9] = r[j + 5] if j + 5 <= N else 0.0
        last = j
        r = rn; p = pn; sg = sn; pv_center = pvc
        if smax_row > smax:
            break
        if ntrap > 0:
            trapped_seen = True
        if trapped_seen:
            rows_after += 1
            if rows_after > stop_after_trap:
                break
        if N - (j + 4) < 8:
            break
    return Rs[:k_snap], Ps[:k_snap], Ss[:k_snap], Us[:k_snap], Ru[:k_snap], diag[:last + 1]


def initial_cone(A, pot, V_w, w, N, h, profile="tanh"):
    """Data on the outgoing cone u = 0: r = v/2, p = A psi(v), s from the vv constraint s_v = r p_v^2 (f = 1 at the centre)."""
    v = np.arange(N + 1) * h
    pf, pt = pot["psi_f"] * A, pot["psi_t"] * A
    p = pt + (pf - pt) * 0.5 * (1.0 - np.tanh((v - V_w) / w))
    r = v / 2.0
    pv = np.gradient(p, h, edge_order=2)
    integrand = r * pv ** 2
    s = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * h)])
    return v, r, p, s


def run(A=0.1, q=1.0, kappa=24.0, V_w=2.0, w=0.2, h=0.01, vmax=None, nrows=None, snap_every=10, smax=14.0, stop_after_trap=60, rho_in=3.0, extra=None):
    pot = make_potential(q, kappa, rho_in)
    if vmax is None:
        vmax = V_w + 10.0 * w + 2.0
    N = int(round(vmax / h))
    v, r0, p0, s0 = initial_cone(A, pot, V_w, w, N, h)
    if nrows is None:
        nrows = N - 8
    Rs, Ps, Ss, Us, Ru, diag = evolve(r0, p0, s0, h, A, kappa, pot["s"], pot["d"], nrows, snap_every, smax, stop_after_trap)
    return dict(pot=pot, A=A, V_w=V_w, w=w, h=h, v=v, R=Rs, P=Ps, S=Ss, U=Us, Ru=Ru, diag=diag, N=N)


if __name__ == "__main__":
    import time
    t0 = time.time()
    out = run(A=0.05, V_w=2.0, w=0.2, h=0.02)
    print("done", time.time() - t0, out["diag"].shape)
    d = out["diag"]
    print("last row u=%.3f, min r_v=%.4f, max s=%.3f, trapped pts=%d" % (d[-1, 0], d[:, 1].min(), d[:, 4].max(), d[:, 5].max()))
