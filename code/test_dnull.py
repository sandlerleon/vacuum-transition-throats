# -*- coding: utf-8 -*-
"""Validation of dnull.py on exact solutions: pure de Sitter (rho = 3) and pure anti-de Sitter (rho = -3) with a regular centre.
Misner-Sharp mass must equal rho r^3 / 3 everywhere; the error must fall as h^2."""
import sys
import numpy as np
import dnull as D


def ms_error(vac, h, urows=1.2):
    """vac = +1: false vacuum everywhere (psi_f), -1: true vacuum everywhere."""
    pot = D.make_potential(1.0, 24.0, 3.0)
    A = 0.1
    vmax = 3.0
    N = int(round(vmax / h))
    v = np.arange(N + 1) * h
    p0 = np.full(N + 1, (pot["psi_f"] if vac > 0 else pot["psi_t"]) * A)
    r0 = v / 2.0
    s0 = np.zeros(N + 1)
    nrows = int(urows / h)
    Rs, Ps, Ss, Us, Ru, diag = D.evolve(r0, p0, s0, h, A, 24.0, pot["s"], pot["d"], nrows, 1, 50.0, 10 ** 9)
    rho = 3.0 * vac
    err = 0.0
    for k in range(Rs.shape[0] - 1):
        j = int(round(Us[k] / h))
        r = 0.5 * (Rs[k] + Rs[k + 1]); s = 0.5 * (Ss[k] + Ss[k + 1])        # row midpoints (u_j + h/2)
        rv_ = 0.5 * ((Rs[k][2:] - Rs[k][:-2]) + (Rs[k + 1][2:] - Rs[k + 1][:-2])) / (2 * h)
        for i in range(j + 6, min(N - 3, j + 200)):
            rv = rv_[i - 1]
            ru = Ru[k, i]
            f = np.exp(s[i])
            m = r[i] * (1 + 4 * ru * rv / f)
            err = max(err, abs(m - rho * r[i] ** 3 / 3.0) / (1 + abs(rho * r[i] ** 3 / 3.0)))
    return err


if __name__ == "__main__":
    for vac in (+1, -1):
        e = [ms_error(vac, h) for h in (0.04, 0.02, 0.01)]
        print("vacuum", "dS" if vac > 0 else "AdS", "max rel. Misner-Sharp error:", ["%.2e" % x for x in e],
              "orders:", ["%.2f" % np.log2(e[i] / e[i + 1]) for i in range(2)])
