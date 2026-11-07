# -*- coding: utf-8 -*-
"""Classical bounce of a closed universe with a false-vacuum (de Sitter) energy and dust: the 'neck' of global de Sitter space as a cold, null-energy-respecting
bridge in TIME (not in space), and what the observed spatial curvature does to it.   8 pi G = 1, c = 1.

Closed FRW (curvature radius a):  adot^2 = f(a) = V_f a^2/3 + C/a - 1,   C = rho_m a^3 / 3  (conserved dust parameter, = 2GM in length),   R_dS = sqrt(3/V_f).
A bounce at a_min > 0 on the branch that comes in from large a exists iff f has two positive roots a1 < a2, i.e. iff C < C_c = 2 R_dS / (3 sqrt 3)  (the Eddington-Lemaitre
condition 9 Lambda G^2 M^2 < 1); the bounce is at the larger root a2, with density rho_neck = 3 / a2^2 (<= a few V_f): a cold bounce, and it needs the strong energy condition
to fail (vacuum energy), not the null energy condition.

    python closed_bounce.py          (writes closed_bounce.json)
"""
import json
import os

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
C_CRIT = 2.0 / (3.0 * np.sqrt(3.0))            # in units of R_dS


def f(a, C):
    """adot^2 in units R_dS = 1 (V_f = 3)."""
    return a * a + C / a - 1.0


def bounce_radius(C):
    """Outer root a2 of f (R_dS = 1); None if C >= C_c."""
    if C >= C_CRIT:
        return None
    a_min = (C / 2.0) ** (1.0 / 3.0)
    return brentq(lambda a: f(a, C), a_min, 10.0)


def observed_ratio(omega_k, omega_m=0.3110, omega_l=0.6889):
    """C / C_c for today's matter in a closed universe of spatial curvature Omega_k:  C = Omega_m / (H0 Omega_k^{3/2}),  R_dS = 1 / (H0 sqrt(Omega_L))."""
    return omega_m * np.sqrt(omega_l) / (omega_k ** 1.5) / C_CRIT


def main():
    out = {"C_crit_over_RdS": float(C_CRIT)}
    rows = []
    for x in (0.0, 0.1, 0.3, 0.5, 0.8, 0.95, 0.99):
        C = x * C_CRIT
        a2 = bounce_radius(C) if x > 0 else 1.0
        rows.append(dict(C_over_Ccrit=x, a_bounce_over_RdS=float(a2), rho_neck_over_Vf=float(1.0 / a2 ** 2)))
    out["bounce_table"] = rows
    out["observed"] = [dict(omega_k=k, C_over_Ccrit=float(observed_ratio(k))) for k in (0.001, 0.004, 0.01, 0.1, 0.5, 0.8)]
    out["omega_k_for_critical"] = float(brentq(lambda k: observed_ratio(k) - 1.0, 0.01, 5.0))
    json.dump(out, open(os.path.join(HERE, "closed_bounce.json"), "w"), indent=1)
    print("C_c = %.4f R_dS" % C_CRIT)
    for r in rows:
        print("C/C_c %.2f: bounce at a = %.3f R_dS, rho_neck = %.3f V_f" % (r["C_over_Ccrit"], r["a_bounce_over_RdS"], r["rho_neck_over_Vf"]))
    for o in out["observed"]:
        print("Omega_k = %.3f: C/C_c = %.3g" % (o["omega_k"], o["C_over_Ccrit"]))
    print("Omega_k at which C = C_c: %.3f" % out["omega_k_for_critical"])


if __name__ == "__main__":
    main()
