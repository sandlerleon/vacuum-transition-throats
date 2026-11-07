# -*- coding: utf-8 -*-
"""Verification with the ORIGINAL Phase 1 code (closed_cycle, closed_stability): an elliptic (linearly stable) exact cold cycle of the single-field model, and its extent in C."""
import json, os
import numpy as np
from scipy.optimize import brentq
import closed_cycle as CC, closed_orbits as CO, closed_stability as CS

HERE = os.path.dirname(os.path.abspath(__file__))


def orbits_at(q, rc, mu, x, lo=-3.0, hi=-0.1, n=300):
    m = CC.Model(q=q, rc=rc, mu=mu, C=x * CC.C_CRIT)
    roots, grid, vals = CO.find_orbits(m, lo=lo, hi=hi, n=n)
    return m, roots


if __name__ == "__main__":
    out = []
    for x in np.round(np.arange(0.30, 0.40, 0.005), 3):
        m, roots = orbits_at(1.0, 0.2, 20.0, x)
        for ls in roots:
            s = 10 ** ls * m.v
            try:
                J = CS.jacobian(m, s, 0.0)
                nm = CS.neck_map(m, s, 0.0)
            except Exception:
                J = None
            if J is None or nm is None:
                continue
            ev = np.linalg.eigvals(J)
            rec = dict(C_over_Ccrit=float(x), log10_seed_over_v=float(ls), phi0_over_v=float(10 ** ls), period=nm[2], return_phi_err=float(abs(nm[0] - s) / m.v),
                       trace=float(np.trace(J)), det=float(np.linalg.det(J)), multipliers_abs=[float(abs(e)) for e in ev], elliptic=bool(abs(np.trace(J)) < 2 and abs(np.linalg.det(J) - 1) < 1e-4))
            out.append(rec)
            print("C/Cc=%.3f phi0/v=%.5f T=%.4f trace=%+.6f det=%.6f |mu|=%s %s" % (x, 10 ** ls, nm[2], rec["trace"], rec["det"], ",".join("%.6f" % z for z in rec["multipliers_abs"]), "ELLIPTIC" if rec["elliptic"] else ""), flush=True)
    json.dump(out, open(os.path.join(HERE, "single_elliptic.json"), "w"), indent=1)
