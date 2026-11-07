# -*- coding: utf-8 -*-
"""Phase 2E: does a weaker dust-scalar coupling (larger rho_crit = rc V_f, M^2 = rc V_f / mu^2) give elliptic cold cycles with milder inhomogeneous growth?   python phase2e_coupling_scan.py <rc> <mu>"""
import json, os, sys
import numpy as np
from scipy.optimize import brentq
import phase2e_perturb_floquet as PF
from phase2d_cycle_validation import Cycle

HERE = os.path.dirname(os.path.abspath(__file__))
rc, mu = float(sys.argv[1]), float(sys.argv[2])
out = []
for x in np.round(np.arange(0.10, 0.96, 0.05), 3):
    c = Cycle(x, q=1.0, rc=rc, mu=mu)
    zs = np.geomspace(0.01, 0.9, 70)
    vals = [c.g_turn(z * c.v) for z in zs]
    for i in range(len(zs) - 1):
        if vals[i] is not None and vals[i + 1] is not None and vals[i] * vals[i + 1] < 0:
            try:
                z = brentq(lambda s: c.g_turn(s * c.v), zs[i], zs[i + 1], xtol=1e-13)
            except Exception:
                continue
            Jn = c.jac_neck(z * c.v, 0.0)
            if Jn is None: continue
            tr = float(np.trace(Jn)); det = float(np.linalg.det(Jn))
            rec = dict(rc=rc, mu=mu, x=float(x), phi0_over_v=float(z), trace=tr, det=det, elliptic=bool(abs(tr) < 2 and abs(det - 1) < 1e-3), period=float(c.flow_to_bounce(z * c.v, 0.0)[0]))
            if rec["elliptic"]:
                bg0 = np.array([c.neck_a(z * c.v, 0.0), 0.0, z * c.v, 0.0]); T = rec["period"]
                inh = {}
                for n in (3, 6, 12, 60):
                    M, _ = PF.Pert(c, n).monodromy(bg0, T)
                    inh[n] = float(max(abs(np.linalg.eigvals(M))))
                rec["inhomogeneous_max_multiplier"] = inh
                g2 = 6 * c.C * (z * c.v) ** 2 / (c.M2 * bg0[0] ** 3 * (2 * c.M2 + (z * c.v) ** 2))
                rec["fifth_force_rate_at_neck"] = float(np.sqrt(g2))
            out.append(rec)
            print("rc=%g mu=%g C/Cc=%.2f phi0/v=%.4f trace=%+.4f det=%.5f T=%.3f %s" % (rc, mu, x, z, tr, det, rec["period"], ("ELLIPTIC inh=%s rate=%.2f" % ({k: "%.3g" % v for k, v in rec["inhomogeneous_max_multiplier"].items()}, rec["fifth_force_rate_at_neck"])) if rec["elliptic"] else ""), flush=True)
json.dump(out, open(os.path.join(HERE, "phase2e_coupling_scan_rc%g_mu%g.json" % (rc, mu)), "w"), indent=1)
print("done")
