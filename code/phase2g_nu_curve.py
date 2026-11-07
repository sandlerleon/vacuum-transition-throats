import json, os
import numpy as np
import phase2f_model as F
HERE = os.path.dirname(os.path.abspath(__file__))
m = F.Model2(0.35); cd = F.cycle_data(m, 0.2229)
out = []
for nu in [2.002, 2.01, 2.05, 2.15, 2.3, 2.5, 2.75, 3.0, 3.25, 3.5, 4.0, 5.0]:
    try:
        P = F.Pert2(m, 3); P.k2 = nu * nu - 1.0
        M = P.monodromy(cd["bg0"], cd["T"]); ev = np.linalg.eigvals(M)
        g = float(max(abs(ev))); gmin = float(min(abs(ev)))
        # also check constraint-propagation sanity: det
        out.append(dict(nu=nu, growth=g, shrink=gmin, det=float(np.linalg.det(M))))
        print("nu=%.3f (k^2=%.3f): max|mu|=%.4g  min|mu|=%.4g det=%.5f" % (nu, nu * nu - 1, g, gmin, np.linalg.det(M)), flush=True)
    except Exception as e:
        print("nu=%.3f: failed %s" % (nu, e)); out.append(dict(nu=nu, failed=str(e)))
json.dump(out, open(os.path.join(HERE, "phase2g_nu_curve.json"), "w"), indent=1)
