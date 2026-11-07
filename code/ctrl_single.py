import json, numpy as np
from scipy.optimize import brentq
import bridge_reservoir as BR, closed_cycle as CC
B = json.load(open('closed_stability_breather.json')); o = [b for b in B if b['mu']==20.0 and b['rc']==0.2][0]
res = []
for x in (0.32, 0.33, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44):
    m = BR.TwoField(o['q'], o['rc'], o['mu'], x*CC.C_CRIT, 20.0, 0.0)
    phis = np.geomspace(0.005, 0.7, 90)
    vals = []
    for pf in phis:
        r = BR.residual(m, pf*m.v, 0.0); vals.append(np.nan if r is None else r[0])
    for i in range(len(phis)-1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i+1]) and vals[i]*vals[i+1] < 0:
            f = lambda z: BR.residual(m, z*m.v, 0.0)[0]
            try: pf = brentq(f, phis[i], phis[i+1], xtol=1e-13)
            except Exception: continue
            fl = BR.floquet(m, (pf*m.v, 0.0, 0.0, 0.0), h=1e-6)
            if fl is None: print("C/Cc=%.2f phi0/v=%.4f Floquet failed" % (x, pf)); continue
            mult = sorted(fl['multipliers_abs'], reverse=True)
            print("C/Cc=%.2f phi0/v=%.4f T=%.4f |mult| (phi-sector) %s" % (x, pf, fl['period'], ",".join("%.5g" % z for z in mult)), flush=True)
            res.append(dict(C_over_Ccrit=x, phi0_over_v=float(pf), period=fl['period'], multipliers_abs=mult))
json.dump(res, open('ctrl_single_fine.json','w'), indent=1)
