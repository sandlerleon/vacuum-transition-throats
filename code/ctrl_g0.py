import json, numpy as np
from scipy.optimize import brentq
import bridge_reservoir as BR, bridge_reservoir_followup as BF, closed_cycle as CC
m0 = BF.model(20.0, 0.0)
c0 = 0.04522531
# 1-D scan of phi_dot at the first turnaround vs phi0 for the decoupled reservoir with the same chi0
phis = np.geomspace(0.01, 0.6, 60)
vals = []
for pf in phis:
    r = BR.residual(m0, pf * m0.v, c0)
    vals.append(np.nan if r is None else r[0])
roots = []
for i in range(len(phis) - 1):
    if np.isfinite(vals[i]) and np.isfinite(vals[i+1]) and vals[i] * vals[i+1] < 0:
        f = lambda x: BR.residual(m0, x * m0.v, c0)[0]
        try:
            roots.append(brentq(f, phis[i], phis[i+1], xtol=1e-13))
        except Exception: pass
print("g=0, chi0=%.5f: orbits (phi_dot=0 at turnaround) at phi0/v =" % c0, roots)
for pf in roots:
    fl = BR.floquet(m0, (pf * m0.v, 0.0, c0, 0.0), h=1e-6)
    if fl is None:
        print("  phi0/v=%.5f: Floquet failed (return map not defined for a perturbation)" % pf); continue
    print("  phi0/v=%.5f T=%.5f return err %.1e |mult| %s" % (pf, fl["period"], fl["state_return_error"], ",".join("%.5g" % z for z in fl["multipliers_abs"])))
# also pure single field (chi0 = 0) orbits across the same phi range
print("chi0 = 0 (no reservoir):")
vals0 = []
for pf in phis:
    r = BR.residual(m0, pf * m0.v, 0.0); vals0.append(np.nan if r is None else r[0])
for i in range(len(phis) - 1):
    if np.isfinite(vals0[i]) and np.isfinite(vals0[i+1]) and vals0[i] * vals0[i+1] < 0:
        f = lambda x: BR.residual(m0, x * m0.v, 0.0)[0]
        pf = brentq(f, phis[i], phis[i+1], xtol=1e-13)
        fl = BR.floquet(m0, (pf * m0.v, 0.0, 0.0, 0.0), h=1e-6)
        if fl is None:
            print("  phi0/v=%.5f: Floquet failed" % pf); continue
        print("  phi0/v=%.5f T=%.5f |mult| %s" % (pf, fl["period"], ",".join("%.5g" % z for z in fl["multipliers_abs"])))
