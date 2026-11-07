# -*- coding: utf-8 -*-
"""Phase 2G-A: perturbation spectrum on the actual cap [0, chi0] with a COMOVING shell (consistent with zero scalar and dust flux through the shell).
Radial modes Q_nu = sin(nu chi)/sin(chi) (regular at the pole); boundary conditions at the shell from the no-flux requirement: d_chi(delta phi) = 0 and dust normal velocity U = v Q' = 0, i.e. Q_nu'(chi0) = 0
  <=>  tan(nu chi0) = nu tan(chi0).   k^2 = nu^2 - 1 (nu = 1 is the homogeneous mode).  Floquet multipliers of each cap mode over one cycle (generalised equations of phase2f_model).
For chi0 = pi the roots are the integers (the closed three-sphere of Phase 2E)."""
import json, os
import numpy as np
from scipy.optimize import brentq
import phase2f_model as F

HERE = os.path.dirname(os.path.abspath(__file__))


def cap_roots(chi0, numax):
    f = lambda nu: nu * np.cos(nu * chi0) * np.sin(chi0) - np.sin(nu * chi0) * np.cos(chi0)
    grid = np.linspace(1.0005, numax, int(numax * 400))
    v = [f(x) for x in grid]
    roots = []
    for i in range(len(grid) - 1):
        if v[i] * v[i + 1] < 0:
            roots.append(brentq(f, grid[i], grid[i + 1], xtol=1e-13))
    return roots


m = F.Model2(0.35)
cd = F.cycle_data(m, 0.2229)
print("cycle: T=%.6f trace=%.5f" % (cd["T"], cd["trace"]))
nec = json.load(open(os.path.join(HERE, "phase2g_comoving_nec.json")))
out = []
for chi_frac in (0.55, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99):
    chi0 = chi_frac * np.pi
    roots = cap_roots(chi0, 40.0)
    rec = dict(chi0_over_pi=chi_frac, nu=roots, growth=[])
    for nu in roots:
        if abs(nu - 2.0) < 2e-3:
            rec["growth"].append(None); continue
        P = F.Pert2(m, 3); P.k2 = nu * nu - 1.0
        M = P.monodromy(cd["bg0"], cd["T"])
        rec["growth"].append(float(max(abs(np.linalg.eigvals(M)))))
    out.append(rec)
    g = [x for x in rec["growth"] if x is not None]
    print("chi0=%.2f pi: %d modes (nu<=40); lowest nu = %s; growth of lowest 6 = %s; max growth %.3g; min growth %.3g" % (chi_frac, len(roots), [round(x, 3) for x in roots[:6]],
          [None if x is None else round(x, 2) for x in rec["growth"][:6]], max(g), min(g)), flush=True)
    json.dump(out, open(os.path.join(HERE, "phase2g_cap_spectrum.json"), "w"), indent=1)
