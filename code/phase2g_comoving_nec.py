# -*- coding: utf-8 -*-
"""Phase 2G-A, step 0: which shell is consistent with a Neumann (no-flux) condition for the scalar field at the shell?  For a homogeneous interior, d_n phi = sinh(beta) phi_dot, so a no-flux
condition forces beta = 0 (a COMOVING shell).  Scan comoving shells around the stable cycle (C/Cc = 0.35): sigma > 0 and sigma + p >= 0 over the whole cycle?"""
import json, os
import numpy as np
from scipy.integrate import solve_ivp
import bridge_junction as BJ
from phase2d_cycle_validation import Cycle

HERE = os.path.dirname(os.path.abspath(__file__))
c = Cycle(0.35)
p = c.orbit(0.2229, span=0.05)
bg0 = [c.neck_a(p, 0.0), 0.0, p, 0.0]
T = c.flow_to_bounce(p, 0.0)[0]
sol = solve_ivp(c.f, (0, T), bg0, method="DOP853", rtol=1e-13, atol=1e-15, dense_output=True)
t = np.linspace(0, T, 4000); Y = sol.sol(t)
a, ad = Y[0], Y[1]
add = np.array([c.f(0, Y[:, i])[1] for i in range(len(t))])
rows = []
best = None
for q in (1.0, 0.0, -0.5, -0.9):
    for chi_frac in np.linspace(0.505, 0.995, 99):
        chi0 = chi_frac * np.pi
        sn = np.sin(chi0)
        R = a * sn
        m_max = float(np.min(R * (sn ** 2 + (ad * sn) ** 2 + q * R ** 2)))
        if m_max <= 0: continue
        for mf in np.linspace(0.02, 0.98, 25):
            res = BJ.shell(mf * m_max, q, chi0, a, ad, add)
            if res is None: continue
            sg, pr, Rr = res
            ok_s = bool(np.all(sg > 0)); nec = float(np.mean((sg > 0) & (sg + pr >= 0)))
            rows.append((q, chi_frac, mf, ok_s, nec, float(np.min(sg + pr) / np.max(sg))))
            if ok_s and (best is None or rows[-1][5] > best[5]): best = rows[-1]
full = [r for r in rows if r[3] and r[4] == 1.0]
print("combinations tried:", len(rows), "| sigma>0 all cycle:", sum(1 for r in rows if r[3]), "| sigma>0 AND NEC all cycle:", len(full))
print("best margin min(sigma+p)/max(sigma) among sigma>0:", best)
for q in (1.0, 0.0, -0.5, -0.9):
    qq = [r for r in rows if r[0] == q and r[3]]
    if qq: print("  q=%g: sigma>0 combos %d, max NEC fraction %.3f, best margin %.4f" % (q, len(qq), max(r[4] for r in qq), max(r[5] for r in qq)))
json.dump(dict(n_tried=len(rows), n_sigma_pos=sum(1 for r in rows if r[3]), n_nec_full=len(full), best=best, nec_full=full[:20]), open(os.path.join(HERE, "phase2g_comoving_nec.json"), "w"), indent=1)
