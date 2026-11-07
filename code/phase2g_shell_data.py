# -*- coding: utf-8 -*-
"""Phase 2G: (1) matter content (sigma(t), p(t), w = p/sigma) of representative comoving NEC shells around the stable cycle; (2) field continuity: range of the interior phi at the shell versus the
exterior vacuum value, and the thin domain-wall tension sigma_wall(phi_in) = int_{phi_in}^{v} sqrt(2 (V - V_vac)) dphi for comparison."""
import json, os
import numpy as np
from scipy.integrate import solve_ivp, quad
import bridge_junction as BJ
from phase2d_cycle_validation import Cycle

HERE = os.path.dirname(os.path.abspath(__file__))
c = Cycle(0.35)
p0 = c.orbit(0.2229, span=0.05)
bg0 = [c.neck_a(p0, 0.0), 0.0, p0, 0.0]
T = c.flow_to_bounce(p0, 0.0)[0]
sol = solve_ivp(c.f, (0, T), bg0, method="DOP853", rtol=1e-13, atol=1e-15, dense_output=True)
t = np.linspace(0, T, 4000); Y = sol.sol(t)
a, ad, ph = Y[0], Y[1], Y[2]
add = np.array([c.f(0, Y[:, i])[1] for i in range(len(t))])
out = dict(shells=[])
for chi_frac in (0.52, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95):
    chi0 = chi_frac * np.pi; q = 1.0
    sn = np.sin(chi0); R = a * sn
    m_max = float(np.min(R * (sn ** 2 + (ad * sn) ** 2 + q * R ** 2)))
    best = None
    for mf in np.linspace(0.02, 0.98, 49):
        sg, pr, Rr = BJ.shell(mf * m_max, q, chi0, a, ad, add)
        if np.all(sg > 0) and np.all(sg + pr >= 0):
            marg = float(np.min(sg + pr) / np.max(sg))
            if best is None or marg > best[0]: best = (marg, mf, sg, pr)
    if best is None:
        print("chi0=%.2fpi: no NEC comoving shell" % chi_frac); continue
    marg, mf, sg, pr = best
    w = pr / sg
    rec = dict(chi0_over_pi=chi_frac, m_over_mmax=float(mf), m=float(mf * m_max), sigma_range=[float(sg.min()), float(sg.max())], w_range=[float(w.min()), float(w.max())], nec_margin=marg, R_range=[float(R.min()), float(R.max())])
    out["shells"].append(rec)
    print("chi0=%.2fpi m/mmax=%.2f: sigma in [%.4f, %.4f], p/sigma in [%.4f, %.4f], NEC margin %.3f, R in [%.3f, %.3f]" % (chi_frac, mf, sg.min(), sg.max(), w.min(), w.max(), marg, R.min(), R.max()))
# field continuity (q = 0: exterior vacuum phi = v, V(v) = 0)
q0 = Cycle(0.35, q=0.0)
v = c.v   # note: Cycle(...,q=1) v; the q = 1 model is the one cycled; exterior vacuum value is its true minimum
Vvac = float(c.V(v))
phr = [float(ph.min()), float(ph.max())]
out["phi_range_over_v"] = [phr[0] / v, phr[1] / v]
wall = lambda pin: quad(lambda x: np.sqrt(max(2 * (c.V(x) - Vvac), 0.0)), pin, v)[0]
out["wall_tension"] = [float(wall(abs(phr[0]))), float(wall(abs(phr[1])))]
out["V_vac"] = Vvac
print("phi in [%.4f, %.4f] v ; exterior vacuum phi = v (V_vac = %.3f); thin-wall tension estimate %.3f-%.3f" % (phr[0] / v, phr[1] / v, Vvac, out["wall_tension"][0], out["wall_tension"][1]))
json.dump(out, open(os.path.join(HERE, "phase2g_shell_data.json"), "w"), indent=1)
