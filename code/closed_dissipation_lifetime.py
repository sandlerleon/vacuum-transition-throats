# -*- coding: utf-8 -*-
"""Phase 1B (continued): lifetime of the cold cycle (number of consecutive necks) against the dissipation strength, over a broad range of seeds (no fine tuning),
for the small-excursion orbit family (mu = 20, rho_crit = 0.2 V_f, C = 0.3 C_c) and a second family.  Radiation is evolved and enters the Friedmann equation.
Writes closed_dissipation_lifetime.json."""
import json, os, sys
import numpy as np
import closed_cycle as CC, closed_dissipation as CD

HERE = os.path.dirname(os.path.abspath(__file__))
FAMS = [dict(q=1.0, rc=0.2, mu=20.0, x=0.3), dict(q=1.0, rc=0.4, mu=10.0, x=0.6)]
GAMMAS = (0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0)
out = []
for fam in FAMS:
    for g in GAMMAS:
        m = CC.Model(q=fam["q"], rc=fam["rc"], mu=fam["mu"], C=fam["x"] * CC.C_CRIT, ups=g)
        life = []
        best = None
        for ls in np.linspace(-2.6, -0.3, 300):
            y0 = CD.start_state(m, 10 ** ls * m.v)
            if y0 is None:
                continue
            B, T, end, res, sol = CD.necks(m, y0, 150.0, rtol=1e-11)
            life.append(len(B))
            if best is None or len(B) > best[0]:
                best = (len(B), float(ls), [float(b[1][4] * b[1][0] ** 4) for b in B], end, float(sol.t[-1]))
        life = np.array(life)
        rec = dict(family=fam, gamma=g, n_seeds=int(len(life)), max_bounces=int(life.max()), frac_ge3=float(np.mean(life >= 3)), frac_ge5=float(np.mean(life >= 5)), frac_ge10=float(np.mean(life >= 10)),
                   median=float(np.median(life)), best_log10_seed=best[1], best_S=best[2][:200], best_end=best[3], best_t_end=best[4])
        out.append(rec)
        print("mu=%g rc=%g x=%.1f G=%-5g | max %3d bounces | frac>=3 %.2f >=5 %.2f >=10 %.2f | best seed %.3f ends %s at t=%.1f" %
              (fam["mu"], fam["rc"], fam["x"], g, rec["max_bounces"], rec["frac_ge3"], rec["frac_ge5"], rec["frac_ge10"], rec["best_log10_seed"], rec["best_end"], rec["best_t_end"]), flush=True)
json.dump(out, open(os.path.join(HERE, "closed_dissipation_lifetime.json"), "w"), indent=1)
