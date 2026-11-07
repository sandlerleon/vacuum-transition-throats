# -*- coding: utf-8 -*-
"""Phase 3C.3: do the magnetized Bianchi-I bounce solutions have enough field to regularise the effective-fluid inertia?  MHD criterion at the bounce: rho_B > |p_spin| = B n^2 (perpendicular fast mode real when b^2 > 2 B n^2, b^2 = 2 rho_B => rho_B > B n^2)
and v_A^2 = b^2/(w + b^2) with w = n (1 - 2 B n)."""
import json, os
import numpy as np
import phase3c2_bianchi as P

HERE = os.path.dirname(os.path.abspath(__file__))
Bj = json.load(open(os.path.join(HERE, "phase3c2_bianchi.json")))
out = []
for r in Bj["critical"]:
    if r["rhoB0"] <= 0 or r["B_critical"] is None: continue
    for mult in (1.5, 3.0, 10.0):
        m = P.Model(B=mult * r["B_critical"], v0sq=r["v0sq"], rhoB0=r["rhoB0"])
        o = P.run(m, r["shear0"])
        if o["outcome"] != "bounce": continue
        ye = o["sol"].y_events[0][0]
        rho, px, py, parts = m.parts(ye[0], ye[1])
        V = np.exp(ye[0] + 2 * ye[1]); n = m.N0 / V
        w = n * (1 - 2 * m.B * n); pspin = m.B * n * n
        rec = dict(case=r["case"], B_over_Bc=mult, V_bounce=float(V), rho_B_over_pspin=float(parts["mag"] / pspin), w_over_n=float(w / n), n_over_nstar=float(n * 2 * m.B), vA2=float(2 * parts["mag"] / (w + 2 * parts["mag"])), perp_fast_real=bool(parts["mag"] > pspin))
        out.append(rec)
        print("%-44s B/Bc=%4.1f  V_b=%.3g  n/n*=%.3g  rho_B/|p_spin|=%.3g  v_A^2=%.3g  perpendicular fast mode real: %s" % (r["case"][:44], mult, V, rec["n_over_nstar"], rec["rho_B_over_pspin"], rec["vA2"], rec["perp_fast_real"]))
json.dump(out, open(os.path.join(HERE, "phase3c3_bounce_fields.json"), "w"), indent=1)
