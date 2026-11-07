# -*- coding: utf-8 -*-
"""Phase 2F: worst-case inhomogeneous growth max_n G(n) (n = 3 .. 3 mu) along the fixed-trace mu-family, to see whether any mu gives bounded inhomogeneous perturbations."""
import json, os
import numpy as np
import phase2f_model as F

HERE = os.path.dirname(os.path.abspath(__file__))
fam = json.load(open(os.path.join(HERE, "phase2f_muscan_up_rc0.2_valid.json")))
sel = [fam[0]] + [r for r in fam if r["mu"] > t][:1] if False else None
targets = [20, 40, 80, 160, 230]
out = []
for t in targets:
    r = min(fam, key=lambda q: abs(q["mu"] - t))
    m = F.Model2(r["x"], mu=r["mu"], rc=0.2)
    cd = F.cycle_data(m, r["phi0_over_v"], span=0.04)
    ns = sorted(set([3, 4, 5, 6, 8, 10, 12, 16, 20, 30, 45, 60, 90, 120, 180, 240, 360, 480, 600]))
    ns = [n for n in ns if n <= 3 * r["mu"] + 10]
    gr = F.growth(m, cd, tuple(ns))
    worst = max(gr.items(), key=lambda kv: kv[1])
    out.append(dict(mu=r["mu"], T=cd["T"], growth={str(k): v for k, v in gr.items()}, worst_n=int(worst[0]), worst_G=float(worst[1]), min_G=float(min(gr.values()))))
    print("mu=%.0f: worst G = %.3g at n=%d; min over n = %.3g; G(n)= %s" % (r["mu"], worst[1], worst[0], min(gr.values()), ", ".join("%d:%.2g" % (k, v) for k, v in gr.items())), flush=True)
json.dump(out, open(os.path.join(HERE, "phase2f_worstcase.json"), "w"), indent=1)
