# -*- coding: utf-8 -*-
"""Figure of the Phase 1B results: lifetime of the cold cycle against the dissipation strength, and growth of the entropy parameter S = rho_r a^4 per bounce."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
L = json.load(open(os.path.join(HERE, "closed_dissipation_lifetime.json")))
fine = json.load(open(os.path.join(HERE, "closed_dissipation_fine.json")))
fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.6))
for fam, col in (({"mu": 20.0, "rc": 0.2}, "#1f77b4"), ({"mu": 10.0, "rc": 0.4}, "#d62728")):
    rows = [r for r in L if r["family"]["mu"] == fam["mu"] and r["family"]["rc"] == fam["rc"]]
    g = [r["gamma"] for r in rows]
    lab = "mu=%g, rho_crit=%g V_f" % (fam["mu"], fam["rc"])
    ax[0].plot(g, [r["max_bounces"] for r in rows], "o-", color=col, label=lab)
    ax[1].plot(g, [100 * r["frac_ge3"] for r in rows], "s-", color=col)
ax[0].set_xlabel("dissipation G (H_Lambda)"); ax[0].set_ylabel("most bounces over 300 seeds"); ax[0].legend(fontsize=7)
ax[1].set_xlabel("dissipation G (H_Lambda)"); ax[1].set_ylabel("% of seeds with >= 3 bounces")
for r in fine["rows"]:
    if r["S"] and r["gamma"] in (0.001, 0.01, 0.1, 0.2):
        ax[2].plot(range(1, len(r["S"]) + 1), r["S"], "o-", label="G=%g" % r["gamma"])
ax[2].set_yscale("log"); ax[2].set_xlabel("bounce number"); ax[2].set_ylabel("S = rho_r a^4 (R_dS units)"); ax[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "..", "figures", "FigC2_dissipation.png"), dpi=150)
print("ok")
