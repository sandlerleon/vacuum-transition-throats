# -*- coding: utf-8 -*-
"""FigE1: (a) inhomogeneous linear growth per cycle versus radial harmonic n; (b) homogeneous stability window of the cold cycle."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
F = json.load(open(os.path.join(HERE, "phase2e_perturb_floquet.json")))
M = json.load(open(os.path.join(HERE, "phase2e_perturb_mechanism.json")))
V = json.load(open(os.path.join(HERE, "phase2d_cycle_validation.json")))
W = json.load(open(os.path.join(HERE, "phase2e_window.json")))
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6))
a = ax[0]
ns = [r["n"] for r in F["modes"]]; mx = [r["max_abs"] for r in F["modes"]]
a.loglog(ns, mx, "o-", color="C3", label="full system (scalar + coupled dust + gravity)")
cn = [r["n"] for r in M["control"]]; cm = [max(r["multipliers_abs"]) for r in M["control"]]
a.loglog(cn, cm, "s--", color="C2", label="control: dust-scalar coupling off in the perturbations")
a.axhline(np.exp(M["estimate_log_growth"]), color="C1", ls=":", label="quasi-static fifth-force estimate, e^%.1f" % M["estimate_log_growth"])
a.axhline(1.0, color="k", lw=0.8, label="homogeneous cycle (n = 1): |mu| = 1")
a.set_xticks([3, 5, 10, 20, 40, 60]); a.set_xticklabels(["3", "5", "10", "20", "40", "60"]); a.minorticks_off(); a.set_xlabel("radial harmonic n (k^2 = n^2 - 1)"); a.set_ylabel("largest Floquet multiplier per cycle")
a.set_title("(a) Inhomogeneous linear stability of the cycle", fontsize=10); a.legend(fontsize=7, loc="lower right")
b = ax[1]
br = [r for r in V["branch"] if r["found"]]
xs = [r["x"] for r in br]; tr = [r["trace"] for r in br]
b.plot(xs, tr, "o-", color="C0", ms=4, label="neck-map trace along the main branch")
for r in W["points"]:
    b.plot([r["x"]], [r["trace"]], "ko", ms=4)
    b.annotate(r["name"], (r["x"], r["trace"]), textcoords="offset points", xytext=(4, 6), fontsize=7)
b.axhspan(-2, 2, color="C2", alpha=0.12, label="elliptic: |trace| < 2")
b.axvline(V["flip"]["C_over_Ccrit_star"], color="r", ls="--", lw=0.9, label="trace = -2 at C/Cc = %.4f" % V["flip"]["C_over_Ccrit_star"])
b.axvline(0.3600, color="gray", ls="--", lw=0.9, label="neck ceases to exist (~0.360)")
b.set_xlim(0.33, 0.364); b.set_ylim(-4.2, 0.6)
b.set_xlabel("C / Cc"); b.set_ylabel("trace of the neck-map Jacobian"); b.set_title("(b) Homogeneous stability window", fontsize=10); b.legend(fontsize=7, loc="lower right")
plt.tight_layout(); plt.savefig(os.path.join(HERE, "..", "figures", "FigE1_phase2e.png"), dpi=160); print("saved")
