# -*- coding: utf-8 -*-
"""FigD1: Phase 2D two-ended geometry around the stable cold cycle."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "phase2d_bridge.json")))
A = D["analysis"][0]
T = D["T"]
m, rh, kappa = A["m"], A["r_h"], A["kappa"]
Y = A["Y_diag"]
t = np.array(Y["t"]); a = np.array(Y["a"]); R = np.array(Y["R"]); s = np.array(Y["s"]); rs = np.array(Y["rstar"])
fig, ax = plt.subplots(2, 2, figsize=(11.5, 9))

# (a) schematic Penrose diagram of the kept spacetime (Kruskal regions I, II, III, IV)
P = ax[0, 0]
P.set_title("(a) Causal structure (schematic, not to scale)", fontsize=10)
pi2 = np.pi / 2
def XY(u, v): return (v - u), (u + v)             # u = atan U, v = atan V
# outline: I, III diamonds and the black-hole/white-hole triangles
poly = [(np.pi, 0), (pi2, pi2), (-pi2, pi2), (-np.pi, 0), (-pi2, -pi2), (pi2, -pi2), (np.pi, 0)]
P.plot([p_[0] for p_ in poly], [p_[1] for p_ in poly], color="0.55", lw=1)
P.plot([-pi2, pi2], [pi2, pi2], color="r", lw=2.2); P.plot([-pi2, pi2], [-pi2, -pi2], color="r", lw=2.2)
P.plot([-pi2, pi2], [pi2, -pi2], color="k", lw=1.1); P.plot([-pi2, pi2], [-pi2, pi2], color="k", lw=1.1)
ss = np.linspace(-3.0, 3.0, 300); r0 = 0.3
uS = np.arctan(np.exp(r0 + ss)); vS = -np.arctan(np.exp(r0 - ss))
xs, ys = XY(uS, vS)
P.plot(xs, ys, color="C0", lw=2.4, label="moving shell (periodic)")
P.fill_betweenx(ys, -np.pi * np.ones_like(xs) + np.abs(ys) * 0.0 - 0.0, xs, color="C0", alpha=0.12, label="closed cyclic region (replaces the left end)")
for k, s0 in enumerate(np.linspace(-1.5, 1.2, 7)):
    u0 = np.arctan(np.exp(r0 + s0)); v0 = -np.arctan(np.exp(r0 - s0))
    vend = pi2 - u0
    vv = np.linspace(v0, vend, 40)
    P.plot(*XY(u0 + 0 * vv, vv), color="C1", lw=0.9, label="radial null rays leaving the cap" if k == 0 else None)
P.text(1.9, 0.0, "I\n(parent\nuniverse)", ha="center", fontsize=8); P.text(0, 1.0, "II\n(black hole)", ha="center", fontsize=8)
P.text(0, -1.05, "IV (white hole)", ha="center", fontsize=8); P.text(-1.45, 0.0, "III", ha="center", fontsize=8)
P.annotate("throat: bifurcation sphere\nr_h = %.4f (hidden)" % rh, xy=(0, 0), xytext=(0.3, -0.45), fontsize=8, arrowprops=dict(arrowstyle="->"))
P.text(0.0, -1.75, "No causal curve joins the cap and region I", ha="center", fontsize=8, style="italic")
P.set_aspect("equal"); P.axis("off"); P.legend(loc="lower left", fontsize=7, frameon=False, bbox_to_anchor=(0.0, -0.17))

# (b) spheres along the slice, in time
B = ax[0, 1]
B.set_title("(b) Equator a(t), shell R(t) and throat r_h", fontsize=10)
B.semilogy(t / T, a, color="C3", label="equator of the cap (maximal sphere), a(t)")
B.semilogy(t / T, R, color="C0", label="shell R(t)")
B.axhline(rh, color="k", ls="--", label="throat radius r_h = %.4f" % rh)
B.set_xlabel("t / T"); B.set_ylabel("areal radius"); B.legend(fontsize=7, loc="center right")
B.text(0.02, 0.1, "ordering a > R > r_h at all times; shell/throat >= %.0f" % A["shell_over_throat_min"], transform=B.transAxes, fontsize=8)

# (c) region III in static coordinates: shell and rays
C = ax[1, 0]
C.set_title("(c) Region III, static coordinates: rays fall to the horizon", fontsize=10)
C.plot(s, np.log(R / rh), color="C0", lw=1.8, label="shell")
inv = lambda r_star: brentq(lambda Rr: Rr + m * np.log(abs(Rr / m - 1.0)) - r_star, m * (1 + 1e-12), 5.0)
for k, ray in enumerate(A["rays"][::3]):
    if not ray["crossed"]: continue
    s_c, rs_c = ray["s_cross"], ray["rstar_cross"]
    ss2 = np.linspace(s_c, s_c + 0.36, 120)
    Rr = []
    for z in ss2:
        try: Rr.append(inv(rs_c - (z - s_c)))
        except Exception: Rr.append(np.nan)
    C.plot(ss2, np.log(np.array(Rr) / rh), color="C1", lw=0.8, label="ray from the cap" if k == 0 else None)
C.set_xlabel("Killing time s"); C.set_ylabel("ln(R / r_h)"); C.legend(fontsize=7); C.set_xlim(0, 4.8); C.set_ylim(0, 5.0)

# (d) neck-map portrait around the elliptic fixed point
Dd = ax[1, 1]
Dd.set_title("(d) Neck-map portrait of the C/Cc = 0.35 cycle", fontsize=10)
for k, pts in enumerate(D["portrait"]):
    if pts:
        pts = np.array(pts); Dd.plot(pts[:, 0], pts[:, 1], ".", ms=2, color="C%d" % (k % 10))
Dd.plot([0], [0], "k+", ms=10)
Dd.set_xlabel("delta phi / v"); Dd.set_ylabel("phi_dot / v")
plt.tight_layout(); plt.savefig(os.path.join(HERE, "..", "figures", "FigD1_phase2d.png") if os.path.isdir(os.path.join(HERE, "..", "figures")) else os.path.join(HERE, "FigD1_phase2d.png"), dpi=160)
print("saved")
