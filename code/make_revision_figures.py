# -*- coding: utf-8 -*-
"""Figures added in revision 2: Fig5_causal.png (causal classification: what the numerics establish and what is conjecture) and Fig7_wall_extrapolation.png.   python make_revision_figures.py"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, FancyArrowPatch

import dnull as D

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(os.path.dirname(HERE), "figures")
WALL = json.load(open(os.path.join(HERE, "wall_extrapolation.json")))
plt.rcParams.update({"font.size": 8})

# ---------------------------------------------------------------- Fig 5
o = D.run(A=0.1, V_w=4.0, w=0.1, h=0.01, snap_every=2, smax=12.0, stop_after_trap=20)
R, Ru, P, h, N = o["R"], o["Ru"], o["P"], o["h"], o["N"]
nrow = R.shape[0]; pot = o["pot"]; mid = 0.5 * 0.1 * (pot["psi_f"] + pot["psi_t"])
img = np.full((nrow, N + 1), np.nan); wall = np.full(nrow, np.nan); last = 0
for k in range(nrow):
    j = 2 * k; r = R[k]
    if not np.isfinite(r[j + 1:]).all():
        break
    last = k
    for i in range(j + 3, N - 2):
        rv = (r[i + 1] - r[i - 1]) / (2 * h); ru = Ru[k, i]
        img[k, i] = 0 if (ru < 0 and rv > 0) else (1 if (ru < 0 and rv < 0) else (2 if (ru > 0 and rv > 0) else 3))
    seg = np.arange(j + 3, N - 2); s = np.sign(P[k][seg] - mid); idx = np.where(s[:-1] * s[1:] < 0)[0]
    if idx.size:
        wall[k] = (seg[idx[0]] + 0.5) * h
u = o["U"]; v = np.arange(N + 1) * h
fig, ax = plt.subplots(1, 2, figsize=(9.4, 4.1), gridspec_kw=dict(width_ratios=[1.15, 1.0]))
a = ax[0]
cmap = ListedColormap(["#dddddd", "#d62728", "#1f77b4", "#ff7f0e"])
a.pcolormesh(v, u[:last + 1], img[:last + 1], cmap=cmap, vmin=-0.5, vmax=3.5, shading="auto")
a.plot(wall[:last + 1], u[:last + 1], "k--", lw=1.0)
ul = u[last]
a.axhline(ul, color="k", lw=1.2, ls=":")
a.axhspan(ul, ul * 1.55, color="#eeeeee", zorder=0)
a.text(0.35, ul * 1.28, "not computed (code stops here)", fontsize=7.5, color="0.25", va="center")
a.annotate("", xy=(2.9, 0.045), xytext=(2.3, 0.0), arrowprops=dict(arrowstyle="-|>", color="k", lw=0.9)); a.text(2.35, 0.055, "future\ncone", fontsize=6.5)
a.set_xlim(0, v[-1]); a.set_ylim(0, ul * 1.55); a.set_xlabel("v"); a.set_ylabel("u")
a.set_title("(a) numerical classification, $R_w=2R_{dS}$ (evolved domain only)", fontsize=8.5)
a.legend(handles=[Patch(color="#dddddd", label="normal"), Patch(color="#d62728", label="trapped"), Patch(color="#1f77b4", label="anti-trapped"), Patch(color="#ff7f0e", label="$r_u>0>r_v$"),
                  Patch(color="#eeeeee", label="no information")], fontsize=6.5, loc="lower right")
a.text(0.45, 0.02, "core (false vacuum)", fontsize=7, color="k")
a.text(4.85, 0.02, "exterior", fontsize=7)
# panel b: schematic with established / conjectured parts
b = ax[1]; b.set_xlim(0, 10); b.set_ylim(0, 10); b.axis("off")
b.set_title("(b) what is established, what is conjectured", fontsize=8.5)
b.fill([0.8, 5.0, 5.0, 0.8], [0.6, 0.6, 4.6, 4.6], color="#cfe0f3", alpha=0.9)
b.fill([5.0, 9.2, 9.2, 5.0], [0.6, 0.6, 4.6, 4.6], color="#f5c6c6", alpha=0.9)
b.plot([5.0, 5.0], [0.6, 4.6], "k--", lw=1.2); b.text(5.0, 4.3, "wall", fontsize=7, ha="center")
b.text(2.9, 2.6, "anti-trapped\nfalse-vacuum core\n(persists in the runs)", ha="center", fontsize=7)
b.text(7.1, 2.6, "exterior:\nnormal, then trapped", ha="center", fontsize=7)
b.plot([0.8, 9.2], [4.6, 4.6], "k:", lw=1.3); b.text(0.9, 4.8, "last evolved row (singularity not reached)", fontsize=7)
b.fill([0.8, 9.2, 9.2, 0.8], [4.6, 4.6, 9.4, 9.4], color="#eeeeee")
b.text(5.0, 7.6, "CONJECTURE, not computed:\nglobal topology of the continuation\n(second asymptotic region = bridge,\nor closes on itself); no slice through a throat;\nno daughter universe established", ha="center", va="center", fontsize=7.2,
       bbox=dict(boxstyle="round", fc="white", ec="0.5"))
b.annotate("", xy=(8.4, 3.55), xytext=(2.0, 3.55), arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.0, ls="--"))
b.text(5.2, 3.75, "core can signal the exterior", fontsize=6.8, color="#333333", ha="center")
b.annotate("", xy=(2.0, 1.6), xytext=(8.4, 1.6), arrowprops=dict(arrowstyle="-|>", color="#b00000", lw=1.0))
b.plot([5.0], [1.6], "x", color="#b00000", ms=8); b.text(5.2, 1.15, "no exterior-to-core causal curve ($du\\,dv\\geq0$)", fontsize=6.8, color="#b00000", ha="center")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig5_causal.png"), dpi=170); plt.close(fig)

# ---------------------------------------------------------------- Fig 7
fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.5))
pts = WALL["deviation_points"]; sl = WALL["deviation_slope_vs_w"]
wv = np.array([p["w"] for p in pts]); dv = np.array([p["dev"] for p in pts]); jit = np.random.default_rng(1).normal(0, 0.004, len(wv))
ax[0].plot(wv + jit, dv, "o", color="#1f77b4", ms=3.5, alpha=0.7, label="runs with a trapped surface and a normal core")
xx = np.linspace(0, 0.22, 50); ax[0].plot(xx, sl["slope_median"] * xx + sl["intercept_at_w0_median"], "k-", lw=1.0, label="linear fit, extrapolated to w → 0")
ax[0].axhspan(-sl["resolution_floor_median"], sl["resolution_floor_median"], color="0.85", label="resolution floor ±h/r_h (median)")
ax[0].axhline(0, color="k", lw=0.5); ax[0].set_xlim(0, 0.22)
ax[0].set_xlabel("wall width w / R_dS (physical: 1e-30)"); ax[0].set_ylabel("r_trapped / r_h(m_ext) − 1"); ax[0].legend(fontsize=6.5, loc="lower left")
ax[0].set_title("(a) horizon-radius test against wall width", fontsize=8.5)
byA = WALL["by_amplitude"]
sg = np.logspace(-3, 0.1, 200); series = sg ** 2 * (0.25 - 1.0 / 12.0) / 2.0
ax[1].loglog(sg, series, "k-", lw=1.0, label="$|m_c/m_c^{(0)}-1|$ from Eq. (9), $q=1$")
ax[1].loglog([byA[k]["sigma_hat_sim"] for k in byA], [abs(byA[k]["m_c_shift"]) for k in byA], "o", color="#1f77b4", label="simulated tensions (exact maximum)")
ax[1].axvline(2.709e-30, color="#d62728", lw=1.0)
ax[1].set_xlim(1e-4, 2); ax[1].set_ylim(1e-8, 1)
ax[1].text(1.3e-4, 2e-8, "physical tension $2.7\\times10^{-30}$:\ncorrection $6\\times10^{-61}$ (off scale)", fontsize=7, color="#d62728")
ax[1].set_xlabel("$\\hat\\sigma=\\sigma R_{dS}/M_{Pl}^2$"); ax[1].set_ylabel("relative shift of the critical mass"); ax[1].legend(fontsize=7, loc="upper left")
ax[1].set_title("(b) critical-mass correction from finite tension", fontsize=8.5)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig7_wall_extrapolation.png"), dpi=170); plt.close(fig)
print("figures written")
