# -*- coding: utf-8 -*-
"""Figures of the manuscript from validation.json, collapse_scan.json and thinwall.json (and two fresh runs for the sign-map panels).

    python make_figures.py          (writes ../figures/Fig1_...Fig5_*.png)
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

import dnull as D
import thinwall as TW

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 150})
S = json.load(open(os.path.join(HERE, "collapse_scan.json")))
V = json.load(open(os.path.join(HERE, "validation.json")))
T = json.load(open(os.path.join(HERE, "thinwall.json")))
COL = {0.05: "#1f77b4", 0.1: "#d62728", 0.2: "#2ca02c"}


def fig1():
    """Apparent-horizon radius against the Schwarzschild-AdS horizon of the measured mass; convergence with h."""
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.8))
    for A in sorted(COL):
        rs = [(r["r_h_SAdS"], r["first_trapped"]["r"]) for r in S if r["A"] == A and r["first_trapped"] is not None]
        ax[0].plot([a for a, b in rs], [b for a, b in rs], "o", ms=4, color=COL[A], label="A = %.2f" % A)
    xs = np.linspace(0.1, 4.0, 10)
    ax[0].plot(xs, xs, "k-", lw=0.8)
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel("horizon radius of the measured mass, $r_h$ ($r_h^3+r_h=m_{ext}$)"); ax[0].set_ylabel("radius of the first trapped point, $r$")
    ax[0].legend(); ax[0].set_title("(a) all %d runs with a trapped surface" % V["scan"]["n_trapped"])
    hs = [t["h"] for t in V["trapped_convergence"]]
    rr = [t["r"] for t in V["trapped_convergence"]]
    ax[1].plot(hs, rr, "o-", color="#d62728")
    ax[1].set_xlabel("grid step $h$"); ax[1].set_ylabel("radius of the first trapped point")
    ax[1].set_title("(b) convergence ($R_w$ = 1, w = 0.1, A = 0.1)")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig2_horizon_check.png")); plt.close(fig)


def fig2():
    """Outcome map in core radius and exterior mass."""
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    for c in V["classes"]:
        if not c["T"]:
            ax.plot(c["R_w"], c["m_over_mc"], "+", color="0.4", ms=7)
        elif c["AT"]:
            ax.plot(c["R_w"], c["m_over_mc"], "s", color=COL[c["A"]], ms=5, mfc=COL[c["A"]])
        else:
            ax.plot(c["R_w"], c["m_over_mc"], "o", color=COL[c["A"]], ms=5, mfc="none")
    ax.axvline(1.0, color="k", ls="--", lw=0.8); ax.axhline(1.0, color="k", ls=":", lw=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("core radius at the wall, $R_w / R_{dS}$"); ax.set_ylabel("exterior mass, $m_{ext} / m_c$")
    ax.plot([], [], "+", color="0.4", label="inconclusive (run ended first)"); ax.plot([], [], "o", mfc="none", color="gray", label="trapped surface only")
    ax.plot([], [], "s", color="gray", label="trapped surface and expanding (anti-trapped) core")
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig3_outcome_map.png")); plt.close(fig)


def quadrant_map(Vw, ax, title):
    o = D.run(A=0.1, V_w=Vw, w=0.1, h=0.01, snap_every=2, smax=12.0, stop_after_trap=20)
    R, Ru, h = o["R"], o["Ru"], o["h"]
    N = o["N"]
    nrow = R.shape[0]
    img = np.full((nrow, N + 1), np.nan)
    for k in range(nrow):
        j = 2 * k
        r = R[k]
        if not np.isfinite(r[j + 1:]).all():
            break
        for i in range(j + 3, N - 2):
            rv = (r[i + 1] - r[i - 1]) / (2 * h)
            ru = Ru[k, i]
            img[k, i] = 0 if (ru < 0 and rv > 0) else (1 if (ru < 0 and rv < 0) else (2 if (ru > 0 and rv > 0) else 3))
    cmap = ListedColormap(["#dddddd", "#d62728", "#1f77b4", "#ff7f0e"])
    u = o["U"]; v = np.arange(N + 1) * h
    ax.pcolormesh(v, u[:nrow], img, cmap=cmap, vmin=-0.5, vmax=3.5, shading="auto")
    ax.axvline(Vw, color="k", lw=0.8, ls="--")
    ax.set_xlabel("v"); ax.set_ylabel("u"); ax.set_title(title, fontsize=9)
    ax.set_xlim(0, v[-1]); ax.set_ylim(0, u[nrow - 1])


def fig3():
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.9))
    quadrant_map(1.6, ax[0], "(a) core radius 0.8 $R_{dS}$ (trapped surface only)")
    quadrant_map(4.0, ax[1], "(b) core radius 2.0 $R_{dS}$ (trapped surface and expanding core)")
    from matplotlib.patches import Patch
    ax[1].legend(handles=[Patch(color="#dddddd", label="normal"), Patch(color="#d62728", label="trapped"), Patch(color="#1f77b4", label="anti-trapped (expanding)"),
                          Patch(color="#ff7f0e", label="$r_u>0, r_v<0$")], fontsize=7, loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig4_sign_maps.png")); plt.close(fig)


def fig4():
    """Thin-wall allowed region in the (R, m) plane."""
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.9))
    for a, sigma in zip(ax, (0.9, 0.3)):
        R = np.linspace(0.0, 1.0, 4001)
        mp = TW.m_pm(R, sigma, 1.0, 3.0, +1); mm = TW.m_pm(R, sigma, 1.0, 3.0, -1)
        a.fill_between(R, mm, mp, color="#f4b6b6", label="forbidden")
        a.plot(R, mp, "k-", lw=0.8); a.plot(R, mm, "k-", lw=0.8)
        mc = TW.m_crit(sigma, 1.0)
        a.axhline(mc, color="#1f77b4", ls="--", lw=1.0, label="$m_c$ = %.3f" % mc)
        a.axvline(1.0, color="gray", lw=0.6)
        a.set_xlabel("wall radius $R / R_{dS}$"); a.set_ylabel("exterior mass $m / R_{dS}$"); a.set_ylim(0, 2.3); a.set_xlim(0, 1.2)
        a.set_title("tension $\\sigma R_{dS}/M_{Pl}^2$ = %.1f, q = 1" % sigma, fontsize=9); a.legend(fontsize=7, loc="upper left")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig1_thin_wall.png")); plt.close(fig)


def fig5():
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.8))
    q = [r["q"] for r in T["mass_table"]]
    ax[0].plot(q, [r["M_c_over_matter_in_Hubble_sphere"] for r in T["mass_table"]], "o-", color="#d62728", label="vs matter within one Hubble radius")
    ax[0].plot(q, [r["deficit_over_matter_in_RdS_sphere"] for r in T["mass_table"]], "s--", color="#1f77b4", label="vs matter within a sphere of radius $R_{dS}$")
    ax[0].set_xlabel("depth of the true vacuum, q"); ax[0].set_ylabel("critical mass $M_c$ / matter mass (today)")
    ax[0].set_title("(a) vacuum-energy deficit against matter"); ax[0].legend(fontsize=7)
    a = [r["a"] for r in T["matter_over_deficit_vs_a"]]; d = [r["ratio"] for r in T["matter_over_deficit_vs_a"]]
    ax[1].plot(a, d, "-", color="#1f77b4")
    ax[1].axhline(1.0, color="k", lw=0.6, ls=":")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("scale factor $a$ of the parent (matter density falls as $a^{-3}$)"); ax[1].set_ylabel("matter density / $(1+q)V_f$")
    ax[1].set_title("(b) q = 1: matter in a sphere of radius $R_{dS}$ against $M_c$")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "Fig5_assembly.png")); plt.close(fig)


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5()
    print("figures written to", FIG)
