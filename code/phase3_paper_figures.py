# -*- coding: utf-8 -*-
"""Figures for the Phase 3 manuscript, drawn from the JSON outputs of the branch phase3-recycling.   python phase3_paper_figures.py"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(os.path.dirname(HERE), "figures")
L = lambda f: json.load(open(os.path.join(HERE, f), encoding="utf-8"))
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 150, "savefig.bbox": "tight"})
C1, C2, C3, C4 = "#1f4e79", "#c0392b", "#2e7d32", "#7b5ea7"

# ------------------------------------------------------------------ Fig 1: barotropic bounce criterion (unpolarized massless ECD gas)
E0 = L("phase3d0_ecd.json")["C_baseline"]["unpolarized massless (g=2)"]
Cc, Bb, nb = E0["C"], E0["B"], E0["n_bounce_analytic"]
x = np.linspace(1e-3, 1.25, 600)
n = x * nb
rho = Cc * n ** (4 / 3) - Bb * n ** 2
p = Cc / 3 * n ** (4 / 3) - Bb * n ** 2
mu = 4 / 3 * Cc * n ** (1 / 3) - 2 * Bb * n
dpdn = 4 / 9 * Cc * n ** (1 / 3) - 2 * Bb * n
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
sc = Cc * nb ** (4 / 3)
ax[0].plot(x, rho / sc, color=C1, label=r"$\rho$"); ax[0].plot(x, p / sc, color=C2, label=r"$p$"); ax[0].axhline(0, color="k", lw=0.6)
for xv, lab, c in ((E0["nsp_over_nb"], r"$dp/dn=0$", C4), (E0["nstar_over_nb"], r"$\mu=0$", C3), (1.0, "bounce", "k")):
    ax[0].axvline(xv, color=c, ls="--", lw=0.9); ax[0].text(xv + 0.01, 0.62, lab, rotation=90, fontsize=7, color=c, va="top")
ax[0].set_xlabel(r"$n/n_b$"); ax[0].set_ylabel(r"energy density, pressure  $/\,C n_b^{4/3}$"); ax[0].set_ylim(-0.7, 0.7); ax[0].legend(loc="lower left", fontsize=8); ax[0].set_title("(a) mean-field gas, $\\rho=Cn^{4/3}-Bn^2$", fontsize=9)
ax[1].plot(x, mu / (Cc * nb ** (1 / 3)), color=C3, label=r"$\mu=d\rho/dn$"); ax[1].plot(x, dpdn / (Cc * nb ** (1 / 3)), color=C4, label=r"$dp/dn$"); ax[1].axhline(0, color="k", lw=0.6)
ax[1].fill_between(x, -1.4, 1.4, where=dpdn < 0, color=C4, alpha=0.10, label="negative compressibility")
ax[1].set_xlabel(r"$n/n_b$"); ax[1].set_ylabel(r"in units of $Cn_b^{1/3}$"); ax[1].set_ylim(-0.9, 1.0); ax[1].legend(fontsize=7, loc="lower left"); ax[1].set_title("(b) chemical potential and compressibility", fontsize=9)
fig.tight_layout(w_pad=2.2); fig.savefig(os.path.join(FIG, "FigP3_1_barotropic.png")); plt.close(fig)

# ------------------------------------------------------------------ Fig 2: viscous relaxation of the anisotropy threshold
H = L("phase3h0_viscous_spin_bounce.json")["rows"]
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
for s0, c in ((0.05, C1), (0.2, C2)):
    rows = [r for r in H if abs(r["s0"] - s0) < 1e-12 and r["c2"] == 1.0 and r["K"] > 0 and r["K"] <= 1.0]
    base = [r for r in H if abs(r["s0"] - s0) < 1e-12 and r["K"] == 0][0]["B_c"]
    Ks = [0.0] + [r["K"] for r in rows]; Bs = [base] + [r["B_c"] for r in rows]
    ax[0].plot(Ks, np.array(Bs) / base, "o-", color=c, label=r"initial shear $s_0=%.2f$ ($B_c^{\rm inv}=%.4g$)" % (s0, base))
    r25 = [r for r in H if abs(r["s0"] - s0) < 1e-12 and r["c2"] == 0.25][0]
    ax[0].plot([r25["K"]], [r25["B_c"] / base], "s", color=c, mfc="white")
ax[0].set_yscale("log"); ax[0].set_xlabel(r"Knudsen number $K=\tau|\Theta|$"); ax[0].set_ylabel(r"$B_c\,/\,B_c^{\rm inviscid}$"); ax[0].legend(fontsize=7); ax[0].set_title("(a) critical spin strength, $c^2=1$ (squares: $c^2=0.25$)", fontsize=9)
rows = [r for r in H if abs(r["s0"] - 0.05) < 1e-12 and r["K"] > 0 and r["K"] <= 1.0]
ax[1].bar([("K=%g, c²=%g" % (r["K"], r["c2"])) for r in rows], [r["S_over_rhod_max"] for r in rows], color=C3)
ax[1].axhline(0.5, color="k", ls="--", lw=0.8); ax[1].text(-0.4, 0.52, "working bound $|\\pi|<0.5\\,\\rho_{\\rm dust}$", fontsize=7)
ax[1].set_ylabel(r"max $|\pi|/\rho_{\rm dust}$"); ax[1].set_ylim(0, 0.8); ax[1].tick_params(axis="x", labelsize=7, rotation=30); ax[1].set_title("(b) size of the dissipative stress", fontsize=9)
fig.tight_layout(w_pad=2.2); fig.savefig(os.path.join(FIG, "FigP3_2_viscous.png")); plt.close(fig)

# ------------------------------------------------------------------ Fig 3: benchmark geometry (black bounce): curvature, tides, transit time
T4 = L("phase3d4_propertime.json"); F2 = L("phase3d2_freefall.json")
fig, ax = plt.subplots(1, 3, figsize=(7.8, 2.9))
def K(r, a, m=1.0):
    R = np.sqrt(r ** 2 + a ** 2); f = 1 - 2 * m / R; Rp = r / R; Rpp = a ** 2 / R ** 3; fp = 2 * m * r / R ** 3; fpp = 2 * m * (1 / R ** 3 - 3 * r ** 2 / R ** 5)
    A_ = fpp / 2; B_ = fp * Rp / (2 * R); C_ = -(fp * Rp / (2 * R) + f * Rpp / R); D_ = (1 - f * Rp ** 2) / R ** 2
    return 4 * A_ ** 2 + 8 * B_ ** 2 + 8 * C_ ** 2 + 4 * D_ ** 2
for aa, c in ((1.0, C1), (3.0, C2)):
    r = np.linspace(0, 8 * aa, 400); ax[0].plot(r / aa, K(r, aa) / K(0.0, aa), color=c, label="$a/m=%g$" % aa)
ax[0].set_yscale("log"); ax[0].set_xlabel("$r/a$"); ax[0].set_ylabel("$K(r)/K(0)$"); ax[0].legend(fontsize=7); ax[0].set_title("(a) Kretschmann profile", fontsize=9)
sc = F2["C_tidal_scaling"]; am = np.array([float(k) for k in sc]); K11 = np.array([v["K11_abs_max"] for v in sc.values()]); Kth = np.array([v["Kth_abs_max"] for v in sc.values()])
ax[1].loglog(am, K11, "o-", color=C1, label="radial"); ax[1].loglog(am, Kth, "s-", color=C2, label="transverse"); ax[1].loglog(am, 1 / am ** 3, "k:", lw=0.8, label=r"$m/a^3$"); ax[1].loglog(am, 2 / am ** 3, "k--", lw=0.8, label=r"$2m/a^3$")
ax[1].set_xticks([0.25, 0.5, 1, 1.9]); ax[1].set_xticklabels(["0.25", "0.5", "1", "1.9"]); ax[1].minorticks_off(); ax[1].set_xlabel("$a/m$"); ax[1].set_ylabel(r"max tidal eigenvalue $[m^{-2}]$"); ax[1].legend(fontsize=7); ax[1].set_title("(b) tides at the throat", fontsize=9)
ams = [0.5, 1.0, 1.9, 3.0, 6.0]
for E_, c in ((1, C1), (5, C2)):
    ax[2].plot(ams, [T4["B_transit_maps"]["a/m=%.1f eps=0.1 E=%g L/Lc=0" % (a_, E_)]["dtau"] for a_ in ams], "o-", color=c, label="$E=%d$" % E_)
ax[2].set_yscale("log"); ax[2].set_xlabel("$a/m$"); ax[2].set_ylabel(r"transit proper time $[m]$"); ax[2].legend(fontsize=7); ax[2].set_title(r"(c) region $K\geq0.1\,K(0)$", fontsize=9)
fig.tight_layout(w_pad=2.2); fig.savefig(os.path.join(FIG, "FigP3_3_benchmark.png")); plt.close(fig)

# ------------------------------------------------------------------ Fig 4: Einstein-scalar-Gauss-Bonnet closed bounces; frame diagnostic
E1 = L("phase3e1_corrections.json"); S0 = L("phase3e0_string.json")
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
bn = [r for r in E1["A_scan"] if r["k"] == 1 and r["lam"] < 0 and r["outcome"] == "bounce"]
ax[0].scatter([r["eps_T_max"] for r in bn], [-r["cT2_min"] for r in bn], s=10, color=C2)
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].axvline(0.3, color="k", ls="--", lw=0.8); ax[0].text(0.34, 0.97, "controlled regime ("+"left of line)", fontsize=7, va="top", transform=ax[0].get_xaxis_transform())
ax[0].set_xlabel(r"tensor-sector strength $\epsilon_T$ at the bounce"); ax[0].set_ylabel(r"$-\min c_T^2$ along the path"); ax[0].set_title("(a) closed GB bounces ($\\lambda\\xi''<0$): all have $c_T^2<0$", fontsize=9)
D = S0["D_summary"]["-1"]
labs = ["string frame", "same data,\nno GB term", "Einstein frame"]; vals = [D["string_frame_bounces"], D["string_bounces_also_without_GB"], D["einstein_frame_bounces"]]
ax[1].bar(labs, vals, color=[C1, C4, C3]); ax[1].set_ylabel("bounces in 160 runs ($k=-1$)"); ax[1].set_title("(b) frame diagnostic, tree-level heterotic model", fontsize=9)
for i, v in enumerate(vals): ax[1].text(i, v + 2, str(v), ha="center", fontsize=8)
fig.tight_layout(w_pad=2.2); fig.savefig(os.path.join(FIG, "FigP3_4_gb.png")); plt.close(fig)

# ------------------------------------------------------------------ Fig 5: scale-free window for free quantum fields
G3 = L("phase3d3_goldilocks.json")
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
la = np.linspace(0, 40, 200)
kq = np.log10(G3["B_constants"]["k_qi_static"]); kinf = np.log10(G3["B_constants"]["k_qi_infinity"]); kty = np.log10(G3["B_constants"]["k_typical"])
ax[0].axhline(kq, color=C2, label=r"required: $N/a^2\geq5.59$ (quantum inequality)"); ax[0].axhline(kinf, color=C2, ls="--", lw=0.8, label=r"$\geq8.4$ ($a/m	o\infty$)"); ax[0].axhline(kty, color=C2, ls=":", label=r"$\gtrsim10^3$ (typical state)")
ax[0].axhline(0.0, color=C1, label=r"allowed: $N/a^2\leq1$ (species cutoff)"); ax[0].fill_between(la, 0.0, kq, color="0.5", alpha=0.25); ax[0].text(20, kq / 2, "empty window", ha="center", va="center", fontsize=8)
ax[0].set_xlim(0, 40); ax[0].set_ylim(-0.6, 3.4); ax[0].set_xlabel(r"$\log_{10}(a/\ell_P)$"); ax[0].set_ylabel(r"$\log_{10}(N/a^2)$"); ax[0].legend(fontsize=6.5, loc="upper left"); ax[0].set_title("(a) species count versus throat size", fontsize=9)
ecd = G3["B_ecd_vs_cutoff"]
labs = ["unpol.\n(spinodal)", "unpol.\n(bounce)", "Weyl\n(spinodal)", "Weyl\n(bounce)"]
vals = [ecd["unpolarized massless (g=2)"]["kF_spinodal"], ecd["unpolarized massless (g=2)"]["kF_bounce"], ecd["single-helicity Weyl (g=1)"]["kF_spinodal"], ecd["single-helicity Weyl (g=1)"]["kF_bounce"]]
ax[1].bar(labs, vals, color=[C4, C1, C4, C1]); ax[1].axhline(1.0, color="k", ls="--", lw=0.8); ax[1].text(-0.4, 1.05, r"one-species cutoff $M_P$", fontsize=7)
ax[1].set_ylabel(r"Fermi momentum $k_F/M_P$"); ax[1].set_title("(b) Einstein–Cartan–Dirac gas", fontsize=9); ax[1].tick_params(axis="x", labelsize=7)
fig.tight_layout(w_pad=2.2); fig.savefig(os.path.join(FIG, "FigP3_5_window.png")); plt.close(fig)
print("figures written to", FIG)
