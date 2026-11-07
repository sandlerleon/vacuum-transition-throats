# -*- coding: utf-8 -*-
"""Wall-thickness extrapolation analysis (revision 2).  Which conclusions of the numerical experiments are insensitive to the wall, and which cannot be extrapolated to the physical wall?

In the experiments both the thickness and the tension of the wall scale with the field amplitude A: thickness / R_dS = A sqrt(2/kappa), tension of order A sqrt(kappa/2); A -> 0 is the thin, light limit.  The three amplitudes
A = 0.05, 0.1, 0.2 span a factor 4 in thickness and tension, against 28 and 30 orders of magnitude to the physical wall.  So the question is answered by (i) quantities fixed by theorem for any thickness, (ii) the measured
sensitivity of the numerical outputs to A within that factor 4, compared with the resolution floor of the scheme, and (iii) what is left unconstrained.
Outputs: wall_extrapolation.json.
"""
import json
import os

import numpy as np

import thinwall as TWm

HERE = os.path.dirname(os.path.abspath(__file__))
TW = json.load(open(os.path.join(HERE, "thinwall.json"), encoding="utf-8"))
VAL = json.load(open(os.path.join(HERE, "validation.json"), encoding="utf-8"))
rows = json.load(open(os.path.join(HERE, "collapse_scan.json"), encoding="utf-8"))
cls = {(c["A"], c["w"], c["R_w"]): c for c in VAL["classes"]}
out = {}

# (i) quantities fixed independently of the thickness
sig_phys = TW["sigma_model_hat"]; thick_phys = TW["thickness_over_RdS"]
out["physical_wall"] = dict(sigma_hat=sig_phys, thickness_over_RdS=thick_phys, gravity_parameter=TW["gravity_parameter_model"])
out["theorems"] = [
    "Proposition 1 (no minimum of the areal radius along a null ray): uses only the null energy condition of the scalar and the dust; independent of wall thickness and tension.",
    "Anti-trapped points in the initial core exist iff the core radius exceeds R_dS (de Sitter interior): a statement about the initial data, independent of the wall except for an offset of the order of its width.",
    "m_c(sigma) = (1+q) R_dS - sigma^2 [1/4 - 1/(6(1+q))] + O(sigma^4): the correction is 1e-60 at the physical tension, and a few per cent at the simulated tensions.",
]

# (ii) measured sensitivity
byA = {}
for A in (0.05, 0.1, 0.2):
    sel = [r for r in rows if r["A"] == A and r["first_trapped"] is not None and not cls[(A, r["w"], r["R_w"])]["AT"]]
    dev = np.array([r["first_trapped"]["r"] / r["r_h_SAdS"] - 1.0 for r in sel])
    floor = np.array([r["h"] / r["r_h_SAdS"] for r in sel])
    sg = sel[0]["sigma_est"]; th = sel[0]["thickness_over_RdS"]
    mc_sim = TWm.m_crit(sg, 1.0); mc_lim = 2.0
    ATsel = [c for c in VAL["classes"] if c["A"] == A and c["AT"]]
    byA[str(A)] = dict(n_trapped_normal_core=len(sel), thickness_over_RdS=th, sigma_hat_sim=sg, median_abs_dev=float(np.median(np.abs(dev))), max_abs_dev=float(np.max(np.abs(dev))), mean_signed_dev=float(np.mean(dev)),
                       resolution_floor_median=float(np.median(floor)), resolution_floor_max=float(np.max(floor)),
                       m_c_at_sim_tension=mc_sim, m_c_shift=mc_sim / mc_lim - 1.0,
                       min_m_over_mc_expanding_core_thin_limit=float(min(c["m_over_mc"] for c in ATsel)) if ATsel else None,
                       min_m_over_mc_expanding_core_finite_tension=float(min(c["m_over_mc"] * mc_lim / mc_sim for c in ATsel)) if ATsel else None)
    print("A = %.2f: thickness %.4f R_dS, tension %.3f, |dev| median %.4f max %.4f (floor h/r_h median %.4f, max %.4f); m_c(sigma_sim) shift %.4f ; min m/m_c over expanding cores: %.3f (thin limit) %.3f (finite tension)" %
          (A, th, sg, byA[str(A)]["median_abs_dev"], byA[str(A)]["max_abs_dev"], byA[str(A)]["resolution_floor_median"], byA[str(A)]["resolution_floor_max"], byA[str(A)]["m_c_shift"],
           byA[str(A)]["min_m_over_mc_expanding_core_thin_limit"], byA[str(A)]["min_m_over_mc_expanding_core_finite_tension"]))
out["by_amplitude"] = byA
# trend of the horizon-radius deviation with A: slope with bootstrap interval (three points, medians of runs) -- tests whether any thickness dependence is visible above the noise
rng = np.random.default_rng(4)
Avals = np.array([0.05, 0.1, 0.2]); slopes = []
dev_runs = {A: np.array([r["first_trapped"]["r"] / r["r_h_SAdS"] - 1.0 for r in rows if r["A"] == A and r["first_trapped"] is not None and not cls[(A, r["w"], r["R_w"])]["AT"]]) for A in (0.05, 0.1, 0.2)}
for _ in range(4000):
    med = [np.median(rng.choice(dev_runs[A], size=len(dev_runs[A]), replace=True)) for A in (0.05, 0.1, 0.2)]
    slopes.append(np.polyfit(Avals, med, 1)[0])
lo, hi = np.percentile(slopes, [2.5, 97.5])
out["deviation_slope_vs_A"] = dict(slope_median=float(np.median(slopes)), ci95=[float(lo), float(hi)], consistent_with_zero=bool(lo <= 0 <= hi))
print("deviation of the first trapped radius from r_h(SAdS) against A: slope %.3f, 95%% interval [%.3f, %.3f]: %s" % (np.median(slopes), lo, hi, "consistent with no dependence on A" if lo <= 0 <= hi else "dependence visible"))

# the wall width w (a separate parameter of the initial data, 0.05-0.2 R_dS) against the same deviation
byw = {}
allsel = [r for r in rows if r["first_trapped"] is not None and not cls[(r["A"], r["w"], r["R_w"])]["AT"]]
for w in (0.05, 0.1, 0.2):
    d = np.array([r["first_trapped"]["r"] / r["r_h_SAdS"] - 1.0 for r in allsel if r["w"] == w])
    byw[str(w)] = dict(n=int(len(d)), median_signed=float(np.median(d)), median_abs=float(np.median(np.abs(d))), max_abs=float(np.max(np.abs(d))))
out["by_width"] = byw
wv = np.array([r["w"] for r in allsel]); dv = np.array([r["first_trapped"]["r"] / r["r_h_SAdS"] - 1.0 for r in allsel])
sl = []; ic = []
for _ in range(4000):
    ii = rng.integers(0, len(wv), len(wv)); p_ = np.polyfit(wv[ii], dv[ii], 1); sl.append(p_[0]); ic.append(p_[1])
lo_w, hi_w = np.percentile(sl, [2.5, 97.5]); lo_i, hi_i = np.percentile(ic, [2.5, 97.5])
floor_med = float(np.median([r["h"] / r["r_h_SAdS"] for r in allsel]))
out["deviation_slope_vs_w"] = dict(slope_median=float(np.median(sl)), ci95=[float(lo_w), float(hi_w)], visible=bool(not (lo_w <= 0 <= hi_w)), intercept_at_w0_median=float(np.median(ic)), intercept_ci95=[float(lo_i), float(hi_i)], resolution_floor_median=floor_med,
                                   intercept_consistent_with_zero=bool(lo_i <= 0 <= hi_i))
out["deviation_points"] = [dict(w=float(w_), dev=float(d_)) for w_, d_ in zip(wv, dv)]
print("linear extrapolation of the deviation to w -> 0: %.4f, 95%% interval [%.4f, %.4f] (resolution floor %.4f)" % (np.median(ic), lo_i, hi_i, floor_med))
print("by wall width (median signed dev): %s ; slope of the signed deviation against w: %.3f, 95%% interval [%.3f, %.3f]: %s" % ({k_: round(v["median_signed"], 4) for k_, v in byw.items()}, np.median(sl), lo_w, hi_w, "a width dependence is visible (the widest wall places the first trapped point a few per cent inside the horizon radius)" if not (lo_w <= 0 <= hi_w) else "no width dependence visible"))
# expanding-core threshold in the core radius against the wall width (bracket)
th = []
for t in VAL["thresholds"]:
    th.append(dict(A=t["A"], w=t["w"], bracket=[t["R_w_no_AT"], t["R_w_AT"]], bracket_width=t["R_w_AT"] - t["R_w_no_AT"]))
out["core_threshold_brackets"] = th
print("expanding-core threshold brackets in R_w/R_dS by wall width:", {(t["A"], t["w"]): t["bracket"] for t in th})
# (iii) what cannot be extrapolated
out["extrapolation_factors"] = dict(thickness_ratio_sim_to_phys=float(0.0144337567 / thick_phys), tension_ratio_sim_to_phys=float(0.2309401076758503 / sig_phys),
                                    grid_points_needed_to_resolve_physical_wall=float(1.0 / thick_phys), grid_points_used=float(0.0144337567 / 0.00625))
print("simulated / physical thickness: %.2g ; tension: %.2g ; grid points per physical thickness needed over one R_dS: %.2g" % (out["extrapolation_factors"]["thickness_ratio_sim_to_phys"], out["extrapolation_factors"]["tension_ratio_sim_to_phys"], out["extrapolation_factors"]["grid_points_needed_to_resolve_physical_wall"]))
out["classification"] = dict(
    insensitive_by_theorem=["no minimum of the areal radius along null rays (energy condition)", "thin-wall critical mass m_c -> (1+q) R_dS as sigma -> 0, correction ~ sigma^2", "existence of anti-trapped points in the initial core iff its radius exceeds R_dS (offset ~ wall width)"],
    insensitive_within_resolution=["formation of a trapped surface for m_ext < m_c (all runs)", "radius of the first trapped point against the horizon radius of the measured mass (|dev| median below the resolution floor; no trend with A at 95%)", "classification of every expanding-core run as m_ext > m_c at either the limit or the finite simulated tension"],
    not_extrapolable=["the dynamics of a wall 1e-30 R_dS thick (neither resolved nor tested)", "any quantity at the 1e-3 level or finer in the horizon radius (below the resolution floor)", "wall stability and energy exchange with the exterior", "the fate of the matter after the code stops at the singularity", "non-spherical or quantum wall dynamics"])
json.dump(out, open(os.path.join(HERE, "wall_extrapolation.json"), "w"), indent=1)
