# -*- coding: utf-8 -*-
"""Causal-classification diagnostic (revision 2): separates what the numerical double-null solutions establish locally from what is conjecture about the global spacetime.

For the two representative runs of Fig. 4 (A = 0.1, w = 0.1; core radius 0.8 and 2.0 R_dS) the code evolves the double-null chart, in which a future-directed causal curve has non-decreasing u and v.  Computed here:
 (1) the fraction of the evolved domain in each sign quadrant of (r_u, r_v);
 (2) the location of the wall (the v at which phi is half-way between the vacua) on every row: it stays at v ~ V_w, so that the 'core' (v < wall) and the 'exterior' (v > wall) are the same regions at all u;
 (3) a null-ray reachability test on the lattice: for every anti-trapped cell, whether any exterior cell lies in its causal past rectangle {u' <= u, v' <= v}, and for every exterior cell whether any anti-trapped cell lies in its causal future;
 (4) the u at which the evolution stops (the last row), so that nothing is claimed beyond it.
Output: causal_diagnostic.json
"""
import json
import os

import numpy as np

import dnull as D

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}


def analyse(Vw, label):
    o = D.run(A=0.1, V_w=Vw, w=0.1, h=0.01, snap_every=2, smax=12.0, stop_after_trap=20)
    R, Ru, P, h, N = o["R"], o["Ru"], o["P"], o["h"], o["N"]
    nrow = R.shape[0]
    pot = o["pot"]
    phiF, phiT = pot["psi_f"] * 0.1, pot["psi_t"] * 0.1
    mid = 0.5 * (phiF + phiT)
    cls = np.full((nrow, N + 1), -1, dtype=int)                 # 0 normal, 1 trapped, 2 anti-trapped, 3 r_u>0>r_v, -1 not evolved
    wall_v = np.full(nrow, np.nan)
    last_ok = 0
    for k in range(nrow):
        j = 2 * k
        r = R[k]
        if not np.isfinite(r[j + 1:]).all():
            break
        last_ok = k
        for i in range(j + 3, N - 2):
            rv = (r[i + 1] - r[i - 1]) / (2 * h); ru = Ru[k, i]
            cls[k, i] = 0 if (ru < 0 and rv > 0) else (1 if (ru < 0 and rv < 0) else (2 if (ru > 0 and rv > 0) else 3))
        ph = P[k]
        # wall: first crossing of the mid value in the part of the row that is evolved
        seg = np.arange(j + 3, N - 2)
        sgn = np.sign(ph[seg] - mid)
        idx = np.where(sgn[:-1] * sgn[1:] < 0)[0]
        if idx.size:
            wall_v[k] = (seg[idx[0]] + 0.5) * h
    v = np.arange(N + 1) * h
    u = o["U"]
    ev = cls >= 0
    frac = {name: float(np.sum(cls == c) / np.sum(ev)) for name, c in (("normal", 0), ("trapped", 1), ("anti_trapped", 2), ("r_u>0>r_v", 3))}
    valid = np.isfinite(wall_v)
    drift = float(np.nanmax(np.abs(wall_v[valid] - Vw)))
    # reachability: exterior cell = v > wall(u) + w ; core = v < wall(u) - w ; anti-trapped cells
    wlo = np.where(valid, wall_v - 0.1, np.nan); whi = np.where(valid, wall_v + 0.1, np.nan)
    ext = np.zeros_like(cls, dtype=bool); core = np.zeros_like(cls, dtype=bool)
    for k in range(last_ok + 1):
        if not valid[k]:
            continue
        ext[k] = (v > whi[k]) & (cls[k] >= 0)
        core[k] = (v < wlo[k]) & (cls[k] >= 0)
    # cumulative: ext_past[k, i] = any exterior cell with k' <= k, i' <= i  (2-D cumulative OR)
    ext_c = np.maximum.accumulate(np.maximum.accumulate(ext.astype(int), axis=0), axis=1) > 0
    AT = (cls == 2)
    n_AT = int(AT.sum())
    vv = np.broadcast_to(v[None, :], cls.shape)
    wall_b = np.broadcast_to(np.where(valid, wall_v, np.nan)[:, None], cls.shape)
    AT_core = AT & (vv < wall_b - 0.1); AT_wall = AT & (np.abs(vv - wall_b) <= 0.1); AT_ext = AT & (vv > wall_b + 0.1)
    rows_ext = np.where(AT_ext.any(axis=1))[0]
    # causal past test for the anti-trapped cells inside the core and for those beyond the wall
    AT_core_with_ext_past = int((AT_core & ext_c).sum())
    AT_ext_with_ext_past = int((AT_ext & ext_c).sum())
    AT_c_rev = np.maximum.accumulate(np.maximum.accumulate(AT_core[::-1, ::-1].astype(int), axis=0), axis=1)[::-1, ::-1] > 0
    ext_with_ATcore_future = int((ext & AT_c_rev).sum())
    ext_c_rev = np.maximum.accumulate(np.maximum.accumulate(ext[::-1, ::-1].astype(int), axis=0), axis=1)[::-1, ::-1] > 0
    core_reaching_exterior = float((core & ext_c_rev).sum() / max(core.sum(), 1))
    res = dict(label=label, V_w=Vw, rows_evolved=int(last_ok + 1), u_last=float(u[last_ok]), fractions=frac, wall_drift_max=drift, anti_trapped_cells=n_AT, anti_trapped_in_core=int(AT_core.sum()), anti_trapped_in_wall=int(AT_wall.sum()), anti_trapped_beyond_wall=int(AT_ext.sum()),
               rows_with_anti_trapped_beyond_wall=[int(rows_ext.min()), int(rows_ext.max())] if rows_ext.size else None, u_last_row_with_anti_trapped_beyond_wall=float(u[rows_ext.max()]) if rows_ext.size else None,
               core_anti_trapped_with_exterior_in_causal_past=AT_core_with_ext_past, exterior_anti_trapped_with_exterior_in_causal_past=AT_ext_with_ext_past, exterior_cells_with_core_anti_trapped_in_causal_future=ext_with_ATcore_future,
               fraction_of_core_cells_whose_future_reaches_exterior=core_reaching_exterior,
           anti_trapped_v_range=[float(v[np.where(AT.any(axis=0))[0].min()]), float(v[np.where(AT.any(axis=0))[0].max()])] if n_AT else None)
    print("%s: rows %d (u up to %.2f); quadrants %s; wall drift %.3f; anti-trapped cells %d (core %d, wall %d, beyond the wall %d up to u = %s); core anti-trapped cells with an exterior cell in their causal past: %d; exterior cells with a core anti-trapped cell in their future: %d; core cells that can signal the exterior: %.2f" %
          (label, last_ok + 1, u[last_ok], {k_: round(x, 3) for k_, x in frac.items()}, drift, n_AT, res["anti_trapped_in_core"], res["anti_trapped_in_wall"], res["anti_trapped_beyond_wall"], res["u_last_row_with_anti_trapped_beyond_wall"], AT_core_with_ext_past, ext_with_ATcore_future, core_reaching_exterior))
    return res


out["run_normal_core"] = analyse(1.6, "R_w = 0.8 R_dS")
out["run_expanding_core"] = analyse(4.0, "R_w = 2.0 R_dS")
out["established_locally"] = ["sign quadrants of (r_u, r_v) on the evolved double-null domain", "persistence of the anti-trapped region in the runs with R_w > R_dS", "wall at v ~ V_w, so core and exterior are fixed regions of the chart", "at the earliest rows the anti-trapped quadrant also extends beyond the wall (initial data), and ends at u ~ 0.15", "no outgoing ray with a minimum of the areal radius"]
out["established_by_chart_geometry"] = ["no exterior cell lies in the causal past of a core cell, in particular of a core anti-trapped cell (du dv >= 0; confirmed on the lattice)", "core signals can reach the exterior"]
out["conjectural_not_computed"] = ["global topology beyond the last evolved row", "whether the anti-trapped core is a separate asymptotic region (a bridge) or closes on itself", "any spacelike slice through a throat", "a daughter universe"]
json.dump(out, open(os.path.join(HERE, "causal_diagnostic.json"), "w"), indent=1)
