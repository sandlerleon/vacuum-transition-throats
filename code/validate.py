# -*- coding: utf-8 -*-
"""Validation and summary of the numerical experiment: exact-solution tests, convergence of the first trapped point, and statistics of collapse_scan.json.
Writes validation.json.

    python validate.py
"""
import json
import os

import numpy as np

import collapse_scan as CS
import dnull as D
import test_dnull as TD
import thinwall as TW

HERE = os.path.dirname(os.path.abspath(__file__))


def first_trapped(A, Vw, w, h):
    o = D.run(A=A, V_w=Vw, w=w, h=h, snap_every=1, smax=12.0)
    R, Ru, hh = o["R"], o["Ru"], o["h"]
    for k in range(R.shape[0]):
        r = R[k]
        if not np.isfinite(r[k + 1:]).all():
            break
        rv = np.zeros_like(r); rv[1:-1] = (r[2:] - r[:-2]) / (2 * hh)
        idx = [i for i in range(k + 4, len(r) - 2) if rv[i] < 0 and Ru[k, i] < 0]
        if idx:
            i = idx[0]
            return dict(h=h, u=k * hh, v=i * hh, r=float(r[i]))
    return None


def null_ray_minima(S):
    """Numerical counterpart of Proposition 1: count, over every row of every run, the places where r_v changes sign from negative to positive going outward
    (a minimum of the areal radius along an outgoing null ray; thresholds 1e-3 on |r_v| exclude noise)."""
    cnt, rows = 0, 0
    for r0 in S:
        h = min(0.01, r0["w"] / 8)
        o = D.run(A=r0["A"], V_w=r0["V_w"], w=r0["w"], h=h, snap_every=1, smax=12.0, stop_after_trap=40)
        R, N = o["R"], o["N"]
        for k in range(R.shape[0]):
            r = R[k]
            if not np.isfinite(r[k + 1:]).all():
                break
            rows += 1
            rv = (r[2:] - r[:-2]) / (2 * h)
            for i in range(k + 6, N - 10):
                if rv[i - 1] < -1e-3 and rv[i] > 1e-3:
                    cnt += 1
    return dict(runs=len(S), rows=rows, minima=cnt)


def main():
    out = {}
    # ---- exact solutions: Misner-Sharp mass of de Sitter and anti-de Sitter with a regular centre
    tests = {}
    for vac, name in ((+1, "dS"), (-1, "AdS")):
        hs = (0.04, 0.02, 0.01)
        e = [TD.ms_error(vac, h, urows=0.6) for h in hs]
        tests[name] = dict(h=list(hs), max_rel_error=[float(x) for x in e], orders=[float(np.log2(e[i] / e[i + 1])) for i in range(2)])
    out["vacuum_tests"] = tests
    # ---- convergence of the first trapped point (core radius 1, wall width 0.1, A = 0.1)
    out["trapped_convergence"] = [first_trapped(0.1, 2.0, 0.1, h) for h in (0.02, 0.01, 0.005)]
    # ---- statistics of the scan
    S = json.load(open(os.path.join(HERE, "collapse_scan.json")))
    dev, rows = [], []
    for r in S:
        ft = r["first_trapped"]
        if ft is not None and r["r_h_SAdS"]:
            dev.append(ft["r"] / r["r_h_SAdS"] - 1.0)
        rows.append(r)
    dev = np.array(dev)
    out["null_ray_minima"] = null_ray_minima(S)
    out["scan"] = dict(n_runs=len(S), n_trapped=int(sum(1 for r in S if r["first_trapped"] is not None)), n_no_trapped=int(sum(1 for r in S if r["first_trapped"] is None)),
                       dev_median=float(np.median(dev)), dev_max_abs=float(np.max(np.abs(dev))), dev_p95=float(np.percentile(np.abs(dev), 95)))
    per_A = {}
    for A in sorted({r["A"] for r in S}):
        d = [r["first_trapped"]["r"] / r["r_h_SAdS"] - 1.0 for r in S if r["A"] == A and r["first_trapped"] is not None]
        per_A[str(A)] = dict(median_abs=float(np.median(np.abs(d))), max_abs=float(np.max(np.abs(d))), sigma_est=float([r["sigma_est"] for r in S if r["A"] == A][0]),
                             thickness=float([r["thickness_over_RdS"] for r in S if r["A"] == A][0]))
    out["per_A"] = per_A
    # outcome table: T only vs T + anti-trapped core, against the core radius and m_ext / m_c (m_c = 2 R_dS for q = 1)
    mc = 2.0
    cls = []
    for r in S:
        T = r["first_trapped"] is not None
        AT = r["at_points_last_row"] > 0
        cls.append(dict(A=r["A"], w=r["w"], R_w=r["R_w"], m_ext=r["m_ext"], m_over_mc=r["m_ext"] / mc, R_M=r["R_M"], T=T, AT=AT))
    out["classes"] = cls
    # thresholds: largest R_w without and smallest R_w with an anti-trapped core, per (A, w)
    th = []
    for A in sorted({c["A"] for c in cls}):
        for w in sorted({c["w"] for c in cls}):
            sub = sorted([c for c in cls if c["A"] == A and c["w"] == w], key=lambda c: c["R_w"])
            no = [c for c in sub if not c["AT"]]
            yes = [c for c in sub if c["AT"]]
            th.append(dict(A=A, w=w, R_w_no_AT=no[-1]["R_w"] if no else None, m_no_AT=no[-1]["m_over_mc"] if no else None,
                           R_w_AT=yes[0]["R_w"] if yes else None, m_AT=yes[0]["m_over_mc"] if yes else None))
    out["thresholds"] = th
    json.dump(out, open(os.path.join(HERE, "validation.json"), "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("vacuum_tests", "trapped_convergence", "scan", "per_A")}, indent=1))


if __name__ == "__main__":
    main()
