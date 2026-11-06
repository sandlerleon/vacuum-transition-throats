# -*- coding: utf-8 -*-
"""Additional validation for the revision: constraint-residual convergence of dynamical scalar-wall solutions, convergence of the null expansions, a detailed
classification of the runs that ended without a trapped surface, the creation of anti-trapped regions in runs whose initial core is normal, and the comparison of
other exterior vacua (q = -0.5, 0, 2). Writes validation_extra.json.

    python validate_extra.py
"""
import json
import os

import numpy as np

import collapse_scan as CS
import dnull as D

HERE = os.path.dirname(os.path.abspath(__file__))


def constraint_rms(A, Vw, w, h, u_lo=0.05, u_hi=0.28, margin=0.4, v_hi=3.0):
    """RMS (absolute, units R_dS = 1) of the two constraints C_vv = r_vv - s_v r_v + (r/2) p_v^2 and C_uu = r_uu - s_u r_u + (r/2) p_u^2 over (u, v) in [u_lo, u_hi] x [u + margin, v_hi]
    (away from the centre, the outer boundary and the horizon)."""
    o = D.run(A=A, V_w=Vw, w=w, h=h, snap_every=1, smax=12.0, stop_after_trap=0)
    R, P, S = o["R"], o["P"], o["S"]
    cvv, cuu, svv, suu = [], [], [], []
    kmax = int(round(u_hi / h))
    for k in range(max(1, int(round(u_lo / h))), min(kmax, R.shape[0] - 2)):
        for i in range(int(round((k * h + margin) / h)), int(round(v_hi / h))):
            r, p, s = R[k], P[k], S[k]
            rv = (r[i + 1] - r[i - 1]) / (2 * h); rvv = (r[i + 1] - 2 * r[i] + r[i - 1]) / h ** 2
            sv = (s[i + 1] - s[i - 1]) / (2 * h); pv = (p[i + 1] - p[i - 1]) / (2 * h)
            cvv.append(rvv - sv * rv + 0.5 * r[i] * pv ** 2)
            svv.append(max(abs(rvv), abs(sv * rv), abs(0.5 * r[i] * pv ** 2)))
            ru = (R[k + 1][i] - R[k - 1][i]) / (2 * h); ruu = (R[k + 1][i] - 2 * R[k][i] + R[k - 1][i]) / h ** 2
            su = (S[k + 1][i] - S[k - 1][i]) / (2 * h); pu = (P[k + 1][i] - P[k - 1][i]) / (2 * h)
            cuu.append(ruu - su * ru + 0.5 * r[i] * pu ** 2)
            suu.append(max(abs(ruu), abs(su * ru), abs(0.5 * r[i] * pu ** 2)))
    cvv, cuu, svv, suu = np.array(cvv), np.array(cuu), np.array(svv), np.array(suu)
    return dict(h=h, n=int(len(cvv)), rms_vv=float(np.sqrt(np.mean(cvv ** 2))), rms_uu=float(np.sqrt(np.mean(cuu ** 2))), max_vv=float(np.abs(cvv).max()), max_uu=float(np.abs(cuu).max()),
                scale_vv=float(np.sqrt(np.mean(svv ** 2))), scale_uu=float(np.sqrt(np.mean(suu ** 2))))


def expansions_at(A, Vw, w, h, u0=0.20, v_lo=0.8, v_hi=3.0, dv=0.1):
    """r_u and r_v at u = u0 on a fixed set of v values (the lattices of the three resolutions share these points)."""
    o = D.run(A=A, V_w=Vw, w=w, h=h, snap_every=1, smax=12.0, stop_after_trap=0)
    R, Ru = o["R"], o["Ru"]
    k = int(round(u0 / h))
    vs = np.arange(v_lo, v_hi + 1e-9, dv)
    idx = [int(round(v / h)) for v in vs]
    r = R[k]
    rv = np.array([(r[i + 1] - r[i - 1]) / (2 * h) for i in idx])
    # r_u by a central difference in u (rows k +- 1)
    ru = np.array([(R[k + 1][i] - R[k - 1][i]) / (2 * h) for i in idx])
    return vs, np.array([r[i] for i in idx]), ru, rv


def expansion_convergence(A, Vw, w, hs=(0.02, 0.01, 0.005), u0=0.20):
    res = [expansions_at(A, Vw, w, h, u0=u0) for h in hs]
    out = dict(h=list(hs))
    for name, j in (("r", 1), ("r_u", 2), ("r_v", 3)):
        d1 = float(np.abs(res[0][j] - res[1][j]).max()); d2 = float(np.abs(res[1][j] - res[2][j]).max())
        out[name] = dict(diff_h1_h2=d1, diff_h2_h3=d2, ratio=float(d1 / d2) if d2 > 0 else None)
    out["sign_agreement"] = bool(np.all(np.sign(res[1][2]) == np.sign(res[2][2])) and np.all(np.sign(res[1][3]) == np.sign(res[2][3])))
    return out


def classify_inconclusive(S):
    rows = []
    for r in S:
        if r["first_trapped"] is not None:
            continue
        rec = dict(A=r["A"], w=r["w"], R_w=r["R_w"], m_ext=r["m_ext"], r_h=r["r_h_SAdS"], ended_by=r["ended_by"], u_end=r["u_end"], min_rv=r["min_rv_last_row"], max_sigma=r["max_sigma_last_row"])
        for fac in (2, 4):
            h = min(0.01, r["w"] / 8) / fac
            a = CS.analyse(r["A"], r["V_w"], r["w"], h)
            rec["refined_h_%g" % h] = dict(first_trapped=a["first_trapped"], ended_by=a["ended_by"], u_end=a["u_end"], max_sigma=a["max_sigma_last_row"])
        rows.append(rec)
    return rows


def main():
    S = json.load(open(os.path.join(HERE, "collapse_scan.json")))
    SQ = json.load(open(os.path.join(HERE, "collapse_scan_q.json")))
    out = {}
    out["constraints"] = [constraint_rms(0.1, 2.0, 0.1, h) for h in (0.02, 0.01, 0.005)]
    cs = out["constraints"]
    out["constraint_orders"] = dict(vv=[float(np.log2(cs[i]["rms_vv"] / cs[i + 1]["rms_vv"])) for i in range(2)], uu=[float(np.log2(cs[i]["rms_uu"] / cs[i + 1]["rms_uu"])) for i in range(2)])
    out["constraints_AT"] = [constraint_rms(0.1, 4.0, 0.1, h, u_hi=0.15, v_hi=4.0) for h in (0.02, 0.01, 0.005)]
    ca = out["constraints_AT"]
    out["constraint_orders_AT"] = dict(vv=[float(np.log2(ca[i]["rms_vv"] / ca[i + 1]["rms_vv"])) for i in range(2)], uu=[float(np.log2(ca[i]["rms_uu"] / ca[i + 1]["rms_uu"])) for i in range(2)])
    out["expansions"] = dict(normal_core=expansion_convergence(0.1, 2.0, 0.1), expanding_core=expansion_convergence(0.1, 4.0, 0.1, u0=0.12))
    out["inconclusive"] = classify_inconclusive(S)
    init_norm = [r for r in S if r["at_points_initial_core"] == 0]
    out["creation"] = dict(n_runs=len(S), n_initially_normal_core=len(init_norm), n_created_later=int(sum(1 for r in init_norm if r["at_points_later_rows"] > 0)),
                           n_initially_expanding_core=len(S) - len(init_norm))
    # other exterior vacua
    byq = {}
    for q in sorted({r["q"] for r in SQ}):
        sub = [r for r in SQ if r["q"] == q]
        dev = [r["first_trapped"]["r"] / r["r_h_SAdS"] - 1 for r in sub if r["first_trapped"] is not None and r["r_h_SAdS"]]
        byq[str(q)] = dict(n=len(sub), n_trapped=len(dev), dev_median=float(np.median(np.abs(dev))) if dev else None, dev_max=float(np.max(np.abs(dev))) if dev else None,
                           n_expanding_core=int(sum(1 for r in sub if r["at_points_last_row"] > 0)), n_expanding_core_subcritical=int(sum(1 for r in sub if r["at_points_last_row"] > 0 and r["m_ext"] < (1 + q) * 1.0)),
                           n_no_bh_horizon=int(sum(1 for r in sub if r["r_h_SAdS"] is None)))
    out["other_vacua"] = byq
    json.dump(out, open(os.path.join(HERE, "validation_extra.json"), "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("constraints", "constraint_orders", "constraint_orders_AT", "expansions", "creation", "other_vacua")}, indent=1)[:6000])
    print(json.dumps(out["inconclusive"], indent=1)[:4000])


if __name__ == "__main__":
    main()
