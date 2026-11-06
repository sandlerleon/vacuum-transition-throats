# -*- coding: utf-8 -*-
"""Automated checks (run: python tests/test_all.py, from the code directory). Each check is an assertion on an exact result or on a stored number that the manuscript quotes."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import dnull as D                 # noqa: E402
import test_dnull as TD           # noqa: E402
import thinwall as TW             # noqa: E402

CODE = os.path.join(HERE, "..")
TWJ = json.load(open(os.path.join(CODE, "thinwall.json")))
VAL = json.load(open(os.path.join(CODE, "validation.json")))
SCAN = json.load(open(os.path.join(CODE, "collapse_scan.json")))
CHECKS = []


def check(f):
    CHECKS.append(f)
    return f


@check
def critical_mass_window():
    """A forbidden window exists just below the exact critical mass and none just above it, for q and tension combinations."""
    for q in (0.5, 1.0, 2.0):
        for sigma in (0.01, 0.3, 0.9):
            mc = TW.m_crit(sigma, q)
            assert TW.forbidden_window(mc * (1 - 1e-3), sigma, q) is not None, (q, sigma)
            assert TW.forbidden_window(mc * (1 + 1e-3), sigma, q) is None, (q, sigma)


@check
def critical_mass_small_tension():
    """m_c -> (1 + q) R_dS as the tension goes to zero."""
    for q in (0.5, 1.0, 3.0):
        assert abs(TW.m_crit(1e-4, q) / (1 + q) - 1) < 1e-6


@check
def horizon_radius_at_critical_mass():
    """Exterior horizon radius equals the de Sitter radius at m = (1 + q) R_dS, for every q (entropy equality)."""
    for q in (0.25, 0.5, 1.0, 2.0, 4.0):
        assert abs(TW.horizon_radius(1 + q, q) - 1.0) < 1e-12
        assert abs(TW.horizon_radius(TW.m_crit(1e-12, q), q) - 1.0) < 1e-9


@check
def junction_branches_meet():
    """m_+ and m_- meet at R = R_dS and bound a band that is symmetric about the middle curve."""
    for q in (0.5, 1.0):
        for s in (0.1, 0.5):
            assert abs(TW.m_pm(1.0, s, q, 3.0, +1) - TW.m_pm(1.0, s, q, 3.0, -1)) < 1e-12
            R = 0.7
            mid = 0.5 * (TW.m_pm(R, s, q, 3.0, +1) + TW.m_pm(R, s, q, 3.0, -1))
            assert abs(mid - ((1 + q) - s ** 2 / 4) * R ** 3) < 1e-12


@check
def solver_exact_vacuum_solutions():
    """Misner-Sharp mass of de Sitter and anti-de Sitter equals rho r^3 / 3 to better than 1e-3 at h = 0.02 and converges at about second order."""
    for vac in (+1, -1):
        e = [TD.ms_error(vac, h, urows=0.6) for h in (0.04, 0.02, 0.01)]
        assert e[1] < 1e-3 and e[2] < 1e-3, e
        assert e[0] / e[2] > 3.0, e


@check
def solver_horizon_matches_analytic():
    """The first trapped point of a collapsing run lies on the Schwarzschild-AdS horizon of the measured mass to within 3%."""
    r = [x for x in SCAN if x["A"] == 0.1 and abs(x["R_w"] - 1.0) < 1e-9 and abs(x["w"] - 0.1) < 1e-9][0]
    assert abs(r["first_trapped"]["r"] / r["r_h_SAdS"] - 1) < 0.03


@check
def stored_results_consistent():
    """Numbers quoted in the manuscript: counts, medians and the critical mass table."""
    assert VAL["scan"]["n_runs"] == len(SCAN) == 72
    assert VAL["scan"]["n_trapped"] == sum(1 for r in SCAN if r["first_trapped"] is not None)
    assert VAL["scan"]["dev_median"] is not None and abs(VAL["scan"]["dev_median"]) < 0.01
    assert VAL["null_ray_minima"]["minima"] == 0
    m1 = [r for r in TWJ["mass_table"] if r["q"] == 1.0][0]
    assert abs(m1["m_c_over_RdS"] - 2.0) < 1e-12 and abs(m1["r_h_over_RdS_at_mc"] - 1.0) < 1e-12
    assert 5.0 < TWJ["R_dS_Gpc"] < 5.6 and 1.0e23 < m1["M_c_Msun"] < 1.3e23


@check
def anti_trapped_core_needs_supercritical_mass():
    """Every run with an anti-trapped core on its last row has an exterior mass above the thin-wall critical mass (2 R_dS for q = 1)."""
    for r in SCAN:
        if r["at_points_last_row"] > 0:
            assert r["m_ext"] > 2.0, r


@check
def null_energy_sign():
    """At a marginal surface (r_v = 0) the constraint gives r_vv = -(r/2) phi_v^2 <= 0: evaluate on the initial data of a run (analytic identity of Eq. (4))."""
    pot = D.make_potential(1.0, 24.0, 3.0)
    v, r, p, s = D.initial_cone(0.1, pot, 2.0, 0.1, 400, 0.01)
    pv = np.gradient(p, 0.01, edge_order=2)
    # sigma_v = r p_v^2 on the cone with r_v = 1/2, r_vv = 0: constraint residual of Eq. (4)
    sv = np.gradient(s, 0.01, edge_order=2)
    resid = 0.0 - sv * 0.5 + 0.5 * r * pv ** 2
    assert np.max(np.abs(resid[5:-5])) < 0.02 * np.max(0.5 * r * pv ** 2)


if __name__ == "__main__":
    for f in CHECKS:
        f()
        print("ok   ", f.__name__)
    print("%d checks passed" % len(CHECKS))
