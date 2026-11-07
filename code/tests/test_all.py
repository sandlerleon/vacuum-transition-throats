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


VE = json.load(open(os.path.join(CODE, "validation_extra.json")))


@check
def finite_tension_series():
    """The small-tension expansion of the critical mass (Eq. 9) reproduces the numerically maximized m_+(R) to O(sigma^4)."""
    for q in (0.5, 1.0, 2.0):
        for s in (0.05, 0.1, 0.3):
            assert abs(TW.m_crit(s, q) - TW.m_crit_series(s, q)) < 2 * s ** 4, (q, s)


@check
def horizon_coincidence_finite_tension():
    """At finite tension the exterior horizon radius at m_c is R_dS {1 - sigma^2 [1/4 - 1/(6(1+q))]/(1+3q)} to second order, so the area coincidence holds only in the limit."""
    for q in (0.5, 1.0, 2.0):
        s = 0.2
        rh = TW.horizon_radius(TW.m_crit(s, q), q)
        pred = 1 - s ** 2 * (0.25 - 1 / (6 * (1 + q))) / (1 + 3 * q)
        assert abs(rh - pred) < 2e-3, (q, rh, pred)
        assert rh < 1.0


@check
def exterior_roots_at_critical_mass():
    """R_dS is a root of the exterior horizon equation at m = (1+q) R_dS for every q > -1; for q < 0 it is the cosmological (larger) root."""
    for q in (-0.5, 0.0, 0.5, 1.0, 2.0):
        roots = TW.exterior_roots(1.0 + q, q)
        assert any(abs(x - 1.0) < 1e-9 for x in roots), (q, roots)
        if q < 0:
            assert max(roots) == roots[-1] and abs(roots[-1] - 1.0) < 1e-9 and roots[0] < 1.0


@check
def no_dynamical_creation_of_expanding_core():
    """None of the runs that start with a normal core develops an anti-trapped point at a later row (persistence, not creation)."""
    c = VE["creation"]
    assert c["n_initially_normal_core"] == 42 and c["n_created_later"] == 0 and c["n_initially_expanding_core"] == 30


@check
def constraints_and_expansions_converge():
    """RMS constraint residuals of the dynamical solutions fall with the grid step at order >= 1.3, and the signs of the null expansions agree between resolutions."""
    for key in ("constraint_orders", "constraint_orders_AT"):
        for k in ("vv", "uu"):
            assert min(VE[key][k]) > 1.3, (key, k, VE[key][k])
    assert VE["expansions"]["normal_core"]["sign_agreement"] and VE["expansions"]["expanding_core"]["sign_agreement"]


@check
def refinement_of_unresolved_runs():
    """Six of the seven runs without a trapped surface at the base resolution form one after refinement; the seventh is reported as unresolved."""
    inc = VE["inconclusive"]
    assert len(inc) == 7
    assert sum(1 for r in inc if any(v["first_trapped"] is not None for k, v in r.items() if k.startswith("refined"))) == 6



@check
def closed_generic_trajectories_collapse():
    """Phase 1: every one of 216 trajectories started at the classical neck (exact constraint, no assumed bounce, with and without dissipation) turns around once and then collapses."""
    R = json.load(open(os.path.join(CODE, "closed_cycle.json")))
    assert len(R) == 216
    assert all(sum(1 for t in o["turns"] if t["kind"] == "bounce") == 0 for o in R)
    assert all(len(o["turns"]) == 1 for o in R)


@check
def closed_exact_symmetric_cycles_exist():
    """Phase 1: exact time-symmetric cycles exist (the start state is recovered after a full cycle) and are hyperbolic: some Floquet multiplier has modulus above 1."""
    B = json.load(open(os.path.join(CODE, "closed_stability_breather.json")))
    assert len(B) == 5
    for b in B:
        assert max(abs(abs(complex(*z)) - 1.0) for z in b["multipliers"]) > 1.0
        assert abs(b["fixed_point_error"][0]) < 1e-9 and abs(b["fixed_point_error"][1]) < 1e-9
    W = json.load(open(os.path.join(CODE, "closed_stability.json")))
    best = [w for w in W if w["mu"] == 10.0 and w["rc"] == 0.05 and abs(w["C_over_Ccrit"] - 0.9) < 1e-9][0]
    assert best["excursion"] > 8.0 and abs(best["fixed_point_error"][0]) < 1e-6 and best["bounces_survived_by_relative_seed_error"]["1e-06"] == 0


@check
def phase1b_periodicity_and_dissipation():
    """Phase 1B: the exact orbits return their full state (a, a_dot, phi, phi_dot) after one period; with dissipation the entropy parameter S = rho_r a^4 grows at every bounce
    (no periodic orbit), the perturbation stretching is not reduced, and no run survives more than 40 bounces."""
    D = json.load(open(os.path.join(CODE, "closed_dissipation.json")))
    for p in D["periodicity"]:
        if p["kind"] == "breather":
            assert abs(p["delta"][0]) < 1e-9 and abs(p["delta"][1]) < 1e-9 and abs(p["delta_over_v"][0]) < 1e-8 and abs(p["delta_over_v"][1]) < 1e-8
            assert abs(p["constraint_start"]) < 1e-11 and abs(p["constraint_end"]) < 1e-11
    base = {}
    for r in D["dissipation"]:
        key = (r["kind"], r["mu"], r["rc"], r["C_over_Ccrit"])
        if r["mode"] == "none":
            base[key] = r
        if r["mode"] != "none" and len(r["S"]) >= 2:
            assert all(b >= a for a, b in zip(r["S"], r["S"][1:])), r["mode"]
    for r in D["dissipation"]:
        key = (r["kind"], r["mu"], r["rc"], r["C_over_Ccrit"])
        if r["mode"] != "none" and r["stretch"] is not None and base[key]["stretch"]:
            assert r["stretch"] > 0.95 * base[key]["stretch"]
    L = json.load(open(os.path.join(CODE, "closed_dissipation_lifetime.json")))
    assert max(x["max_bounces"] for x in L) <= 40


@check
def phase2_junction_conservation_and_counts():
    """Phase 2 step 1: the shell stress of the comoving bridge obeys (sigma R^2)' + p (R^2)' = 0 to finite-difference accuracy, sigma > 0 iff m < m_max, and the stored feasibility map has the reported counts."""
    import bridge_junction as BJ
    B = json.load(open(os.path.join(CODE, "closed_stability_breather.json")))
    o = [b for b in B if b["mu"] == 20.0 and b["rc"] == 0.2][0]
    m, t, Y, acc = BJ.orbit_history(o, n=20001)
    a, ad = Y[0], Y[1]
    chi0, q = 0.7 * np.pi, 0.0
    sn = np.sin(chi0)
    R = a * sn
    m_max = float(np.min(R * (sn ** 2 + (ad * sn) ** 2 + q * R ** 2)))
    sg, p, Rr = BJ.shell(0.5 * m_max, q, chi0, a, ad, acc)
    assert np.all(sg > 0)
    lhs = np.gradient(sg * Rr ** 2, t) + p * np.gradient(Rr ** 2, t)
    scale = np.max(np.abs(np.gradient(sg * Rr ** 2, t)))
    assert np.max(np.abs(lhs[50:-50])) < 2e-3 * scale
    assert np.any(BJ.shell(1.05 * m_max, q, chi0, a, ad, acc)[0] <= 0)
    J = json.load(open(os.path.join(CODE, "bridge_junction.json")))
    assert len(J) == 72 and sum(1 for r in J if r["frac_nec"] > 0) == 15


@check
def phase2_leak_identity_and_negative_result():
    """Phase 2 Stage A: K_eff drifts exactly by -int a^2 Q/3 dt; no leak run with K_eff >= 0.9 beats its zero-leak same-dissipation baseline; no run stabilises (finite-time exponent > 0 everywhere)."""
    import bridge_leak as BL
    B = json.load(open(os.path.join(CODE, "closed_stability_breather.json")))
    o = [b for b in B if b["mu"] == 20.0 and b["rc"] == 0.2][0]
    m = BL.LeakModel(q=o["q"], rc=o["rc"], mu=o["mu"], C=o["C_over_Ccrit"] * BL.CC.C_CRIT, gamma=0.1, eps=10.0, law="diffuse")
    seed = 10 ** o["log10_seed_over_v"] * m.v
    sol, Bn, Tn, end = BL.run_leak(m, BL.start(m, seed), 1.2)
    t = np.linspace(0.0, sol.t[-1], 6001)
    Y = sol.sol(t)
    Q = np.array([m.Qleak(Y[0, i], Y[3, i], Y[4, i]) for i in range(len(t))])
    dK = m.keff(Y[:, -1]) - m.keff(Y[:, 0])
    pred = -np.trapezoid(Y[0] ** 2 * Q / 3.0, t)
    assert abs(dK - pred) < 2e-4 * max(1.0, abs(pred)), (dK, pred)
    R = json.load(open(os.path.join(CODE, "bridge_leak.json")))
    assert len(R) == 122
    for r in R:
        base = [x for x in R if x["orbit"] == r["orbit"] and x["gamma"] == r["gamma"] and x["law"] == "none"][0]
        free = [x for x in R if x["orbit"] == r["orbit"] and x["gamma"] == 0.0][0]
        if r["law"] != "none" and r["keff_min"] >= 0.9:
            assert r["bounces"] <= base["bounces"], r
        assert r["bounces"] <= free["bounces"] + 0 or r["keff_min"] < 0.0 or r["gamma"] >= 1.0, r
        if r.get("lyapunov") is not None:
            assert r["lyapunov"] > 0.0, r


@check
def phase2c_moving_shell_and_elliptic_cycle():
    """Phase 2C: (a) the moving-shell junction obeys the flux conservation law (sigma R^2)' + p (R^2)' = -R^2 (rho + P) cosh(beta) sinh(beta) to 1e-5 and the interior is untouched by the shell;
    (b) the single-field model has an ELLIPTIC exact cold cycle at C/Cc = 0.35 (trace in (-2, 2), det 1) -- no reservoir or bridge involved; (c) the stored reservoir orbit has all |mu| = 1 and a tiny constraint residual."""
    import bridge_moving_shell as MS
    import closed_stability as CS
    import closed_cycle as CC2
    B = json.load(open(os.path.join(CODE, "closed_stability_breather.json")))
    o = [b for b in B if b["mu"] == 20.0 and b["rc"] == 0.2][0]
    per = MS.periodic_shell(o, 0.7, 0.2, 1.0, -1.0)
    assert len(per) == 1 and abs(per[0]["chi0_over_pi"] - 0.6738) < 2e-4
    assert per[0]["diag"]["flux_identity_resid"] < 1e-5 and abs(per[0]["det"] - 1.0) < 1e-6 and all(abs(z - 1.0) < 1e-6 for z in per[0]["multipliers_abs"])
    E = json.load(open(os.path.join(CODE, "single_elliptic.json")))
    assert not any(r["elliptic"] for r in E if r["C_over_Ccrit"] <= 0.345 and r["det"] > 0.9999 and r["det"] < 1.0001)
    assert sum(1 for r in E if r["elliptic"] and 0.35 <= r["C_over_Ccrit"] <= 0.36) >= 3
    r35 = [r for r in E if r["elliptic"] and abs(r["C_over_Ccrit"] - 0.35) < 1e-9][0]
    m = CC2.Model(q=1.0, rc=0.2, mu=20.0, C=0.35 * CC2.C_CRIT)
    s = 10 ** r35["log10_seed_over_v"] * m.v
    J = CS.jacobian(m, s, 0.0)
    assert abs(np.trace(J) + 1.7043) < 1e-3 and abs(np.linalg.det(J) - 1.0) < 1e-6 and abs(np.trace(J)) < 2.0
    P = json.load(open(os.path.join(CODE, "bridge_phase2c.json")))
    assert P["orbit"]["max_abs"] - 1.0 < 1e-6 and P["orbit"]["constraint_max"] < 1e-10
    assert all(v["bounces"] >= 300 for v in P["controls"]["coupling_off"].values())


@check
def phase2d_independent_cycle_and_bridge_geometry():
    """Phase 2D: (a) an independent re-implementation reproduces the stable cycle (orbit, period, neck-map trace, monodromy eigenvalues on the unit circle, tolerance convergence) and the multiplier = -1 point
    at C/Cc = 0.3458 with no local 2-cycle (winding number of P^2 - id equals the index of the fixed point); (b) the two-ended geometry: ordering a > R > r_h, sigma > 0 and the NEC all cycle, rays all cross the shell,
    60-period survival under perturbations, interior identical to the no-bridge baseline, constraint residuals; (c) the stiffest shell is only linearly stable (recorded failure at 1e-3)."""
    import phase2d_cycle_validation as PV
    V = json.load(open(os.path.join(CODE, "phase2d_cycle_validation.json")))
    r = [x for x in V["reproduction"] if abs(x["x"] - 0.35) < 1e-9][0]
    assert abs(r["phi0_over_v"] - 0.22294841) < 2e-8 and abs(r["period"] - 0.596356281) < 1e-8 and r["return_err_full_state"] < 1e-12
    assert abs(r["conformal_period"] - r["period"]) < 1e-10
    assert abs(r["neck_map_trace"] + 1.704348) < 1e-4 and abs(r["neck_map_det"] - 1.0) < 1e-6
    for e in r["monodromy_eigs"]:
        assert abs(abs(complex(*e)) - 1.0) < 1e-6
    assert abs(V["flip"]["C_over_Ccrit_star"] - 0.345844) < 1e-5
    c = PV.Cycle(0.35)
    p = c.orbit(0.2229, span=0.05)
    J = c.jac_neck(p, 0.0)
    assert abs(np.trace(J) + 1.704348) < 1e-4 and abs(np.linalg.det(J) - 1.0) < 1e-6
    W = json.load(open(os.path.join(CODE, "phase2d_winding.json")))
    for rec in W:
        idx = -1.0 if rec["dx"] < 0 else 1.0
        assert all(rec["windings"][k] == idx for k in ("0.002", "0.005", "0.01", "0.02"))
    B = json.load(open(os.path.join(CODE, "phase2d_bridge.json")))
    A = B["analysis"][0]
    assert A["ordering_equator_gt_shell_gt_throat"] and A["sigma_pos"] == 1.0 and A["nec_fraction"] == 1.0 and A["shell_over_throat_min"] > 90 and A["rays_crossed"] == 24
    assert A["interior_trapped_fraction"] == 0.0 and A["exterior_trapped_fraction"] == 0.0 and A["shell_outside_horizon_always"]
    L = [x for x in B["long"] if x["shell_index"] == 0]
    assert len(L) == 8 and all(x["completed"] and x["periods_reached"] >= 59.99 and x["interior_vs_alone_diff"] < 1e-8 and x["R_over_rh_min"] > 80 for x in L)
    L4 = [x for x in B["long"] if x["shell_index"] == 4 and x["pert"] >= 1e-3]
    assert L4 and all(not x["completed"] for x in L4)
    K = json.load(open(os.path.join(CODE, "phase2d_convergence.json")))["residuals"]
    assert K["friedmann"] < 1e-12 and K["gamma_in"] < 1e-12 and K["israel_p"] < 1e-12 and K["mass"] < 1e-10 and K["flux"] < 1e-5


@check
def phase2e_global_junction_and_inhomogeneous_instability():
    """Phase 2E: (a) global junction consistency from raw Christoffel-symbol extrinsic curvatures (residuals ~1e-15, all orientation flags true, sigma and p recovered, p = w sigma);
    (b) the shell exchanges energy with the bulk comparable to its own energy per cycle; (c) linear inhomogeneous perturbations (n >= 3) of the cycle grow by >= 14 per cycle (up to ~1.5e5), the derivation
    passes constraint-propagation / trace-equation checks, the growth matches the fifth-force estimate at large n and disappears at large n when the coupling is removed from the perturbations;
    (d) the homogeneous window survives at the 2/5 and 3/8 resonances but its basin shrinks near the end of the window."""
    J = json.load(open(os.path.join(CODE, "phase2e_junction_global.json")))
    assert max(J["residuals"].values()) < 1e-12 and all(J["flags"].values())
    assert 0.5 < J["exchange_over_shell_energy"] < 2.0 and J["exchange_over_cap_energy"] < 0.05
    F = json.load(open(os.path.join(CODE, "phase2e_perturb_floquet.json")))
    assert all(v["constraint_propagation"] < 1e-6 and v["trace_equation"] < 1e-4 for v in F["checks"].values())
    assert all(r["max_abs"] > 14.0 and abs(r["det"] - 1.0) < 1e-4 for r in F["modes"]) and F["modes"][-1]["max_abs"] > 1e5
    Mx = json.load(open(os.path.join(CODE, "phase2e_perturb_mechanism.json")))
    c60 = [r for r in Mx["compare"] if r["n"] == 60][0]["ln_multiplier"]
    assert 0.9 < c60 / Mx["estimate_log_growth"] < 1.0
    assert [max(r["multipliers_abs"]) for r in Mx["control"] if r["n"] == 60][0] < 2.0
    W = json.load(open(os.path.join(CODE, "phase2e_window.json")))
    by = {r["name"]: r for r in W["points"]}
    for name in ("2/5", "3/8", "centre", "near_flip"):
        assert all(min(s["bounces"]) == 150 for s in by[name]["survival"])
    assert min(by["near_end"]["survival"][-1]["bounces"]) == 0


@check
def phase2e_manuscript_builds_from_results():
    """The Phase 2E manuscript is built from the JSON results: title as agreed, abstract <= 250 words, all references cited, and the key numbers in the text match the result files."""
    import phase2e_paper_content as PC
    assert PC.TITLE == "A Nontraversable Einstein–Rosen Bridge to a Nonsingular Cyclic Cosmological Region: Junction Dynamics, Stability, and Causal Structure"
    abstract = [b[1] for b in PC.B if b[0] == "p0"][0]
    assert len(abstract.split()) <= 250 and "14" in abstract
    txt = " ".join(b[1] for b in PC.B if len(b) > 1 and isinstance(b[1], str))
    assert "0.222948413319" in txt or "0.222948413319" in " ".join(str(c) for b in PC.B if b[0] == "table" for row in b[3] for c in row)
    assert "%.1f" % PC.PM["estimate_log_growth"] in txt
    import docx_engine as DE
    cit = DE.Citations(PC.LIB)
    for b in PC.B:
        for x in b[1:]:
            if isinstance(x, str):
                cit.scan(x)
    cit.entries()
    assert os.path.exists(os.path.join(CODE, "..", "manuscript", "Nontraversable_ER_Bridge_Cyclic_Region.docx"))


@check
def phase2f_stabilization_attempts_fail():
    """Phase 2F: the generalised derivation reproduces Phase 2E at K = 0; continuation in the coupling, the potential scale, the sextic hilltop shape and a coupled fluid with pressure never gives bounded inhomogeneous
    perturbations: n = 3 growth stays >= 8.6 per cycle, the worst case over n stays >= 3e3, pressure only removes modes above n_c ~ g a / c_s (WKB criterion validated), the shell exchange is ~92% dust."""
    P = json.load(open(os.path.join(CODE, "phase2f_pressure.json")))
    assert abs(P[0]["growth"]["3"] - 14.3398) < 1e-3 and abs(P[0]["growth"]["60"] / 154820.5 - 1) < 1e-4
    for r in P:
        assert abs(r["growth"]["3"] / 14.34 - 1) < 0.01 and abs(r["growth"]["6"] / 169.0 - 1) < 0.01 and abs(r["growth"]["12"] / 2050.0 - 1) < 0.05
    assert abs(P[-1]["growth"]["240"] - 1.0) < 1e-3 and P[-1]["growth"]["120"] > 300
    CU = json.load(open(os.path.join(CODE, "phase2f_continuation_up_mu20.json"))) + json.load(open(os.path.join(CODE, "phase2f_continuation_down_mu20.json")))
    assert len(CU) > 40 and all(abs(r["trace"] + 1.7043) < 1e-3 and 14.0 <= r["growth"]["3"] <= 17.0 and r["growth"]["60"] > 8e4 for r in CU)
    MU = json.load(open(os.path.join(CODE, "phase2f_muscan_up_rc0.2_valid.json")))
    assert min(r["growth"]["3"] for r in MU) > 8.6 and min(r["growth"]["12"] for r in MU) > 170
    WC = json.load(open(os.path.join(CODE, "phase2f_worstcase.json")))
    assert min(r["worst_G"] for r in WC) > 3e3 and min(r["min_G"] for r in WC) > 8.6
    PO = json.load(open(os.path.join(CODE, "phase2f_potential.json")))
    assert all(12.0 < r["growth"]["3"] < 17.5 and r["growth"]["12"] > 1e3 and r["growth"]["60"] > 3e4 for r in PO["wider"] + PO["narrower"])
    A = json.load(open(os.path.join(CODE, "phase2f_analysis.json")))
    hi = [w for w in A["wkb"] if w["family"] == "pressure" and w["n"] in (120, 240) and w["K"] in (0.04, 0.08)]
    assert hi and all(abs(w["S5"] - w["lnG"]) < 0.6 for w in hi)
    assert A["shell_split"]["dust_fraction"] > 0.9 and abs(A["shell_split"]["net_dust"]) < 1e-10


if __name__ == "__main__":
    for f in CHECKS:
        f()
        print("ok   ", f.__name__)
    print("%d checks passed" % len(CHECKS))
