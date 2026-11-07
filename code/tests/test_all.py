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


@check
def phase2g_cap_spectrum_and_boundary_conditions():
    """Phase 2G-A: no-flux boundary conditions on the finite cap (Q' = 0 at the shell, comoving shell): for the AdS exterior that matches the model, NEC comoving shells exist only for chi0 <~ 0.6 pi;
    the cap spectrum nu_j(chi0) is quantised by tan(nu chi0) = nu tan(chi0); the growth G(nu) is > 1 for every nu > 1 (-> 1 only at nu = 1) and increases for nu > 2, so no boundary condition can stabilise the cap;
    the lowest cap mode keeps G >= 1.9 and the worst case is unchanged."""
    import numpy as np
    N = json.load(open(os.path.join(CODE, "phase2g_comoving_nec.json")))
    assert N["n_nec_full"] > 1000
    S = json.load(open(os.path.join(CODE, "phase2g_shell_data.json")))
    assert max(r["chi0_over_pi"] for r in S["shells"]) <= 0.6 + 1e-9 and all(r["w_range"][1] > r["w_range"][0] + 0.3 for r in S["shells"])
    assert S["phi_range_over_v"][0] > 0.2 and S["phi_range_over_v"][1] < 0.7 and min(S["wall_tension"]) < 0.4 and min(r["sigma_range"][0] for r in S["shells"]) > 2.0
    C = json.load(open(os.path.join(CODE, "phase2g_cap_spectrum.json")))
    for rec in C:
        chi0 = rec["chi0_over_pi"] * np.pi
        for nu in rec["nu"]:
            assert abs(np.tan(nu * chi0) - nu * np.tan(chi0)) < 1e-6 * max(1.0, abs(nu * np.tan(chi0))) or abs(np.cos(nu * chi0)) < 1e-3
        g = [x for x in rec["growth"] if x is not None]
        assert min(g) > 1.85 and max(g) > 7e4
    low = {round(r["chi0_over_pi"], 2): r["growth"][0] for r in C}
    assert abs(low[0.55] - 10.43) < 0.05 and abs(low[0.6] - 7.78) < 0.05
    L = json.load(open(os.path.join(CODE, "phase2g_nu_low.json"))); H = json.load(open(os.path.join(CODE, "phase2g_nu_curve.json")))
    assert all(r["max_abs"] > 1.0 for r in L) and L[0]["max_abs"] < 1.1
    hv = [r["growth"] for r in H if "growth" in r and r["nu"] >= 2.01]
    assert all(b > a for a, b in zip(hv, hv[1:])) and hv[0] > 1.9 and hv[-1] > 80


@check
def phase3_baseline_and_effective_fluid_feasibility():
    """Phase 3A-3C: exact dust collapse baselines (OS, LTB) and the stored scalar scan form no daughter region; the spin-fluid (EC-type) effective fluid is thermodynamically consistent and coincides with the LQC-type
    effective fluid for dust; a homogeneous bounce exists; vacuum-exterior matching fails (R_bounce < m in most cases, mass mismatch in all); sound speed is imaginary below n* and superluminal above;
    perturbation transfer matrices are singular at n* and the spin term adds a k-dependent growth on top of the dust-collapse growth."""
    A = json.load(open(os.path.join(CODE, "phase3a_baseline.json")))
    os_ = A["oppenheimer_snyder"]
    assert abs(os_["t_singular_num"] / os_["t_singular_exact"] - 1) < 1e-6 and abs(os_["t_horizon_num"] / os_["t_horizon_exact"] - 1) < 1e-6 and os_["mass_spread"] < 1e-12 and not os_["reexpands"]
    L = A["ltb"]
    assert L["n_trapped"] == L["n_shells"] and L["central_first"] and not L["shell_crossing_before_collapse"] and L["first_integral_residual"] < 1e-3
    assert A["scalar_scan"]["n_trapped"] == 65 and A["scalar_scan"]["creation"]["n_created_later"] == 0 and A["scalar_scan"]["creation"]["n_initially_normal_core"] == 42
    f1, f2 = A["scalar_fresh"]
    assert f1["anti_trapped_points_max"] == 0 and f2["anti_trapped_points_max"] == 0 and abs(f1["u_first_trapped"] - f2["u_first_trapped"]) <= 1.01 * f1["h"]
    E = json.load(open(os.path.join(CODE, "phase3c_ec_homogeneous.json")))
    assert E["fluid"]["thermodynamic_identity_residual"] < 1e-9 and E["lqc_equivalence_max_diff"] < 1e-8
    assert E["fluid"]["cs2_below_nstar_max"] < 0 and E["fluid"]["cs2_above_nstar_min"] > 1.0 and E["fluid"]["NEC_violated_above_nstar"] and E["fluid"]["NEC_holds_below_nstar"]
    c = E["cycle"]
    assert c["return_error"] < 1e-12 and c["friedmann_residual"] < 1e-9 and c["bounce_above_nstar"]
    mt = E["matching"]
    assert mt["counts"]["horizon_forms_bounce_inside_exterior_horizon"] >= 20 and E["theorem_check"]["n_cases_bounce_inside"] >= 25
    assert all(abs(r["m_eff_max_over_m"] - r["m_eff_min_over_m"]) > 0.01 for r in mt["rows"])
    W = json.load(open(os.path.join(CODE, "phase3c_ec_window.json")))
    for w in W["windows"]:
        for md in w["modes"]:
            if md["n"] <= 6:
                assert abs(md["det"] - 1.0) < 1e-2 and md["ln_growth"] > 10.0
    Pt = json.load(open(os.path.join(CODE, "phase3c_ec_perturb.json")))
    assert abs(np.log(Pt["modes"][0]["max_multiplier"]) - 22.57) < 0.2 and abs(Pt["modes"][0]["det"] - 1.0) > 1e3 and Pt["frac_imaginary_cs"] > 0.9
    K = json.load(open(os.path.join(CODE, "phase3c_ec_control.json")))
    for row in K:
        d = {}
        for cs in row["cases"]:
            d.setdefault(cs["n"], {})[cs["B"] > 0] = cs["ln"]
        extra = [d[n][True] - d[n][False] for n in (3, 6, 12)]
        assert extra[0] < extra[1] < extra[2] and extra[2] > 2.0


@check
def phase3c2_rotation_spin_magnetic_feasibility():
    """Phase 3C.2/3C.3: verified kinematic identities (effective-fluid link, spin-vorticity perfect square, circulation law); exact Bianchi-I axisymmetric model (constraint, conservation, GR limit, analytic bounce volume);
    rotation proxy and magnetic field RAISE the critical spin strength (no help), rotating/magnetized GR without torsion always singular; local MHD dispersion relation of the effective fluid (factorisation, instability for
    every oblique angle, perpendicular stability only above a field threshold); the field at the bounce is far below the regularisation threshold."""
    import numpy as np
    K = json.load(open(os.path.join(CODE, "phase3c2_kinematics.json")))
    assert K["effective_fluid_identity"] == "0" and K["perfect_square_difference"] == "0" and K["omega_dot_over_omega"] == "-2*H_perp"
    assert 1.2 < K["fit_rotation"]["exponent"] < 1.5 and 1.4 < K["fit_magnetic"]["exponent"] < 1.7
    Bj = json.load(open(os.path.join(CODE, "phase3c2_bianchi.json")))
    c = Bj["checks"]
    assert c["constraint_max"] < 1e-8 and c["conservation_rel"] < 1e-5 and c["GR_limit_outcome"] == "singular" and c["sigma2_V2_spread"] < 1e-9 and abs(c["numeric_Vb"] / c["analytic_Vb"] - 1) < 1e-6
    r = {x["case"].strip(): x["B_critical"] for x in Bj["critical"]}
    iso = r["A  nonrotating EC (isotropic start)"]
    rot = [r["C  rotating EC v0^2 = 0.01"], r["C  rotating EC v0^2 = 0.05"], r["C  rotating EC v0^2 = 0.2"]]
    mag = [r["D  magnetized nonrotating, rho_B0 = 0.01"], r["D  magnetized nonrotating, rho_B0 = 0.05"], r["D  magnetized nonrotating, rho_B0 = 0.2"]]
    assert iso < 1e-6 and all(b > 1e3 * iso for b in rot + mag) and rot == sorted(rot) and mag == sorted(mag)
    assert r["E  rotating + magnetized, 0.2 / 0.2"] > max(rot[2], mag[2]) and all(x["outcome"] == "singular" for x in Bj["model_B_no_torsion"])
    M = json.load(open(os.path.join(CODE, "phase3c3_mhd.json")))
    assert M["factorization_check"] == "0" and M["fraction_unstable"] > 0.8
    for cs2, vA2 in M["perp_stable"]:
        assert vA2 > -cs2 / (1 - cs2)
    for sc in M["scan"]:
        if sc["cs2"] < 0 and sc["vA2"] > 0 and abs(sc["psi"] - np.pi / 2) > 1e-6:
            assert min(sc["omega2_over_k2"]) < 0
    F = json.load(open(os.path.join(CODE, "phase3c3_bounce_fields.json")))
    assert len(F) >= 12 and max(x["rho_B_over_pspin"] for x in F) < 0.2 and not any(x["perp_fast_real"] for x in F)


@check
def phase3c4_polarized_spin_cross_terms():
    """Phase 3C.4: polarized-spin cross terms derived from the paper's effective stress (Eq. 77): sign convention fixed by the paper's Raychaudhuri cross term; energy 2 s omega, equatorial pressure s omega, trace 0;
    shear source proportional to omega (omega + kappa s / 2); the reduced homogeneous axisymmetric system is NOT closed (stress-energy conservation mismatch of exactly beta_dot with separately conserved s ~ n and omega ~ a_perp^-2);
    what-if bound: cancelling a fraction f of the rotation-induced stress lowers B_c monotonically."""
    P = json.load(open(os.path.join(CODE, "phase3c4_polarized.json")))
    assert P["eta_fixed_by_paper"] == ["1"]
    assert "2*omega*s" in P["cross_stress"]["T_uu"].replace(" ", "") and P["cross_stress"]["trace"] == "0"
    assert P["shear_source_xx"].replace(" ", "") == "omega*(kappa*s+2*omega)/6"
    assert P["raychaudhuri_perfect_square"].replace(" ", "") in ("(kappa*s-2*omega)**2/2", "(2*omega-kappa*s)**2/2")
    assert P["conservation_required"].replace(" ", "") == "-alpha_dot-3*beta_dot" and P["conservation_mismatch"] == "beta_dot"
    bc = [r["B_critical"] for r in P["what_if"]]
    assert bc == sorted(bc, reverse=True) and bc[0] > 0.03 and bc[-1] < 1e-4


@check
def phase3c5_spin_transport_and_closure():
    """Phase 3C.5: Frenkel-preserving spin transport is forced (lambda = 1) and gives d ln s/dtau = -Theta for any acceleration; [D_a,D_b]f = -omega_ab fdot on an explicit rotating-frame metric (so rest-space
    homogeneity is incompatible with rotation plus time dependence); vector-constraint sign k = 2; the three laws (s ~ n, Helmholtz, cross-term energy balance) force beta_dot = 0; the cancellation surface is
    left at rate -alpha_dot s; the shear is not removed when its source vanishes."""
    C = json.load(open(os.path.join(CODE, "phase3c5_closure.json")))
    assert C["A_frenkel_residual"] == ["0"] * 4 and C["A_s_squared"] == "s**2"
    assert C["B_lambda_forced"] == ["{lambda: 1}"] and C["B_dlns_over_dtau"] == "-Theta"
    assert all(abs(r + 1.0) < 1e-9 for r in C["C_commutator_ratio"]) and len(C["C_commutator_ratio"]) == 3
    assert all(abs(k - 2.0) < 1e-9 for k in C["C_vector_constraint_k"])
    assert C["D_mismatch"] == "beta_dot" and C["D_preservation_derivative"].replace(" ", "") == "-alpha_dot*s"
    assert C["D_shear_fraction_at_alpha_dot_zero"] == "1/4"
    assert all(r["max_abs_deviation"] > 1.0 for r in C["E_bianchi_ratio"] if r["outcome"] == "bounce")
    assert C["E_dust_sigma2V2_spread"] < 1e-9 and C["E_dust_outcome"] == "singular"


@check
def phase3d0_einstein_cartan_dirac_feasibility():
    """Phase 3D.0: ECD contact interaction L_eff = -(3 kappa/16) A_a A^a derived from explicit gamma matrices (|c_K| = 1/4, vector/tensor torsion decouple, R_quad = 6 a.a); Hartree-Fock energies of fermion states all have the
    Weyssenhoff form eps = p = -B n^2 (fully polarized NR gas: 0, Pauli); barotropic bounce necessarily has mu = 0 crossing and negative compressibility (random-EOS check); exact baseline numbers."""
    E = json.load(open(os.path.join(CODE, "phase3d0_ecd.json")))
    assert abs(E["A_cK_abs"] - 0.25) < 1e-12 and E["A_nonaxial_decoupling_residual"] < 1e-12
    assert E["A_R_quadratic_over_aa"] == "6" and E["A_pairing_over_aA"] == "-6"
    assert E["A_Leff_coefficient_times_AaAa"].replace(" ", "") == "-3*kappa/16"
    hf = E["B_hf_eps_over_lam_n2"]
    assert abs(hf["unpolarized massless (g=2, particles)"] + 0.5) < 1e-9 and abs(hf["single-helicity Weyl (g=1)"] + 2.0) < 1e-9
    assert abs(hf["unpolarized nonrelativistic"] + 1.5) < 1e-9 and abs(hf["fully polarized NR"]) < 1e-9
    for v in E["C_baseline"].values():
        assert abs(v["n_bounce_numeric"] / v["n_bounce_analytic"] - 1) < 1e-6 and v["constraint_residual"] < 1e-9
        assert abs(v["nstar_over_nb"] - (2 / 3) ** 1.5) < 1e-9 and abs(v["nsp_over_nb"] - (2 / 9) ** 1.5) < 1e-9
        assert v["t_mu_zero"] is not None and v["t_spinodal"] is not None and v["t_spinodal"] < v["t_to_bounce"]
    r = E["C_random_theorem_check"]
    assert r["bounces"] > 100 and r["mu_cross"] == r["bounces"] and r["neg_compress"] == r["bounces"]


@check
def phase3d1_wave_transmission_and_backreaction_requirements():
    """Phase 3D.1: Simpson-Visser black-bounce benchmark. Radial NEC T_kk = -2 a^2/(kappa R^4) < 0 at every r; geometry regular; E=1 radial geodesics cross the throat for a/m = 3, 2 and 1 (two-way, extremal, one-way through a
    horizon); scalar waves: flux conserved to <1e-5 (massless), box-size converged to 1e-6, wormhole throat is a high-pass filter (T_0 << Schwarzschild greybody at low omega) and ~1 at high omega; massive field reflected near threshold
    for mu a = 0.3; quantum-inequality estimate needs N ~ 5.6 (a/l_P)^2 species."""
    W = json.load(open(os.path.join(CODE, "phase3d1_wave.json")))
    assert all(v["Tkk_radial_max"] < 0 for v in W["A_energy_conditions"].values())
    assert abs(W["A_energy_conditions"]["3.0"]["Tkk_at_throat"] + 2 / (8 * 3.141592653589793 * 9)) < 1e-12
    assert max(W["A_kretschmann_r0"].values()) < 2e3 and all(v["radial_geodesic_crosses"] for v in W["A_causal"].values())
    rows = [r for r in W["B_scattering"] if r["case"] != "Ellis m=0 a=3"]
    assert max(abs(r["flux_residual"]) for r in W["B_scattering"]) < 1e-5
    pick = lambda case, om: [r["T"] for r in W["B_scattering"] if r["case"] == case and r["ell"] == 0 and abs(r["omega"] - om) < 1e-12][0]
    assert pick("wormhole a/m=3", 0.05) < 0.01 < 0.05 < pick("Schwarzschild a=0", 0.05) and pick("wormhole a/m=3", 0.8) > 0.999
    assert 0.8 < pick("black bounce a/m=1", 0.1) / pick("Schwarzschild a=0", 0.1) < 1.3
    conv = [c["T"] for c in W["B_convergence"]["box"]]
    assert max(conv) - min(conv) < 1e-5
    ms = [r for r in W["B_massive"] if abs(r["mu_a"] - 0.3) < 1e-9]
    assert ms[0]["T"] < 0.05 and ms[-1]["T"] > 0.99
    assert abs(W["C_QI_species_coefficient"] - 5.59) < 0.05 and W["C_volume"]["3.0"]["volume_integral_rho_plus_pr"] < -6


@check
def phase3d2_freefall_transport_tidal_and_causal():
    """Phase 3D.2: radial free fall crosses the minimum-radius surface (E^2 drift < 1e-9), bound fall oscillates through a chain of trapped / anti-trapped / exterior regions; critical L for crossing equals the
    Schwarzschild barrier value (4.0 at E=1) for a/m <= 3 and is throat-limited for a/m = 6; angular momentum never helps; tidal tensor finite (K_th finite at the horizon, max ~ m/a^3), transverse Jacobi field equals R(tau) to 1e-10;
    anti-trapped region has positive volume expansion; T_kk integral = -a/4."""
    F = json.load(open(os.path.join(CODE, "phase3d2_freefall.json")))
    for v in F["A_radial"].values():
        assert v["conserved_E_drift"] < 1e-9 and v["tau_first_zero"] is not None and v["R_min_along"] < 3.01
    crit = F["B_critical_L"]
    for k_, v in crit.items():
        assert v["direct_integration"]["L=0.98 Lc"] == "crosses throat" and v["direct_integration"]["L=1.02 Lc"] == "scattered back"
    assert abs(crit["a/m=1.0 E=1.00"]["L_critical"] - 4.0) < 1e-6 and abs(crit["a/m=6.0 E=1.00"]["L_critical"] - 6 * 0.5 ** 0.5 * 1.0) > 0.0
    assert abs(crit["a/m=6.0 E=1.00"]["L_critical"] - 4.2426) < 1e-3
    C = F["C_jacobi"]
    assert C["transverse_vs_R_max_dev"] < 1e-9 and C["zero_crossings_transverse"] == 0 and abs(C["R_min"] - 1.0) < 1e-3
    assert F["C_Kth_at_horizon_a1m1"].replace(" ", "") == "3/32-E**2/16"
    assert abs(F["C_tidal_scaling"]["1.0"]["K11_abs_max"] - 1.0) < 1e-6
    seq = [s_[0] for s_ in F["D_region_sequence"]]
    assert seq[:4] == ["exterior (r>0)", "trapped (R decreasing)", "anti-trapped (R increasing)", "exterior (r<0)"]
    assert all(v[0] > 0 for v in F["D_theta_antitrapped"].values()) and F["D_second_exterior"]["returns_to_horizon_branch"]


@check
def phase3d3_source_census_and_goldilocks_closure():
    """Phase 3D.3: the NEC-violating sector of the benchmark is a phantom scalar with phi' = sqrt(2/kappa) a/R^2 (G^r_r - G^t_t = -2 f a^2/R^4); free-quantum-field Goldilocks map is empty for every a (supply N >= 5.59 a^2 vs
    species-cutoff validity N <= a^2; margin 5.59); ECD spinodal k_F >= 1.0 M_P; no census source passes all hard constraints H1-H7."""
    G = json.load(open(os.path.join(CODE, "phase3d3_goldilocks.json")))
    assert G["A_Grr_minus_Gtt_residual"] == "0" and abs(float(G["A_phi_profile_residual"].split("*")[0])) < 1e-12
    mp = G["B_goldilocks_map"]
    assert all(v["admissible_grid_points"] == 0 for v in mp.values())
    assert abs(mp["QI a/m=6 | cutoff M_P"]["closure_margin_k_over_Mp2"] - 5.59) < 0.05
    assert G["B_slack_fraction_nonempty"]["qi_a_over_m_6"] < 0.05 and G["B_slack_fraction_nonempty"]["qi_a_over_m_inf"] < 0.01
    assert all(v["over_cutoff_MP"] > 0.99 for v in G["B_ecd_vs_cutoff"].values())
    assert G["C_sources_passing_all_hard_constraints"] == []
    assert G["C_census"]["S2 phantom scalar + magnetic NED (Bronnikov-Walia source of the benchmark)"]["H3"] == "FAIL"


@check
def phase3d4_propertime_and_instability_window():
    """Phase 3D.4: invariant curvature profile validated (3D.1 values, Schwarzschild limit); proper transit time through K >= eps K(0) is finite, scales ~ a and falls ~1/E; redshift formula from four-momenta verified
    (observed frequency 1/(E + sqrt(E^2 - f)), regular at the horizon); Killing time diverges logarithmically near the horizon; test-field tachyonic instability exists above xi_c (computed); Killing-time residence has a
    light-crossing floor independent of E, so fast transit cannot beat it."""
    T = json.load(open(os.path.join(CODE, "phase3d4_propertime.json")))
    for k_, (x, y) in T["A_validation"].items():
        assert abs(x / y - 1) < 1e-9
    assert all(v["peak_at_throat"] for v in T["A_profile"].values())
    s = T["B_scaling_dtau_vs_a"]
    assert s["0.5"] < s["1.0"] < s["1.9"] < s["3.0"] < s["6.0"]
    m = T["B_transit_maps"]
    assert m["a/m=3.0 eps=0.1 E=5 L/Lc=0"]["dtau"] < 0.2 * m["a/m=3.0 eps=0.1 E=1 L/Lc=0"]["dtau"]
    hz = T["B_horizon_clock"]
    assert 90 < hz["tau_R0_to_horizon"] < 100 and hz["killing_time_to_r_h_plus_delta"]["1e-06"] - hz["killing_time_to_r_h_plus_delta"]["0.0001"] > 5
    g = T["D_tachyonic_test_field"]
    assert g["wormhole a/m=3"]["xi_critical"] > 1.0 and g["wormhole a/m=6"]["xi_critical"] < g["wormhole a/m=3"]["xi_critical"]
    assert all(r["gamma_killing"] == 0.0 for r in g["wormhole a/m=3"]["rows"] if r["xi"] <= 1.0)
    fl = T["D_killing_time_floor"]
    assert all(v["light_crossing_killing_time"] > 10 for v in fl.values())
    assert all(r["dt_killing"] >= fl[r["case"]]["light_crossing_killing_time"] - 1e-6 for r in T["D_amplification"])


@check
def phase3e0_string_inspired_homogeneous_baseline():
    """Phase 3E.0: heterotic tree-level string-frame minisuperspace (dilaton + Gauss-Bonnet). Identities verified by brute-force curvature; dilaton-vacuum limit reproduced; flat collapse has no bounce (H=0 forces Phi'=0, addot=0);
    closed collapse has no bounce (H=0 would need Phi'^2<0); shear-sector kinetic coefficient 1 - 8cH Phi'; open (k=-1) string-frame bounces exist but all of them also occur without the GB term and none is an Einstein-frame bounce."""
    S = json.load(open(os.path.join(CODE, "phase3e0_string.json")))
    assert S["A_GB_total_derivative_residual"] == "0" and S["A_bruteforce_R_residual"] == "0" and S["A_bruteforce_G_residual"] == "0"
    assert S["B_prebigbang_branches"] == ["1"]
    assert S["B_flat_constraint_H0"].replace(" ", "") == "4*A**3*Pd**2" and S["B_flat_H0_addot"] == ["A*Pdd"] and S["B_flat_H0_Phidd"] == ["3*addot/(2*A)"]
    assert S["B_closed_constraint_H0"].replace(" ", "") == "4*A**3*Pd**2+6*A"
    assert S["C_Kshear_normalised"].replace(" ", "") == "32*H*Pd*c-4"
    D = S["D_summary"]
    assert D["0"]["string_frame_bounces"] == 0 and D["1"]["string_frame_bounces"] == 0
    assert D["-1"]["string_frame_bounces"] > 100 and D["-1"]["string_bounces_also_without_GB"] == D["-1"]["string_frame_bounces"] and D["-1"]["einstein_frame_bounces"] == 0
    assert D["-1"]["gs_max_over_bounces"] > 1.0 and D["-1"]["ghost_runs"] > 0


@check
def phase3e1_field_dependent_gb_nested_throats_and_speeds():
    """Phase 3E.1: for arbitrary xi(phi,chi), V and two canonical scalars, flat collapse has addot/a = -(kappa/2)(phi'^2+chi'^2) <= 0 at H=0 (no bounce); closed bounce needs 8 lam xi'' phi'^2 comparable to kappa a^2 phi'^2 and the
    opposite coupling sign, and then (scan) eps_T >= 4, F_T < 0 (c_T^2 < 0) in every bounce; G^r_r - G^t_t = -2 f R''/R; nested throats need two kinetic-sign changes and worsen the test-field instability; Horndeski-mapped
    G_T, F_T verified, G_T independently derived from Bianchi I."""
    E = json.load(open(os.path.join(CODE, "phase3e1_corrections.json")))
    assert E["A_flat_residual"] == "0" and E["A_addot_over_a_H0_flat"].replace(" ", "") == "kappa*(-Cd**2-Pd**2+2*V(Phi,Chi))/4"
    sm = E["A_scan_summary"]
    assert sm["k=0 lam+"]["bounce"] == 0 and sm["k=0 lam-"]["bounce"] == 0 and sm["k=1 lam+"]["bounce"] == 0
    s1 = sm["k=1 lam-"]
    assert s1["bounce"] > 50 and s1["bounce_epsT_min"] > 1.0 and s1["bounce_with_epsT_below_0p3"] == 0 and s1["bounce_with_ghost_or_gradient_instability"] == s1["bounce"]
    assert E["B_Grr_minus_Gtt"].replace(" ", "") == "-2*f(r)*Derivative(R(r),(r,2))/R(r)" and E["B_ricci_f1_check"] == "0"
    nz = E["B_nested_profile"]["nested (two throats + belly)"]
    assert len(nz["local_minima_of_R"]) == 2 and len(nz["local_maxima_of_R"]) == 1 and len(nz["inflection_points"]) == 2
    assert nz["nec_violating_budget_x_kappa"] > 2 * E["B_nested_profile"]["single throat"]["nec_violating_budget_x_kappa"]
    ti = E["B_test_field_instability"]
    assert ti["nested a=1, beta=1.5, w=0.8"]["xi_critical"] < ti["single throat a=1"]["xi_critical"] and ti["nested a=1, beta=3, w=0.6"]["xi_critical"] < ti["nested a=1, beta=1.5, w=0.8"]["xi_critical"]
    g12 = lambda name: [r["gamma"] for r in ti[name]["rows"] if r["xi"] == 12.0][0]
    assert g12("nested a=1, beta=1.5, w=0.8") > g12("single throat a=1")
    assert E["C_GT_residual"] == "0" and E["C_FT_residual"] == "0" and E["C_GT_bianchi_residual"] == "0"


@check
def phase3f0_causal_fluid_feasibility():
    """Phase 3F.0: (0) the Phase 3E.1 EsGB equations equal the published Guo-Schwarz form exactly (E_N = a^3/2 GS_F, E_a = 3a^2 GS_R + 3a^2/2 GS_F, E_phi = -a^3 GS_S), Noether identity residual 0, coefficient of Hdot at H=0 is 2,
    Hdot(H=0) = -phid^2/2; (1) causal bulk viscosity only lowers the Lemaitre bounce mass M_c; (2) shear viscosity damps Sigma^2 in contraction iff c^2 K > 3(1-w)/(4(1+w)) (12/12 admissible runs, inviscid radiation anisotropy grows 80x);
    (3) the dissipative stress keeps the null energy condition (min ratio > 0.9)."""
    F = json.load(open(os.path.join(CODE, "phase3f0_fluid.json")))
    assert F["0_noether_residuals"][0] == "0" and F["0_EN_over_GSF"] == "A**3/2" and F["0_Ephi_over_GSS"] == "-A**3"
    assert F["0_Ea_decomposition"] == ["{c1: 3*A**2, c2: 3*A**2/2}"] and F["0_coeff_Hd_at_H0"] == "2" and F["0_Hdot_at_H0"].replace(" ", "") == "-Pdot**2/2"
    Mc = F["1_critical_dust_mass"]
    assert abs(Mc["K=0.1 c_b2=0"] - Mc["K=1 c_b2=0"]) < 1e-6 and 14.0 < Mc["K=1 c_b2=0"] < 14.3
    assert all(v < Mc["K=1 c_b2=0"] for k_, v in Mc.items() if "c_b2=0" not in k_ or "c_b2=0." in k_ or "c_b2=1" in k_)
    assert Mc["K=1 c_b2=1"] < Mc["K=1 c_b2=0.25"] < Mc["K=1 c_b2=0"] and Mc["K=0.1 c_b2=1"] > Mc["K=1 c_b2=1"]
    assert F["1_example_run"]["bounce"] and F["1_example_run"]["Pi_min"] >= 0.0
    assert F["2_summary_admissible_agreement"] == F["2_summary"]["admissible_runs"] == 12
    assert F["2_inviscid_baseline"]["0.3333333333333333"]["Sigma2_end"] > 50 * F["2_inviscid_baseline"]["0.3333333333333333"]["Sigma2_start"]
    sh = F["2_shear_scan"]
    assert all(r["nec_ratio_min"] > 0.9 and r["constraint_residual"] < 1e-4 for r in sh)
    assert [r["damped"] for r in sh if r["w"] < 0.5 and r["K"] <= 1.0 and r["c2"] <= 1.0] == [False, False, False, False, False, True]


if __name__ == "__main__":
    for f in CHECKS:
        f()
        print("ok   ", f.__name__)
    print("%d checks passed" % len(CHECKS))
