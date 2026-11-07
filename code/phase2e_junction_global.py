# -*- coding: utf-8 -*-
"""Phase 2E, Test A: global junction consistency, computed from RAW geometry (Christoffel symbols in each side's own coordinates), independent of the shortcut formulas of bridge_moving_shell.py.

Interior (cap):   g_in = diag(-1, a(t)^2, a^2 sin^2 chi, a^2 sin^2 chi sin^2 th)                       coordinates (t, chi, th, ph)
Exterior (kept):  g_out = -f ds^2 + dR^2/f + R^2 dOmega^2,  f = 1 - m/R (+ q R^2)   region III, future-directed Killing time s    coordinates (s, R, th, ph)
Shell worldline: interior (t(tau), chi(tau)), exterior (s(tau), R(tau)).  Unit normals n point FROM the cap TO the kept exterior.
Checks: (1) induced metric h_tt = -1, h_thth = R^2 on both sides; (2) normals unit, orthogonal to u, with n^chi > 0 (inside) and n^R < 0 (outside: toward the throat);
(3) time orientation: dt/dtau > 0 and ds/dtau > 0; (4) K_ab recomputed from Christoffels; S^a_b = -([K^a_b] - delta^a_b [K]); compare S^tau_tau = -sigma and S^th_th = p with the values of the
shell equation of motion (p = w sigma) -- an independent test of the equation of motion and of every orientation; (5) shell stays in region III (R > r_h), kept side contains the throat;
(6) curvature invariants at shell and throat; (7) energy: dust+scalar energy exchanged with the shell per cycle.

    python phase2e_junction_global.py          (writes phase2e_junction_global.json)
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp

import bridge_moving_shell as MS
from phase2d_shell_search import Interior, setup

HERE = os.path.dirname(os.path.abspath(__file__))


def raw_K_inside(a, ad, add, chi, beta, bp):
    """K_ab of the shell from the interior side, from Christoffel symbols of g_in = diag(-1, a^2, a^2 sin^2 chi, ...)."""
    ch, sh = np.cosh(beta), np.sinh(beta)
    u = np.array([ch, sh / a])                                   # (t', chi')
    n_up = np.array([sh, ch / a])                                # outward (increasing chi)
    g = np.diag([-1.0, a * a])
    n_dn = g @ n_up
    # Christoffels for (t, chi) block and theta-theta
    G = np.zeros((2, 2, 2))                                      # G[mu, alpha, beta]
    G[0, 1, 1] = a * ad                                          # Gamma^t_{chi chi}
    G[1, 0, 1] = G[1, 1, 0] = ad / a                             # Gamma^chi_{t chi}
    # derivative of u along the worldline
    ut = sh * bp                                                 # d/dtau cosh(beta)
    uchi = bp * ch / a - sh / a ** 2 * ad * ch                   # d/dtau (sinh(beta)/a)
    A = np.array([ut, uchi]) + np.einsum("mab,a,b->m", G, u, u)
    K_tt = -n_dn @ A                                             # K_{tau tau} = -n_nu A^nu
    Kup_tt = -K_tt                                               # K^tau_tau = g^{tau tau} K_tau tau
    s, c = np.sin(chi), np.cos(chi)
    Gth = np.array([a * ad * s * s, -s * c])                     # Gamma^mu_{theta theta}, mu = (t, chi)
    K_thth = -(Gth @ n_dn)
    R = a * s
    h_tt = u @ g @ u
    return dict(K_tau_tau=Kup_tt, K_th_th=K_thth / R ** 2, n=n_up, u=u, h_tt=h_tt, h_thth=a * a * s * s, n_norm=n_up @ g @ n_up, n_dot_u=n_up @ g @ u, Gamma=n_up[1] * 0 + (n_up @ np.array([ad * s, a * c]) * 1.0))


def raw_K_outside(R, Rd, Rdd, m, q):
    """Same from the exterior side, coordinates (s, R)."""
    f = 1.0 - m / R + q * R * R
    fp = m / R ** 2 + 2.0 * q * R
    F = f + Rd * Rd
    sd = np.sqrt(F) / f
    u = np.array([sd, Rd])
    n_up = np.array([-Rd / f, -np.sqrt(F)])                      # points toward decreasing R
    g = np.diag([-f, 1.0 / f])
    n_dn = g @ n_up
    sdd = (2 * Rd * Rdd + fp * Rd) / (2 * np.sqrt(F) * f) - np.sqrt(F) * fp * Rd / f ** 2
    G = np.zeros((2, 2, 2))
    G[0, 0, 1] = G[0, 1, 0] = fp / (2 * f)
    G[1, 0, 0] = f * fp / 2
    G[1, 1, 1] = -fp / (2 * f)
    A = np.array([sdd, Rdd]) + np.einsum("mab,a,b->m", G, u, u)
    K_tt = -n_dn @ A
    Gth = np.array([0.0, -f * R])
    K_thth = -(Gth @ n_dn)
    return dict(K_tau_tau=-K_tt, K_th_th=K_thth / R ** 2, n=n_up, u=u, h_tt=u @ g @ u, h_thth=R * R, n_norm=n_up @ g @ n_up, n_dot_u=n_up @ g @ u, s_dot=sd, f=f)


if __name__ == "__main__":
    I, yi0, T = setup()
    sh = json.load(open(os.path.join(HERE, "phase2d_shells_q0.json")))["shells"][0]
    q, w, m, chi0 = sh["q"], sh["w"], sh["m_ext"], sh["chi0_over_pi"] * np.pi
    ev = lambda tau, y: y[0] - T
    ev.terminal = True
    sol = solve_ivp(MS.rhs_factory(I, m, q, w), (0.0, 4.0 * T), [0.0, chi0, 0.0] + yi0, method="DOP853", rtol=1e-13, atol=1e-15, max_step=0.001, dense_output=True, events=[ev])
    tau = np.linspace(0, sol.t_events[0][0], 3001)
    Y = sol.sol(tau)
    rh = m                                                         # q = 0
    res = dict(h_tt_in=0.0, h_tt_out=0.0, h_th_match=0.0, n_unit=0.0, n_orth=0.0, K_tt_in=0.0, K_tt_out=0.0, sigma=0.0, p=0.0, p_minus_w_sigma=0.0)
    flags = dict(n_chi_pos=True, n_R_neg=True, t_dot_pos=True, s_dot_pos=True, outside_horizon=True)
    for i in range(len(tau)):
        y = Y[:, i]
        st = MS.shell_state(I, y, m, q, w)
        t_, chi, beta, a, ad = y[:5]
        add = st["add"]
        Kin = raw_K_inside(a, ad, add, chi, beta, st["bp"])
        Kout = raw_K_outside(st["R"], st["Rp"], st["Rpp"], m, q)
        res["h_tt_in"] = max(res["h_tt_in"], abs(Kin["h_tt"] + 1)); res["h_tt_out"] = max(res["h_tt_out"], abs(Kout["h_tt"] + 1))
        res["h_th_match"] = max(res["h_th_match"], abs(Kin["h_thth"] - Kout["h_thth"]))
        res["n_unit"] = max(res["n_unit"], abs(Kin["n_norm"] - 1), abs(Kout["n_norm"] - 1)); res["n_orth"] = max(res["n_orth"], abs(Kin["n_dot_u"]), abs(Kout["n_dot_u"]))
        res["K_tt_in"] = max(res["K_tt_in"], abs(Kin["K_tau_tau"] - (st["bp"] + ad / a * np.sinh(beta))))
        res["K_tt_out"] = max(res["K_tt_out"], abs(Kout["K_tau_tau"] + (2 * st["Rpp"] + m / st["R"] ** 2 + 2 * q * st["R"]) / (2 * np.sqrt(st["F"]))))
        # junction: [X] = out - in ; S^a_b = -([K^a_b] - delta [K])
        dKtt = Kout["K_tau_tau"] - Kin["K_tau_tau"]; dKth = Kout["K_th_th"] - Kin["K_th_th"]
        dK = dKtt + 2 * dKth
        S_tt = -(dKtt - dK); S_th = -(dKth - dK)
        sigma_raw, p_raw = -S_tt, S_th
        res["sigma"] = max(res["sigma"], abs(sigma_raw - st["sigma"])); res["p"] = max(res["p"], abs(p_raw - st["p"]))
        res["p_minus_w_sigma"] = max(res["p_minus_w_sigma"], abs(p_raw - w * sigma_raw))
        flags["n_chi_pos"] &= bool(Kin["n"][1] > 0); flags["n_R_neg"] &= bool(Kout["n"][1] < 0)
        flags["t_dot_pos"] &= bool(Kin["u"][0] > 0); flags["s_dot_pos"] &= bool(Kout["s_dot"] > 0)
        flags["outside_horizon"] &= bool(st["R"] > rh)
    out = dict(residuals={k: float(v) for k, v in res.items()}, flags=flags, m=m, r_h=rh, w=w, chi0_over_pi=sh["chi0_over_pi"])
    # curvature invariants (q = 0): exterior Kretschmann 48 m^2 / R^6; interior FRW Kretschmann 12[(add/a)^2 + (ad/a)^4 ... ] -> 12 ((add/a)^2 + (ad^2+1)^2/a^4)
    R_s = np.array([MS.shell_state(I, Y[:, i], m, q, w)["R"] for i in range(len(tau))])
    out["kretschmann_exterior_at_shell"] = [float(48 * m * m / R_s.max() ** 6), float(48 * m * m / R_s.min() ** 6)]
    out["kretschmann_exterior_at_throat"] = float(48 * m * m / rh ** 6)
    ad_ = Y[4]; a_ = Y[3]; add_ = np.array([I.rhs(0, list(Y[3:, i]))[1] for i in range(len(tau))])
    kin = 12 * ((add_ / a_) ** 2 + ((ad_ ** 2 + 1) / a_ ** 2) ** 2)
    out["kretschmann_interior"] = [float(kin.min()), float(kin.max())]
    # energy exchanged between shell and bulk (dust + scalar) per cycle, relative to the cap's energy
    S = [MS.shell_state(I, Y[:, i], m, q, w) for i in range(len(tau))]
    sig = np.array([s_["sigma"] for s_ in S]); R = np.array([s_["R"] for s_ in S]); rho = np.array([s_["rho"] for s_ in S]); P = np.array([s_["P"] for s_ in S])
    flux = -(R ** 2) * (rho + P) * np.cosh(Y[2]) * np.sinh(Y[2])       # d/dtau-type exchange term of (sigma R^2)' + p (R^2)'
    Es = 4 * np.pi * sig * R ** 2
    exch = np.trapezoid(np.abs(4 * np.pi * flux), tau)
    # cap energy: rho * volume of the cap (Vol = 4 pi a^3 (chi0/2 - sin(2 chi0)/4)), mean over the cycle
    vol = 4 * np.pi * Y[3] ** 3 * (chi0 / 2 - np.sin(2 * chi0) / 4)
    out["shell_energy_range"] = [float(Es.min()), float(Es.max())]
    out["exchanged_energy_per_cycle_abs"] = float(exch)
    out["cap_energy_mean"] = float(np.mean(rho * vol))
    out["exchange_over_cap_energy"] = float(exch / np.mean(rho * vol))
    out["exchange_over_shell_energy"] = float(exch / np.mean(Es))
    print("residuals:", {k: "%.1e" % v for k, v in res.items()})
    print("orientation flags:", flags)
    print("Kretschmann: exterior at shell %s, at throat %.2e, interior %s" % (out["kretschmann_exterior_at_shell"], out["kretschmann_exterior_at_throat"], out["kretschmann_interior"]))
    print("energy: shell %s, exchanged per cycle %.3e, cap %.3e, ratio to cap %.3e, to shell %.3e" % (out["shell_energy_range"], exch, out["cap_energy_mean"], out["exchange_over_cap_energy"], out["exchange_over_shell_energy"]))
    json.dump(out, open(os.path.join(HERE, "phase2e_junction_global.json"), "w"), indent=1)
