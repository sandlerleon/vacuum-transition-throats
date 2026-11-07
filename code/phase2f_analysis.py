# -*- coding: utf-8 -*-
"""Phase 2F items 6 and 8: (1) WKB criteria for the two instability mechanisms, validated against the Floquet results; (2) dust/scalar decomposition of the shell-bulk energy exchange.
Criteria (quasi-static / WKB): fifth force   S_5(n) = int sqrt(max(g^2 - c_s^2 k^2/a^2, 0)) dt ,  g^2 = 6 C phi0^2 / (M^2 a^3 (2 M^2 + phi0^2));
hilltop     S_h(n) = int sqrt(max(-(k^2/a^2 + V'' + rho_*/M^2), 0)) dt ;  ln G(n) ~ S_5 + S_h (upper-bound-like; sums are compared with the numerics)."""
import json, os
import numpy as np
from scipy.integrate import solve_ivp
import phase2f_model as F
import bridge_moving_shell as MS
from phase2d_shell_search import Interior, setup

HERE = os.path.dirname(os.path.abspath(__file__))
res = dict(wkb=[], shell_split=None)


def S_of(m, cd, n, include_fifth=True):
    k2 = n * n - 1.0
    sol = solve_ivp(m.f, (0, cd["T"]), cd["bg0"], method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True)
    tt = np.linspace(0, cd["T"], 6001); Y = sol.sol(tt)
    a, p = Y[0], Y[2]
    nn = 3 * m.C / a ** 3
    Vpp = 2 * m.a2 + 12 * m.a4 * p ** 2 + 30 * m.a6 * p ** 4
    A = 1 + p ** 2 / (2 * m.M2)
    cs2 = (m.gam * m.K * nn ** (m.gam - 1) / (A + m.gam * m.K * nn ** (m.gam - 1) / (m.gam - 1))) if m.K else 0 * a
    g2 = 6 * m.C * p ** 2 / (m.M2 * a ** 3 * (2 * m.M2 + p ** 2))
    S5 = np.trapezoid(np.sqrt(np.maximum(g2 - cs2 * k2 / a ** 2, 0)), tt)
    Sh = np.trapezoid(np.sqrt(np.maximum(-(k2 / a ** 2 + Vpp + nn / m.M2), 0)), tt)
    return float(S5), float(Sh)


# (a) pressure family
P = json.load(open(os.path.join(HERE, "phase2f_pressure.json")))
for r in P:
    m = F.Model2(r["x"], mu=20.0, rc=0.2, K=r["K"], gam=r["gam"])
    cd = F.cycle_data(m, r["phi0_over_v"], span=0.04)
    for n in (3, 30, 60, 120, 240):
        S5, Sh = S_of(m, cd, n)
        res["wkb"].append(dict(family="pressure", K=r["K"], n=n, S5=S5, Sh=Sh, lnG=float(np.log(r["growth"][str(n)] if str(n) in r["growth"] else np.nan))))
        print("K=%-5g n=%3d: S5=%.2f Sh=%.2f  S5+Sh=%.2f  vs ln G = %.2f" % (r["K"], n, S5, Sh, S5 + Sh, res["wkb"][-1]["lnG"]), flush=True)
# (b) mu family (hilltop dominated)
fam = json.load(open(os.path.join(HERE, "phase2f_muscan_up_rc0.2.json")))
for t in (20, 80, 208):
    r = min(fam, key=lambda q: abs(q["mu"] - t))
    m = F.Model2(r["x"], mu=r["mu"], rc=0.2)
    cd = F.cycle_data(m, r["phi0_over_v"], span=0.04)
    for n in (3, 12, 30, 60):
        S5, Sh = S_of(m, cd, n)
        G = F.growth(m, cd, (n,))[n]
        res["wkb"].append(dict(family="mu", mu=r["mu"], n=n, S5=S5, Sh=Sh, lnG=float(np.log(G))))
        print("mu=%-6.1f n=%3d: S5=%.2f Sh=%.2f  S5+Sh=%.2f  vs ln G = %.2f" % (r["mu"], n, S5, Sh, S5 + Sh, np.log(G)), flush=True)
# (c) shell energy exchange decomposition (conservative periodic shell of Phase 2D)
I, yi0, T = setup()
sh = json.load(open(os.path.join(HERE, "phase2d_shells_q0.json")))["shells"][0]
q, w, m_, chi0 = sh["q"], sh["w"], sh["m_ext"], sh["chi0_over_pi"] * np.pi
from scipy.integrate import solve_ivp as sivp
ev = lambda tau, y: y[0] - T; ev.terminal = True
sol = sivp(MS.rhs_factory(I, m_, q, w), (0.0, 4 * T), [0.0, chi0, 0.0] + yi0, method="DOP853", rtol=1e-12, atol=1e-14, max_step=0.001, dense_output=True, events=[ev])
tau = np.linspace(0, sol.t_events[0][0], 4001); Y = sol.sol(tau)
S = [MS.shell_state(I, Y[:, i], m_, q, w) for i in range(len(tau))]
R = np.array([s["R"] for s in S]); rho = np.array([s["rho"] for s in S]); P_ = np.array([s["P"] for s in S])
phid = Y[6]; cs = np.cosh(Y[2]) * np.sinh(Y[2])
rhod = np.array([I.c.rho(Y[3, i], Y[5, i], 0.0) for i in range(len(tau))]) - np.array([I.c.V(Y[5, i]) for i in range(len(tau))])   # dust part: rho - V  at phid = 0 -> A rho_*
flux_dust = 4 * np.pi * R ** 2 * rhod * cs; flux_scal = 4 * np.pi * R ** 2 * phid ** 2 * cs
ed = float(np.trapezoid(np.abs(flux_dust), tau)); es = float(np.trapezoid(np.abs(flux_scal), tau))
res["shell_split"] = dict(dust_gross=ed, scalar_gross=es, dust_fraction=ed / (ed + es), total_gross=ed + es,
                          net_dust=float(np.trapezoid(flux_dust, tau)), net_scalar=float(np.trapezoid(flux_scal, tau)))
print("shell exchange per cycle: dust %.3f, scalar %.3f (dust fraction %.2f); net dust %.2e, net scalar %.2e" % (ed, es, ed / (ed + es), res["shell_split"]["net_dust"], res["shell_split"]["net_scalar"]))
json.dump(res, open(os.path.join(HERE, "phase2f_analysis.json"), "w"), indent=1)
