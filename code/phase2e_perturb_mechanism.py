# -*- coding: utf-8 -*-
"""Phase 2E Test B: attribute the inhomogeneous instability.  (1) analytic quasi-static fifth-force rate  g(t)^2 = 6 C p0^2 / (M2 a^3 (2 M2 + p0^2)) integrated over the cycle and compared with
ln(max Floquet multiplier) at large n; (2) control: perturbation equations with the dust-scalar coupling switched off in the PERTURBATIONS only (M2 -> huge): gravity-only growth."""
import json, os
import numpy as np
from scipy.integrate import solve_ivp
import phase2e_perturb_floquet as PF
from phase2d_cycle_validation import Cycle

HERE = os.path.dirname(os.path.abspath(__file__))
c = Cycle(0.35)
p = c.orbit(0.2229, span=0.05)
bg0 = np.array([c.neck_a(p, 0.0), 0.0, p, 0.0])
T = c.flow_to_bounce(p, 0.0)[0]
sol = solve_ivp(c.f, (0, T), bg0, method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True)
tt = np.linspace(0, T, 4001); Y = sol.sol(tt)
g2 = 6 * c.C * Y[2] ** 2 / (c.M2 * Y[0] ** 3 * (2 * c.M2 + Y[2] ** 2))
rate = np.sqrt(g2)
integral = np.trapezoid(rate, tt)
print("fifth-force estimate: sqrt(g2) in [%.2f, %.2f], integral over one cycle = %.2f (growth factor e^that = %.3g)" % (rate.min(), rate.max(), integral, np.exp(integral)))
out = dict(rate_min=float(rate.min()), rate_max=float(rate.max()), estimate_log_growth=float(integral), compare=[], control=[])
fl = json.load(open(os.path.join(HERE, "phase2e_perturb_floquet.json")))
for rec in fl["modes"]:
    if rec["n"] in (20, 40, 60):
        out["compare"].append(dict(n=rec["n"], ln_multiplier=float(np.log(rec["max_abs"]))))
        print("  n=%d: ln(Floquet multiplier) = %.2f" % (rec["n"], np.log(rec["max_abs"])))
# control: switch off dust-scalar coupling in the perturbation equations only
class PertNoCoupling(PF.Pert):
    def __init__(self, cyc, n):
        super().__init__(cyc, n)
        self.par = (cyc.mu, cyc.lam, cyc.C, 1e14)
for n in (3, 10, 30, 60):
    P = PertNoCoupling(c, n)
    M, _ = P.monodromy(bg0, T)
    ev = np.linalg.eigvals(M)
    out["control"].append(dict(n=n, multipliers_abs=sorted([float(abs(e)) for e in ev], reverse=True)))
    print("  control (perturbation coupling off) n=%d: |multipliers| %s" % (n, ", ".join("%.5g" % abs(e) for e in sorted(ev, key=lambda z: -abs(z)))))
json.dump(out, open(os.path.join(HERE, "phase2e_perturb_mechanism.json"), "w"), indent=1)
