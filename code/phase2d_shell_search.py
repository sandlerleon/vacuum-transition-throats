# -*- coding: utf-8 -*-
"""Phase 2D: periodic moving thin shells around the stable single-field cold cycle (C/Cc = 0.35).  python phase2d_shell_search.py <q>   (writes phase2d_shells_q<q>.json)"""
import json, os, sys
import numpy as np
import bridge_moving_shell as MS
from phase2d_cycle_validation import Cycle

HERE = os.path.dirname(os.path.abspath(__file__))


class Interior:
    """Adapter giving the generic shell code the interface rhs(t, yi), rho_P(yi) for the independent Cycle model; yi = [a, a_dot, phi, phi_dot]."""
    def __init__(self, x=0.35):
        self.c = Cycle(x)
        self.rhs = self.c.f

    def rho_P(self, yi):
        a, b, p, w = yi
        return self.c.rho(a, p, w), 0.5 * w * w - self.c.V(p)


def setup(x=0.35, guess=0.2229):
    I = Interior(x)
    p = I.c.orbit(guess, span=0.05)
    a0 = I.c.neck_a(p, 0.0)
    T = I.c.flow_to_bounce(p, 0.0)[0]
    return I, [a0, 0.0, p, 0.0], T


if __name__ == "__main__":
    q = float(sys.argv[1])
    I, yi0, T = setup()
    out = []
    for m_frac in (0.1, 0.3, 0.5, 0.7, 0.85):
        for w in (-0.5, -0.25, 0.0, 1.0 / 3.0, 1.0):
            res = MS.periodic_shell_generic(I, yi0, T, (0.505, 0.97), m_frac, q, w, ngrid=100)
            for r in res:
                r.update(q=q, m_frac=m_frac, w=w)
                out.append(r)
                d = r["diag"]
                print("q=%g m/mmax=%.2f w=%+.2f chi0=%.4fpi m=%.4f | mult=%s | sigma>0 %.2f NEC %.2f | |beta|max %.2f R[%.3f,%.3f] flux %.1e" % (q, m_frac, w, r["chi0_over_pi"], r["m_ext"],
                      ",".join("%.4f" % z for z in r["multipliers_abs"]) if r["multipliers_abs"] else "-", d["sigma_pos"], d["nec"], d["max_abs_beta"], d["R_min"], d["R_max"], d["flux_identity_resid"]), flush=True)
    print("done q=%g, %d periodic shells" % (q, len(out)))
    json.dump(dict(yi0=yi0, T=T, shells=out), open(os.path.join(HERE, "phase2d_shells_q%g.json" % q), "w"), indent=1)
