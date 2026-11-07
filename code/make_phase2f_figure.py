# -*- coding: utf-8 -*-
"""FigF1: Phase 2F stabilization attempts."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
L = lambda f: json.load(open(os.path.join(HERE, f)))
PF = L("phase2e_perturb_floquet.json"); WC = L("phase2f_worstcase.json"); PR = L("phase2f_pressure.json"); CU = L("phase2f_continuation_up_mu20.json"); MU = L("phase2f_muscan_up_rc0.2_valid.json"); PO = L("phase2f_potential.json")
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
a = ax[0]
ns = [r["n"] for r in PF["modes"]]; a.loglog(ns, [r["max_abs"] for r in PF["modes"]], "k-o", ms=3, label="baseline (mu=20, rho_crit=0.2 V_f)")
for r in WC[1:]:
    k = sorted(int(x) for x in r["growth"]); a.loglog(k, [r["growth"][str(x)] for x in k], "-", label="mu = %.0f" % r["mu"])
for r in PR[2:5]:
    k = sorted(int(x) for x in r["growth"]); a.loglog(k, [r["growth"][str(x)] for x in k], "--", label="pressure c_s^2 = %.2f" % r["cs2_neck"])
a.axhline(1, color="gray", lw=0.7); a.set_xlabel("radial harmonic n"); a.set_ylabel("Floquet growth per cycle G(n)"); a.set_title("(a) Growth versus n", fontsize=10); a.legend(fontsize=6)
b = ax[1]
rc = [r["rc"] for r in CU]; b.semilogx(rc, [r["growth"]["3"] for r in CU], "o-", label="n = 3"); b.semilogx(rc, [r["growth"]["60"] / 1e4 for r in CU], "s-", label="n = 60 (/10^4)")
b.set_xlabel("coupling strength rho_crit / V_f (fixed return-map trace)"); b.set_ylabel("growth per cycle"); b.set_title("(b) Weakening the dust-scalar coupling", fontsize=10); b.legend(fontsize=7)
c = ax[2]
mu = [r["mu"] for r in MU]; c.loglog(mu, [r["growth"]["3"] for r in MU], "o-", label="n = 3"); c.loglog(mu, [r["growth"]["12"] for r in MU], "s-", label="n = 12")
c.loglog([r["mu"] for r in WC], [r["worst_G"] for r in WC], "k^", label="worst case over all n")
c.set_xlabel("potential scale mu (fixed return-map trace)"); c.set_ylabel("growth per cycle"); c.set_title("(c) Changing the hilltop", fontsize=10); c.legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(HERE, "..", "figures", "FigF1_phase2f.png"), dpi=150); print("saved")
