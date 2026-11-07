# -*- coding: utf-8 -*-
"""One-off patch of phase2e_paper_content.py: data loading, revised abstract, new Section 5.3 (attempts to stabilize), discussion edits."""
p = "phase2e_paper_content.py"
s = open(p, encoding="utf-8").read()
if "FIG3" in s:
    raise SystemExit("already patched")
s = s.replace('CV = L("phase2d_convergence.json")',
              'CV = L("phase2d_convergence.json")\nPR = L("phase2f_pressure.json")\nCU = L("phase2f_continuation_up_mu20.json") + L("phase2f_continuation_down_mu20.json")\nMUF = L("phase2f_muscan_up_rc0.2.json")\nPOT = L("phase2f_potential.json")\nAN = L("phase2f_analysis.json")\nFIG3 = "FigF1_phase2f.png"')
a = s.index('Linear spherical perturbations of the interior, however')
b = s.index('    % (g(mode[3]["max_abs"]')
newabs = ('Linear spherical perturbations of the interior, however, grow by factors of %s to %s per cycle, driven by a scalar-mediated force on the coupled dust and by the hilltop instability of the field. "\n'
          '    "Continuation in the coupling and the potential scale, a different hilltop shape and a fluid with pressure do not remove them: the worst-case growth stays above %s. The shell also exchanges energy with the bulk at a rate "\n'
          '    "comparable to its own energy, which needs an unspecified interaction. The construction is an exact, homogeneous, idealized solution, not a robust cosmological model; we state what is demonstrated and what is open."\n')
s = s[:a] + newabs + s[b:]
s = s.replace('% (g(mode[3]["max_abs"], 2), sci(mode[60]["max_abs"], 1).replace(" × 10$^{5}$", "×10⁵")))',
              '% (g(mode[3]["max_abs"], 2), sci(mode[60]["max_abs"], 1).replace(" × 10$^{5}$", "×10⁵"), "3×10³"))')
s = s.replace("A pressure or finite-range force could weaken the coupling-driven part, but nothing here shows that it does, and a scan of the coupling strength at a coarse step found an elliptic cycle only at the value used. ",
              "Section 5.3 tests whether changing the model removes this growth. ")
BS = chr(92)
new53 = r'''add(("h2", "5.3 Attempts to stabilize"))
add(("p", "We re-derived the perturbation equations for a polynomial potential $V=3+a_2\phi^2+a_4\phi^4+a_6\phi^6$ and for a coupled fluid with pressure, a polytrope $P=Kn^{5/3}$ with $\rho_f=An+Kn^{5/3}/(\gamma-1)$, whose motion follows from total energy–momentum conservation; at $K=a_6=0$ they reproduce Table 3. "
          "Each change was followed by continuation along the elliptic orbit at fixed return-map trace ($-1.7043$), solving for $\phi_0$ and adjusting $C/C_c$ at every step, not by a scan of initial conditions (Table 4, Fig. 3). "
          "Changing the dust–scalar coupling from $\rho_{\rm crit}/V_f=%.4f$ to %.2f changes the $n=3$ growth only between %.1f and %.1f, because the orbit retunes $C/C_c$ and the fifth-force rate at the neck stays near 16–24. Raising the hilltop curvature $\mu$ to %.0f reduces the fifth-force part ($g\,T$ falls from 10.7 to 0.5) but widens the tachyonic band to $n\lesssim\mu$, "
          "and the worst growth over all $n$ stays above $3.6\times10^{3}$. A sextic reshaping of the hilltop at fixed curvature and well depth changes the $n=3$ growth by about 20%%. Pressure removes only modes with $n\gtrsim n_c=ga/c_s$: at $c_s^2=0.04$ the modes $n=240$ are stable and $n=120$ is reduced to %.0f, while $n=3$–12 are unchanged, "
          "and the cycle ceases to exist before $c_s^2$ reaches 0.08." % (min(r["rc"] for r in CU), max(r["rc"] for r in CU), min(r["growth"]["3"] for r in CU), max(r["growth"]["3"] for r in CU), MUF[-1]["mu"], [r for r in PR if r["K"] == 0.08][0]["growth"]["120"])))
RNG = lambda rows, k: "%.3g–%.3g" % (min(r["growth"][k] for r in rows), max(r["growth"][k] for r in rows))
POTALL = POT["wider"] + POT["narrower"]
add(("table", "Table 4. Inhomogeneous growth per cycle under each attempt (range of the largest Floquet multiplier along the continuation).",
     ["change (range followed)", "$n=3$", "$n=12$", "$n=60$"],
     [["baseline", g(mode[3]["max_abs"], 3), g(mode[12]["max_abs"], 3), g(mode[60]["max_abs"], 3)],
      ["coupling $\rho_{\rm crit}/V_f$ from 0.0086 to 0.51", RNG(CU, "3"), RNG(CU, "12"), RNG(CU, "60")],
      ["hilltop scale $\mu$ from 20 to 307", RNG(MUF, "3"), RNG(MUF, "12"), RNG(MUF, "60")],
      ["sextic hilltop, $v/v_0$ from 0.98 to 1.74", RNG(POTALL, "3"), RNG(POTALL, "12"), RNG(POTALL, "60")],
      ["pressure, $K=0$ to 0.08 ($c_s^2\le0.04$)", RNG(PR, "3"), RNG(PR, "12"), RNG(PR, "60")]],
     [3.0, 1.3, 1.3, 1.3]))
WK = [w for w in AN["wkb"] if w["family"] == "pressure" and w["n"] in (120, 240) and w["K"] in (0.04, 0.08)]
add(("p", "Two criteria summarize why. In the quasi-static limit the fluid contrast obeys $\ddot D=(g^2-c_s^2k^2/a^2)D$, so mode $n$ grows when $n<n_c=ga/c_s$; with pressure the exponent $\int\sqrt{g^2-c_s^2k^2/a^2}\,dt$ equals the Floquet exponent to within %.2f for $n=120$ and 240 at $K=0.04$ and 0.08, and fails at low $n$ where the quasi-static assumption breaks. "
          "The scalar sector alone is tachyonic, $\ddot{\delta\phi}+3H\dot{\delta\phi}+(k^2/a^2+V''+n/M^2)\delta\phi=0$ with $V''<0$ near the top, for $k<a\mu$; its scalar-only multiplier accounts for most of the growth at large $\mu$ (7.1 of 7.9 in the logarithm at $n=60$, $\mu=208$). "
          "Along the whole family $\mu T$ stays between 12 and 16, so the hilltop exponent does not shrink when $\mu$ grows. A stable cycle would need both the fifth-force and the hilltop band closed; pressure closes only the first and only above $n_c$, and nothing tried closes the second. "
          "These are mechanism-level estimates validated numerically, not rigorous bounds, and the low-$n$ exponent ($\ln G\approx2.2$–$2.7$) is an observed invariant, not a derived one. We did not run nonlinear evolutions because no promising region was found." % max(abs(w["S5"] - w["lnG"]) for w in WK)))
add(("fig", FIG3, "Figure 3. Attempts to stabilize. (a) Growth per cycle versus $n$ for the baseline, for larger hilltop curvature $\mu$, and for a fluid with pressure. (b) Weakening the dust–scalar coupling at fixed return-map trace. (c) Changing the hilltop curvature: growth at $n=3$ and 12 and the worst case over all $n$.", 6.4))
'''
new53 = new53.replace(BS, BS + BS)
marker = 'add(("h2", "5.3 The shell matter"))'
s = s.replace(marker, new53 + 'add(("h2", "5.4 The shell matter"))')
s = s.replace("the scan of the parameter space here is far from complete, and narrow windows with weaker coupling cannot be excluded. ",
              "continuation in the coupling, the potential scale, the hilltop shape and a fluid with pressure (Section 5.3) did not remove the growth, though other orbits, windows and matter models were not explored. ")
s = s.replace("The coupling-driven instability depends on the dust–scalar coupling that the model needs for its transition, so a different mechanism for the field's return to the false vacuum may be required; ",
              "The instability has two sources, the fifth force on the coupled fluid and the hilltop of the potential, so a different mechanism for the field's return to the false vacuum may be required; ")
open(p, "w", encoding="utf-8").write(s)
print("patched")
