# -*- coding: utf-8 -*-
"""Patch phase2e_paper_content.py: new Section 5.5 (finite cap and boundary conditions), a sentence in 5.4, a row in Table 5."""
p = "phase2e_paper_content.py"
s = open(p, encoding="utf-8").read()
if "5.5 Boundary conditions" in s:
    raise SystemExit("already patched")
s = s.replace('FIG3 = "FigF1_phase2f.png"',
              'FIG3 = "FigF1_phase2f.png"\nCAP = L("phase2g_cap_spectrum.json")\nSHG = L("phase2g_shell_data.json")\nNEC0 = L("phase2g_comoving_nec.json")')
sec = r'''add(("h2", "5.5 Boundary conditions on the finite cap"))
add(("p", "Section 5.2 used the radial harmonics of the full three-sphere, whereas the bridge contains a cap bounded by the shell, and the scalar field must obey some condition there. For the homogeneous interior $\partial_n\phi=\sinh\beta\,\dot\phi$, so a no-flux condition holds only for a comoving shell ($\beta=0$); "
          "the moving shells of Table 2 need an absorbing condition for which no field theory is given. For a comoving shell and a field confined to the cap, zero scalar flux and zero normal dust velocity reduce to $Q_\nu'(\chi_0)=0$ for the modes $Q_\nu=\sin(\nu\chi)/\sin\chi$, that is $\tan(\nu\chi_0)=\nu\tan\chi_0$ ($\nu=1$ is the homogeneous mode). "
          "Comoving shells with $\sigma>0$ and $\sigma+p\ge0$ over the cycle exist for many exteriors (%d of 9900 combinations), but for the anti-de Sitter exterior that matches the potential, $V(\pm v)=-V_f$, only for $\chi_0\lesssim0.6\pi$, with $\sigma\approx%.1f$–%.1f and $p/\sigma$ varying between %.2f and %.2f. "
          "On such caps the lowest mode grows by %.1f–%.1f per cycle (at $\chi_0=0.6\pi$ and $0.55\pi$) and the worst case is unchanged at %s. The growth $G(\nu)$ exceeds one for every $\nu>1$, tends to one only at $\nu=1$, and rises monotonically for $\nu>2$; a boundary condition only quantizes $\nu$, so with infinitely many modes at $\nu\gtrsim3$ (where $G>14$) none can stabilize the cap. "
          "A shell with scalar-dependent tension $\sigma(\phi)$ would give $p=-\sigma$ and a Robin condition, but continuity of $\phi$ across the shell conflicts with a static exterior vacuum: the field would have to equal $v$ at the shell while the interior field varies between %.2f$v$ and %.2f$v$, and the tension of a thin wall between them, %.2f–%.2f, is far below the $\sigma\approx2.3$ the shell requires. We did not pursue it." %
          (NEC0["n_nec_full"], min(r["sigma_range"][0] for r in SHG["shells"]), max(r["sigma_range"][1] for r in SHG["shells"]), min(r["w_range"][0] for r in SHG["shells"]), max(r["w_range"][1] for r in SHG["shells"]),
           [r for r in CAP if abs(r["chi0_over_pi"] - 0.6) < 1e-9][0]["growth"][0], [r for r in CAP if abs(r["chi0_over_pi"] - 0.55) < 1e-9][0]["growth"][0], "7.7×10⁴",
           SHG["phi_range_over_v"][0], SHG["phi_range_over_v"][1], min(SHG["wall_tension"]), max(SHG["wall_tension"]))))
'''
sec = sec.replace('%d of 9900', '%d of 9900')
sec = sec.replace(chr(92), chr(92) * 2)
marker = 'add(("h1", "6. Discussion"))'
s = s.replace(marker, sec + marker)
s = s.replace("the energy bookkeeping is satisfied only through an unmodelled interaction.", "the energy bookkeeping is satisfied only through an unmodelled interaction. For a comoving shell with the field confined to the cap the exchange vanishes (Section 5.5).")
bs = chr(92)
old = '"fails (linear, $n' + bs + bs + 'ge3$)"]'
new = '"fails (linear, $n' + bs + bs + 'ge3$); not removed by weaker coupling, other hilltops, pressure or the boundary conditions of the finite cap"]'
assert old in s
s = s.replace(old, new)
open(p, "w", encoding="utf-8").write(s)
print("patched")
