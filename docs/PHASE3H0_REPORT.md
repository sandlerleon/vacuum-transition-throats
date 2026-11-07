# Phase 3H.0: causal shear viscosity relaxes the anisotropy threshold of the spin-torsion bounce (branch `phase3-recycling`)

Reproduce: `code/phase3h0_viscous_spin_bounce.py` -> `code/phase3h0_viscous_spin_bounce.json` (about 2 min); test `phase3h0_viscous_relaxation_of_anisotropy_threshold` (40 checks pass).

Combines Phase 3C.2 (axisymmetric Bianchi I, spin fluid rho_s = p_s = -B n^2, critical spin strength B_c set by the shear) with Phase 3F.0 (Israel-Stewart-type shear stress, eta = c^2 tau n, tau = K/max(|Theta|, 0.01)).

Result: the inviscid B_c reproduces the analytic anisotropy threshold (0.0025 for initial shear s0 = 0.05; 0.04 for s0 = 0.2). With causal viscosity (c^2 <= 1, Knudsen K <= 1, viscous stress below 0.45 of the dust energy density, constraint residual below 1e-4) B_c falls monotonically with K: factor 2.7 (K = 0.1), 12 (K = 0.3), 730 (K = 1, c^2 = 1) for s0 = 0.05; 1.8, 4, 19 for s0 = 0.2.

Approximations: dust heated by viscous work treated as pressureless; viscosity acts on the dust only; truncated Israel-Stewart; homogeneous anisotropic model; K > 1 rows (K = 3) are outside the hydrodynamic regime and unreliable (constraint residual 1.5e-3). The result relaxes the anisotropy threshold; it does not remove the fluid-sector instability of the spin-fluid equation of state (Phase 3C / Section 3.1 of the manuscript).
