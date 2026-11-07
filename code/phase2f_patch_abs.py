import re
p = "phase2e_paper_content.py"
s = open(p, encoding="utf-8").read()
a = s.index("ABSTRACT = (")
b = s.index('add(("p0", ABSTRACT))')
new = '''ABSTRACT = (
    "We construct and test a candidate spacetime in which a closed, nonsingular, cyclically evolving cosmological region is joined across a moving thin shell to a Schwarzschild exterior, "
    "so that a spatial slice has two ends and a throat. The region is the closed Friedmann universe of a scalar field with a double-well potential coupled to dust; for a window of dust content "
    "it has exact bouncing cycles with an elliptic return map, reproduced with an independent code. For $w=-1/2$ and a flat exterior the shell has periodic solutions with positive tension that satisfy "
    "the null energy condition, verified from raw extrinsic curvatures to $10^{-15}$. The cyclic region and the parent universe are spacelike separated: the bridge is not traversable and nothing crosses it. "
    "The data need a pre-existing throat but no white-hole region to the future of a suitable slice. The conservative shell survives 60 cycles under 1%% perturbations. However, linear spherical perturbations of the interior "
    "grow by factors of %s to %s per cycle, through a scalar-mediated force on the coupled dust and the hilltop instability; changing the coupling, the potential scale or shape, or adding fluid pressure does not remove them "
    "(worst case above 3×10³). The shell's energy exchange with the bulk needs an unspecified interaction. The construction is an exact, idealized, homogeneous solution, not a robust cosmological model; "
    "we state what is demonstrated and what is open."
    % (g(mode[3]["max_abs"], 2), sci(mode[60]["max_abs"], 1).replace(" × 10$^{5}$", "×10⁵")))
'''
s = s[:a] + new + s[b:]
open(p, "w", encoding="utf-8").write(s)
