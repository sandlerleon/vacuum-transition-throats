# -*- coding: utf-8 -*-
"""Recomputes only the constraint-residual entries of validation_extra.json with the normalization scales (RMS of the largest term of each constraint)."""
import json, os
import validate_extra as VX
p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "validation_extra.json")
V = json.load(open(p))
V["constraints"] = [VX.constraint_rms(0.1, 2.0, 0.1, h) for h in (0.02, 0.01, 0.005)]
V["constraints_AT"] = [VX.constraint_rms(0.1, 4.0, 0.1, h, u_hi=0.15, v_hi=4.0) for h in (0.02, 0.01, 0.005)]
json.dump(V, open(p, "w"), indent=1)
for k in ("constraints", "constraints_AT"):
    for c in V[k]:
        print(k, c["h"], "rms_vv %.3g scale_vv %.3g rel %.3g | rms_uu %.3g scale_uu %.3g rel %.3g" % (c["rms_vv"], c["scale_vv"], c["rms_vv"] / c["scale_vv"], c["rms_uu"], c["scale_uu"], c["rms_uu"] / c["scale_uu"]))
