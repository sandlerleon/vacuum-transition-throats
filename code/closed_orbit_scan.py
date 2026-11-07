# -*- coding: utf-8 -*-
"""Scan of the parameter space for exact time-symmetric cycles (closed_orbits.py) that complete at least two bounces in double precision, and their excursion a_T / a_neck."""
import json, os, sys
import numpy as np
from scipy.optimize import brentq
import closed_cycle as CC
import closed_orbits as CO

HERE = os.path.dirname(os.path.abspath(__file__))
rows = []
for q in (1.0,):
    for mu in (10.0, 20.0):
        for rc in (0.05, 0.2, 0.4):
            for x in (0.3, 0.6, 0.9):
                m = CC.Model(q=q, rc=rc, mu=mu, C=x * CC.C_CRIT)
                roots, grid, vals = CO.find_orbits(m, lo=-4.5, hi=-0.3, n=70)
                good = []
                for r in roots:
                    o = CC.run(m, seed=10 ** r * m.v, t_max=80.0, max_turns=6, rtol=1e-12)
                    seq = "".join("B" if t["kind"] == "bounce" else "T" for t in o["turns"])
                    a_neck = o["a_neck_initial"]
                    aT = max([t["a"] for t in o["turns"] if t["kind"] == "turnaround"] or [0])
                    nB = seq.count("B")
                    if nB >= 1:
                        good.append(dict(log10_seed=float(r), sequence=seq, bounces=nB, a_neck=a_neck, a_turn_max=aT, excursion=aT / a_neck, period=o["turns"][1]["t"] if len(o["turns"]) > 1 else None))
                rows.append(dict(q=q, mu=mu, rc=rc, C_over_Ccrit=x, n_roots=len(roots), closing=good))
                best = max(good, key=lambda g: g["excursion"]) if good else None
                print("mu=%g rc=%g C/Cc=%.1f: %d roots, %d close at least one bounce; max excursion %s seq %s" % (mu, rc, x, len(roots), len(good), None if not best else round(best["excursion"], 2), None if not best else best["sequence"]), flush=True)
json.dump(rows, open(os.path.join(HERE, "closed_orbit_scan.json"), "w"), indent=1)
