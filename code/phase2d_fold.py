# -*- coding: utf-8 -*-
"""Continue the main orbit branch in phi0 (instead of C) past the apparent end at C/Cc ~ 0.3602: a fold (saddle-node) shows as a maximum of C(phi0)."""
import json, numpy as np
from scipy.optimize import brentq
from phase2d_cycle_validation import Cycle
out=[]
xprev=0.36
for z in np.arange(0.28,0.40,0.004):
    f=lambda x: Cycle(x).g_turn(z*Cycle(x).v)
    lo,hi=xprev-0.004,xprev+0.004
    try:
        flo,fhi=f(lo),f(hi)
        if flo is None or fhi is None or flo*fhi>0:
            lo,hi=xprev-0.012,xprev+0.012
            flo,fhi=f(lo),f(hi)
        if flo is None or fhi is None or flo*fhi>0:
            print("phi0/v=%.3f: no root near C/Cc=%.4f"%(z,xprev)); break
        x=brentq(f,lo,hi,xtol=1e-12)
    except Exception as e:
        print("phi0/v=%.3f: %s"%(z,e)); break
    c=Cycle(x); p=z*c.v; J=c.jac_neck(p,0.0)
    out.append(dict(phi0_over_v=float(z),x=float(x),trace=float(np.trace(J)),det=float(np.linalg.det(J)),period=float(c.flow_to_bounce(p,0.0)[0])))
    print("phi0/v=%.3f  C/Cc=%.6f  trace=%+.5f  det=%.6f  T=%.5f"%(z,x,np.trace(J),np.linalg.det(J),out[-1]['period']),flush=True)
    xprev=x
json.dump(out,open('phase2d_fold.json','w'),indent=1)
