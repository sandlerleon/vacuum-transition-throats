import json, numpy as np
from phase2d_cycle_validation import Cycle
out=[]
g_prev=[(0.3575,0.26508),(0.36,0.28182)]
x=0.36
for k in range(40):
    x=round(x+0.0005,5)
    g=g_prev[-1][1]+(g_prev[-1][1]-g_prev[-2][1])/(g_prev[-1][0]-g_prev[-2][0])*(x-g_prev[-1][0])
    c=Cycle(x); p=None
    for span in (0.02,0.05,0.12):
        p=c.orbit(g,span=span)
        if p is not None: break
    if p is None:
        print("x=%.4f: orbit lost (guess %.4f)"%(x,g)); out.append(dict(x=x,found=False)); break
    J=c.jac_neck(p,0.0); a0=c.neck_a(p,0.0); add=c.f(0,[a0,0,p,0])[1]
    tr=float(np.trace(J)); g_prev.append((x,p/c.v))
    out.append(dict(x=x,found=True,phi0_over_v=float(p/c.v),trace=tr,det=float(np.linalg.det(J)),a_neck=float(a0),add_at_neck=float(add),period=float(c.flow_to_bounce(p,0.0)[0])))
    print("x=%.4f phi0/v=%.5f trace=%+.5f det=%.6f T=%.5f a_neck=%.5f a_ddot(neck)=%+.4f"%(x,p/c.v,tr,np.linalg.det(J),out[-1]['period'],a0,add),flush=True)
json.dump(out,open('phase2d_window_end.json','w'),indent=1)
