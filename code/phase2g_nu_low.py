# one nu per process (timeout in the caller): prints "nu max min"
import sys, numpy as np
import phase2f_model as F
nu = float(sys.argv[1])
m = F.Model2(0.35); cd = F.cycle_data(m, 0.2229)
P = F.Pert2(m, 3); P.k2 = nu * nu - 1.0
M = P.monodromy(cd["bg0"], cd["T"], rtol=1e-9)
ev = np.linalg.eigvals(M)
print("%.4f %.6g %.6g" % (nu, max(abs(ev)), min(abs(ev))))
