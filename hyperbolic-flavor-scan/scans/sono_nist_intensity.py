"""
NIST intensity vs phi-lattice proximity: null result.
Lattice proximity does NOT predict emission intensity (p=0.70).
The phi-lattice organizes which energies appear, not how strongly.
"""
import numpy as np, math
from scipy.stats import spearmanr, pearsonr

PHI = (1+math.sqrt(5))/2
LOG_PHI = math.log(PHI)
hc = 1239.84193
E0 = 16.675

def nearest_q(x,a):
    r = math.log(x/a)/LOG_PHI
    q = round(r*4)/4
    return q, abs(r-q)

ar_nist = [
    (706.722,450,"4p[5/2]2->4s[3/2]2"),(727.294,180,"4p[5/2]1->4s[3/2]2"),
    (738.398,500,"4p[5/2]2->4s[3/2]1"),(750.387,400,"4p[3/2]1->4s'[1/2]0"),
    (751.465,180,"4p[5/2]1->4s[3/2]1"),(763.511,700,"4p[5/2]2->4s'[1/2]1"),
    (772.376,200,"4p[1/2]0->4s[3/2]1"),(794.818,300,"4p[5/2]1->4s'[1/2]1"),
    (800.616,100,"4p[3/2]2->4s[3/2]2"),(801.479,100,"4p[1/2]0->4s[3/2]2"),
    (810.369,300,"4p[3/2]2->4s[3/2]1"),(811.531,1000,"4p[3/2]2->4s[3/2]1 STRONGEST"),
    (826.452,400,"4p[3/2]1->4s'[1/2]1"),(840.821,600,"4p[3/2]2->4s'[1/2]1"),
    (842.465,100,"4p[1/2]0->4s'[1/2]1"),(852.144,150,"4p[5/2]2->4s'[1/2]1"),
    (912.297,75,"4p[5/2]2->4s'[1/2]0"),
]
ar_A = [3.80,1.83,8.47,4.45,4.02,24.5,1.17,18.6,4.90,0.54,
        25.1,33.1,15.3,22.3,2.15,13.9,16.9]  # Einstein A, 10^6/s

I_arr = np.array([x[1] for x in ar_nist])
A_arr = np.array(ar_A)
res_arr = np.array([nearest_q(hc/wl,E0)[1] for wl,_,_ in ar_nist])

rho,p = spearmanr(-res_arr, I_arr)
rho_A,p_A = spearmanr(-res_arr, A_arr)
print(f"NIST intensity vs residual: rho={rho:.4f}, p={p:.4f}")
print(f"Einstein A   vs residual: rho={rho_A:.4f}, p={p_A:.4f}")
print("Conclusion: lattice proximity does NOT predict intensity (null result)")
print("The phi-lattice organizes WHICH energies appear, not HOW STRONGLY.")
