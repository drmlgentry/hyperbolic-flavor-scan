"""
Test corrected encoding q = 5a + 9b + 14c (from symmetric space tangent space dimensions)
vs original q = 8a + 15b + 24c
"""
import numpy as np
from itertools import product

# PDG 2024 fermion masses in GeV
masses = {
    "e":      0.000510999,
    "mu":     0.105658375,
    "tau":    1.77686,
    "u":      0.00216,
    "d":      0.00467,
    "s":      0.0934,
    "c":      1.275,
    "b":      4.18,
    "t":      172.76,
    "W":      80.379,
    "Z":      91.1876,
    "H":      125.25,
}

phi = (1 + 5**0.5) / 2
mref = masses["e"]

def log_phi(m):
    return np.log(m / mref) / np.log(phi)

def q_original(a, b, c):
    return 8*a + 15*b + 24*c

def q_corrected(a, b, c):
    return 5*a + 9*b + 14*c

def mass_from_q(q):
    return mref * phi**(q/4)

print("Corrected encoding q = 5a + 9b + 14c")
print("="*60)

scan_range = range(-120, 121)

for name, m in masses.items():
    target = log_phi(m)
    best_res = 999
    best_abc = None
    best_q = None
    for a in range(-100, 101):
        for b in range(-40, 41):
            for c in range(-30, 31):
                q = q_corrected(a, b, c)
                mpred = mass_from_q(q)
                res = abs(log_phi(mpred) - target)
                if res < best_res:
                    best_res = res
                    best_abc = (a, b, c)
                    best_q = q
    mpred = mass_from_q(best_q)
    print(f"{name:4s}: m={m:.6g} GeV  mpred={mpred:.6g}  residual={best_res:.4f}  abc={best_abc}  q={best_q}")

print("\nOriginal encoding q = 8a + 15b + 24c")
print("="*60)
for name, m in masses.items():
    target = log_phi(m)
    best_res = 999
    best_abc = None
    best_q = None
    for a in range(-100, 101):
        for b in range(-40, 41):
            for c in range(-30, 31):
                q = q_original(a, b, c)
                mpred = mass_from_q(q)
                res = abs(log_phi(mpred) - target)
                if res < best_res:
                    best_res = res
                    best_abc = (a, b, c)
                    best_q = q
    mpred = mass_from_q(best_q)
    print(f"{name:4s}: m={m:.6g} GeV  mpred={mpred:.6g}  residual={best_res:.4f}  abc={best_abc}  q={best_q}")
