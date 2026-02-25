import math
from math import gcd
import matplotlib.pyplot as plt

# ---------------- CONFIG ----------------
N_VALUES = [60, 100, 120, 180]
WEIGHTS = {
    "8-15-24": (8,15,24),
    "7-15-24": (7,15,24),
    "8-14-25": (8,14,25),
    "9-15-23": (9,15,23),
}

TARGETS = {
    "CKM": (65.718, 1.5),
    "PMNS_NO": (177.0, 10.0),
    "PMNS_IO": (285.0, 10.0),
}

# ---------------- FUNCTIONS ----------------
def circ_dist(a,b):
    d = abs(a-b) % 360
    return min(d, 360-d)

def residues(N, L, w):
    r=set()
    w1,w2,w3=w
    amax = L//w1
    bmax = L//w2
    cmax = L//w3
    for a in range(-amax, amax+1):
        for b in range(-bmax, bmax+1):
            for c in range(-cmax, cmax+1):
                if a==b==c==0:
                    continue
                if gcd(gcd(abs(a),abs(b)),abs(c)) != 1:
                    continue
                if w1*abs(a)+w2*abs(b)+w3*abs(c) <= L:
                    r.add((w1*a+w2*b+w3*c) % N)
    return r

def lambda_sat(N,w, Lmax=500):
    for L in range(1, Lmax+1):
        if len(residues(N,L,w)) == N:
            return L
    raise RuntimeError(f"No saturation found for N={N}, weights={w} up to Lmax={Lmax}")

def lambda_adm(N,w,delta,tol, Lmax=500):
    for L in range(1, Lmax+1):
        res = residues(N,L,w)
        for k in res:
            if circ_dist(360*k/N, delta) <= tol:
                return L
    raise RuntimeError(f"No admission found for N={N}, weights={w}, delta={delta} up to Lmax={Lmax}")

# ---------------- TABLE ----------------
rows=[]
for name,w in WEIGHTS.items():
    for N in N_VALUES:
        Lsat = lambda_sat(N,w)
        for sec,(delta,tol) in TARGETS.items():
            Ladm = lambda_adm(N,w,delta,tol)
            eta = Ladm / Lsat
            sec_tex = sec.replace("_", r"\_")
            # IMPORTANT: no stray spaces, one row per line, ends with \\ only
            rows.append(f"{name} & {N} & {sec_tex} & {Ladm} & {Lsat} & {eta:.3f} \\\\")
# IMPORTANT: final newline, but no extra TeX commands
with open(r"..\paper\deformation_eta_table.tex","w", encoding="utf-8", newline="\n") as f:
    f.write("\n".join(rows) + "\n")

# ---------------- FIGURE 1 ----------------
base_w = WEIGHTS["8-15-24"]
cov=[]
for L in range(0, 131):
    cov.append(len(residues(120,L,base_w))/120)

plt.figure()
plt.plot(range(len(cov)), cov, linewidth=2)
plt.xlabel(r"Weighted cutoff $\Lambda$")
plt.ylabel(r"Carrier coverage $|K_\Lambda|/N$")
plt.xlim(0, 120)
plt.ylim(0, 1.02)
plt.tight_layout()
plt.savefig(r"..\paper\fig1_carrier_coverage_N120.png", dpi=300)
plt.close()

# ---------------- FIGURE 2 ----------------
plt.figure()
for name,w in WEIGHTS.items():
    cov=[]
    for L in range(0, 131):
        cov.append(len(residues(120,L,w))/120)
    plt.plot(range(len(cov)), cov, label=name)

plt.xlabel(r"Weighted cutoff $\Lambda$")
plt.ylabel(r"Carrier coverage $|K_\Lambda|/N$")
plt.xlim(0, 120)
plt.ylim(0, 1.02)
plt.legend(title="Weights", fontsize=8)
plt.tight_layout()
plt.savefig(r"..\paper\fig2_deformation_stability_N120.png", dpi=300)
plt.close()

print("All deliverables regenerated:")
print(r"  ..\paper\deformation_eta_table.tex")
print(r"  ..\paper\fig1_carrier_coverage_N120.png")
print(r"  ..\paper\fig2_deformation_stability_N120.png")
