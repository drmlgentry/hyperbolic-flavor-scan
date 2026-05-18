#!/usr/bin/env sage
# unified_fermion_witnesses.sage
# Finds best approximations of m/m_e by norms in:
#   - Eisenstein integers Z[ω] for leptons (Q(√-3))
#   - Z[√17] for quarks (Q(√17))
# Outputs LaTeX table and verification.

from sage.all import *

# ------------------------------------------------------------
# Fermion mass ratios (PDG 2024) relative to electron
# ------------------------------------------------------------
ratios = {
    'mu'   : 206.771,      # muon
    'tau'  : 3477.22,      # tau
    'u'    : 4.2361,       # up
    'd'    : 9.1368,       # down  (approx, actual 9.14)
    's'    : 182.78,       # strange
    'c'    : 2495.11,      # charm
    'b'    : 8180.04,      # bottom
    't'    : 338082.19,    # top
}

# ------------------------------------------------------------
# Search in Eisenstein integers Z[ω], norm N = a² - a*b + b²
# ------------------------------------------------------------
def eisenstein_witness(target, bound=200):
    best_a, best_b = 0, 0
    best_norm = 0
    best_err = float('inf')
    for a in range(-bound, bound+1):
        for b in range(-bound, bound+1):
            norm = a*a - a*b + b*b
            if norm == 0: continue
            err = abs(norm - target) / target
            if err < best_err:
                best_err = err
                best_norm = norm
                best_a, best_b = a, b
    return best_a, best_b, best_norm, best_err

# ------------------------------------------------------------
# Search in Z[√17], norm N = |a² - 17*b²|
# ------------------------------------------------------------
def sqrt17_witness(target, bound=200):
    best_a, best_b = 0, 0
    best_norm = 0
    best_err = float('inf')
    for a in range(-bound, bound+1):
        for b in range(-bound, bound+1):
            norm = abs(a*a - 17*b*b)
            if norm == 0: continue
            err = abs(norm - target) / target
            if err < best_err:
                best_err = err
                best_norm = norm
                best_a, best_b = a, b
    return best_a, best_b, best_norm, best_err

# ------------------------------------------------------------
# LaTeX table generation
# ------------------------------------------------------------
def latex_table():
    lines = []
    lines.append("\\begin{tabular}{lcccccc}")
    lines.append("\\toprule")
    lines.append("Particle & $m/m_e$ & Field & Norm $N$ & Witness $(a,b)$ & Error \\\\")
    lines.append("\\midrule")
    
    # Leptons (Eisenstein)
    for p in ['mu', 'tau']:
        target = ratios[p]
        a,b,norm,err = eisenstein_witness(target, bound=300)
        lines.append(f"{p} & {target:.2f} & $\\mathbb{{Z}}[\\omega]$ & ${norm}$ & $({a},{b})$ & ${err*100:.3f}\\%$ \\\\")
    
    # Quarks (Z[√17])
    for p in ['u','d','s','c','b','t']:
        target = ratios[p]
        a,b,norm,err = sqrt17_witness(target, bound=300)
        lines.append(f"{p} & {target:.2f} & $\\mathbb{{Z}}[\\sqrt{{17}}]$ & ${norm}$ & $({a},{b})$ & ${err*100:.3f}\\%$ \\\\")
    
    lines.append("\\bottomrule")
    lines.append("\\end{tabular}")
    return "\n".join(lines)

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
print("="*70)
print("UNIFIED ARITHMETIC WITNESSES FOR FERMION MASS RATIOS")
print("="*70)
print()
print("Leptons (Eisenstein integers Z[ω], norm a²-ab+b²):")
for p in ['mu','tau']:
    a,b,norm,err = eisenstein_witness(ratios[p], bound=300)
    print(f"  {p:3s}: target={ratios[p]:8.2f} → N({a:3d},{b:3d}) = {norm:5d}  error={err*100:.4f}%")

print("\nQuarks (Z[√17], norm |a²-17b²|):")
for p in ['u','d','s','c','b','t']:
    a,b,norm,err = sqrt17_witness(ratios[p], bound=300)
    print(f"  {p:3s}: target={ratios[p]:8.2f} → N({a:3d},{b:3d}) = {norm:6d}  error={err*100:.4f}%")

print("\n" + "="*70)
print("LaTeX table:")
print("="*70)
print(latex_table())