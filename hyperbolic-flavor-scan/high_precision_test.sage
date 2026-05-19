#!/usr/bin/env sage
"""
high_precision_test.sage
========================
Definitive test: compute |lambda_b^2/lambda_a|^2 at 1000-bit precision
and compare directly to mb/mc = 4.18/1.27.

If they agree to 900+ bits: the relation is EXACT (theorem)
If they disagree at bit 50: the relation is APPROXIMATE (coincidence)

Also fix the eigenvalue branch issue from the previous script.

Run: conda run -n sage sage high_precision_test.sage
"""

import snappy

print("="*65)
print("HIGH PRECISION DEFINITIVE TEST")
print("="*65)
print()

# Use 1000-bit precision
prec = 1000
CC = ComplexField(prec)
RR = RealField(prec)

M_P = snappy.OrientableClosedCensus[1]
rho_P = M_P.polished_holonomy()

def get_mat_hp(rho, word, prec=1000):
    """Get normalized holonomy matrix at high precision."""
    CC_hp = ComplexField(prec)
    mat = matrix(CC_hp, [[CC_hp(rho(word)[i,j]) for j in range(2)]
                          for i in range(2)])
    d = mat.det()
    return mat / d.sqrt()

def dominant_eig_hp(rho, word, prec=1000):
    """
    Get dominant eigenvalue with correct branch selection.
    Dominant = |lambda| > 1 for loxodromic elements.
    """
    CC_hp = ComplexField(prec)
    mat = get_mat_hp(rho, word, prec)
    eigs = mat.eigenvalues()
    # Select eigenvalue with |lambda| > 1
    lam_dom = max(eigs, key=lambda x: abs(x))
    lam_sub = min(eigs, key=lambda x: abs(x))
    
    # Verify: |lam_dom| * |lam_sub| should = 1 (det = 1)
    product = abs(lam_dom) * abs(lam_sub)
    
    return lam_dom, lam_sub, product

print("Computing at 1000-bit precision...")
print()

lam_a_dom, lam_a_sub, prod_a = dominant_eig_hp(rho_P, 'a')
lam_b_dom, lam_b_sub, prod_b = dominant_eig_hp(rho_P, 'b')

print(f"Generator a:")
print(f"  |lambda_dom| = {abs(lam_a_dom)}")
print(f"  |lambda_sub| = {abs(lam_a_sub)}")
print(f"  |lam_dom|*|lam_sub| = {prod_a}  (should be 1.000...)")
print(f"  ell_a = 2*log|lam_dom| = {2*log(abs(lam_a_dom))}")
print()
print(f"Generator b:")
print(f"  |lambda_dom| = {abs(lam_b_dom)}")
print(f"  |lambda_sub| = {abs(lam_b_sub)}")
print(f"  |lam_dom|*|lam_sub| = {prod_b}  (should be 1.000...)")
print(f"  ell_b = 2*log|lam_dom| = {2*log(abs(lam_b_dom))}")
print()

ell_a = 2*log(abs(lam_a_dom))
ell_b = 2*log(abs(lam_b_dom))
print(f"Canonical values:")
print(f"  ell_a = {RR(ell_a)}")
print(f"  ell_b = {RR(ell_b)}")
print()

print("="*65)
print("CALIBRATION: tr(aa on m006) vs 3-sqrt(17)")
print("="*65)
print()

M_C = snappy.OrientableClosedCensus[43]
rho_C = M_C.polished_holonomy()
mat_aa_C = get_mat_hp(rho_C, 'aa', prec)
tr_aa_C = mat_aa_C.trace()

three_minus_sqrt17 = RR(3) - RR(17).sqrt()

re_tr = real(tr_aa_C)
print(f"Re(tr(aa on m006))  = {re_tr}")
print(f"3 - sqrt(17)        = {three_minus_sqrt17}")
diff_calib = abs(re_tr - three_minus_sqrt17)
print(f"Difference          = {diff_calib}")
print(f"Relative error      = {diff_calib/abs(three_minus_sqrt17)}")
print()

n_correct_bits_calib = -log(diff_calib/abs(three_minus_sqrt17))/log(RR(2))
print(f"Calibration: agreement to ~{float(n_correct_bits_calib):.0f} bits")
print()

print("="*65)
print("THE DEFINITIVE TEST: |lambda_b^2/lambda_a|^2 vs mb/mc")
print("="*65)
print()

ratio = lam_b_dom^2 / lam_a_dom
ratio_sq = abs(ratio)^2

mb_mc_pdg = RR(4.18) / RR(1.27)

print(f"|lambda_b^2/lambda_a|^2 = {ratio_sq}")
print()
print(f"mb/mc (PDG 2024)        = {mb_mc_pdg}")
print()
diff = abs(ratio_sq - mb_mc_pdg)
rel_err = diff / mb_mc_pdg
print(f"Absolute difference     = {diff}")
print(f"Relative error          = {rel_err}")
print(f"Error in percent        = {float(rel_err)*100}%")
print()

n_correct_bits = -log(rel_err)/log(RR(2))
print(f"Agreement to ~{float(n_correct_bits):.1f} bits")
print()

print("INTERPRETATION:")
if float(n_correct_bits) > 500:
    print("  EXACT RELATION: agrees to 500+ bits")
    print("  |lambda_b^2/lambda_a|^2 = mb/mc is a THEOREM")
elif float(n_correct_bits) > float(n_correct_bits_calib) * 0.8:
    print(f"  MATCHES CALIBRATION PRECISION ({float(n_correct_bits_calib):.0f} bits)")
    print("  Cannot distinguish exact from approximate at this precision")
    print("  The relation may or may not be exact")
else:
    print(f"  APPROXIMATE: agrees to only {float(n_correct_bits):.0f} bits")
    print(f"  vs calibration precision of {float(n_correct_bits_calib):.0f} bits")
    print("  The match is NUMERICAL COINCIDENCE, not exact")
print()

print("="*65)
print("EXTENDED PRECISION TABLE")
print("="*65)
print()
print("Showing first 60 significant digits of each:")
print()
print(f"  ratio^2 = {ratio_sq:.60f}")
print(f"  mb/mc   = {mb_mc_pdg:.60f}")
print()

# Find where they first differ
r_str = f"{ratio_sq:.200f}"
m_str = f"{mb_mc_pdg:.200f}"
first_diff = None
for i,(a,b) in enumerate(zip(r_str, m_str)):
    if a != b:
        first_diff = i
        break
if first_diff:
    print(f"  First digit difference at position {first_diff}")
    print(f"  ratio^2: ...{r_str[max(0,first_diff-3):first_diff+5]}...")
    print(f"  mb/mc:   ...{m_str[max(0,first_diff-3):first_diff+5]}...")
print()

print("="*65)
print("CROSS-CHECK: mtau/mmu on m006")
print("="*65)
print()

lam_a_C_dom, _, _ = dominant_eig_hp(rho_C, 'a', prec)
ell_a_C = 2*log(abs(lam_a_C_dom))
exp_3a = exp(3*ell_a_C)
mtau_mmu = RR(1.7769)/RR(0.10566)

print(f"exp(3*ell_a on m006) = {exp_3a}")
print(f"mtau/mmu (PDG)       = {mtau_mmu}")
diff2 = abs(exp_3a - mtau_mmu)
rel_err2 = diff2/mtau_mmu
print(f"Relative error       = {float(rel_err2)*100}%")
n_bits2 = -log(rel_err2)/log(RR(2))
print(f"Agreement to ~{float(n_bits2):.1f} bits")
print()

print("Also checking m003 mtau/mmu:")
lam_b_P_dom = lam_b_dom
ell_b_P = ell_b
exp_2a_b = exp(2*ell_a + ell_b_P)
print(f"exp(2*ell_a + ell_b on m003) = {exp_2a_b}")
print(f"mtau/mmu (PDG)               = {mtau_mmu}")
diff3 = abs(exp_2a_b - mtau_mmu)
rel_err3 = diff3/mtau_mmu
print(f"Relative error               = {float(rel_err3)*100}%")
n_bits3 = -log(rel_err3)/log(RR(2))
print(f"Agreement to ~{float(n_bits3):.1f} bits")
print()

print("="*65)
print("FINAL VERDICT")
print("="*65)
print()
print(f"Calibration (known exact tr(aa)=3-sqrt(17)): {float(n_correct_bits_calib):.0f} bits")
print(f"mb/mc test:                                  {float(n_correct_bits):.1f} bits")
print(f"mtau/mmu on m006:                            {float(n_bits2):.1f} bits")
print(f"mtau/mmu on m003:                            {float(n_bits3):.1f} bits")
print()
print("If mb/mc bits << calibration bits: APPROXIMATE")
print("If mb/mc bits ~ calibration bits: PRECISION LIMITED (inconclusive)")
print("If mb/mc bits >> calibration bits: EXACT (strong evidence)")
