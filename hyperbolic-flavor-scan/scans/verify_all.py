"""
verify_all.py -- Reproduce key numerical results from Papers 5-8.
Usage: python scans/verify_all.py
"""
import snappy, numpy as np

def check(label, computed, expected, tol_pct=1.0):
    err = 100*abs(computed-expected)/abs(expected)
    status = "PASS" if err < tol_pct else "FAIL"
    print(f"  [{status}] {label}: {computed:.6f} vs {expected:.6f} ({err:.3f}%)")
    return status == "PASS"

print("="*60)
print("Hyperbolic Flavor Geometry -- Verification")
print("Manifolds: m003=Meyerhoff (idx 1), m006 (idx 43)")
print("="*60)

M3 = snappy.OrientableClosedCensus[1]   # Meyerhoff manifold
M6 = snappy.OrientableClosedCensus[43]  # m006
G3, G6 = M3.fundamental_group(), M6.fundamental_group()

def phi(G, w):
    L = np.array(G.SL2C(w), dtype=complex)
    ev = np.linalg.eigvals(L)
    lam = ev[np.argmax(np.abs(ev))]
    return float(np.degrees(np.angle(lam)))

def modlam(G, w):
    L = np.array(G.SL2C(w), dtype=complex)
    return float(np.max(np.abs(np.linalg.eigvals(L))))

print(f"\nm003 (Meyerhoff, vol={float(M3.volume()):.4f}, H1={M3.homology()}):")
check("180-phi(aa) ~ delta_CKM=68.0 -- wait, this is m006")

print(f"\nm006 (vol={float(M6.volume()):.4f}, H1={M6.homology()}):")
check("180-phi(aa) ~ delta_CKM=68.0", 180-abs(phi(G6,"aa")), 68.0, tol_pct=1.0)
check("phi(b) ~ 89.16 deg", abs(phi(G6,"b")), 89.156, tol_pct=0.5)

print(f"\nm003 (Meyerhoff):")
check("|lam(bbbb)|/|lam(bAbA)| ~ mb/mc=3.2913",
      modlam(G3,"bbbb")/modlam(G3,"bAbA"), 3.2913, tol_pct=0.1)
check("180-phi(aabb) ~ theta23_nu=49.1",
      180-abs(phi(G3,"aabb")), 49.1, tol_pct=3.0)
cp_sum = phi(G3,"aa") - phi(G3,"ab") + phi(G3,"aB")
# normalize to [0,360]
while cp_sum < 0: cp_sum += 360
while cp_sum > 360: cp_sum -= 360
check("phi_aa-phi_ab+phi_aB ~ delta_CP=197", cp_sum, 197.0, tol_pct=5.0)

print("\nDone.")
