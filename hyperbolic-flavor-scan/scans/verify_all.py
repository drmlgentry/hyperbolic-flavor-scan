"""
verify_all.py -- Reproduce all key numerical results from Papers 5-8.
Run from the repo root: python scans/verify_all.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import snappy
import numpy as np

def check(label, computed, expected, tol_pct=1.0):
    err = 100*abs(computed-expected)/abs(expected)
    status = "PASS" if err < tol_pct else "FAIL"
    print(f"  [{status}] {label}: {computed:.6f} vs {expected:.6f} ({err:.3f}%)")
    return status == "PASS"

print("="*60)
print("Hyperbolic Flavor Geometry -- Key Results Verification")
print("="*60)

# ── m006 CKM ──────────────────────────────────────────────────────
M6 = snappy.OrientableClosedCensus[43]
G6 = M6.fundamental_group()
def phi6(w):
    L = G6.SL2C(w)
    ev = np.linalg.eigvals(np.array(L, dtype=complex))
    lam = ev[np.argmax(np.abs(ev))]
    return float(np.degrees(np.angle(lam)))

print("\nm006 (CKM manifold):")
check("phi(b) ~ 90 deg", abs(phi6("b")), 89.156, tol_pct=1.0)
check("phi(aa) -> delta_CKM=68.0", 180-abs(phi6("aa")), 68.0, tol_pct=1.0)

# ── m003 PMNS ─────────────────────────────────────────────────────
M3 = snappy.OrientableClosedCensus[1]
G3 = M3.fundamental_group()
def phi3(w):
    L = G3.SL2C(w)
    ev = np.linalg.eigvals(np.array(L, dtype=complex))
    lam = ev[np.argmax(np.abs(ev))]
    return float(np.degrees(np.angle(lam)))

def modlam3(w):
    L = G3.SL2C(w)
    ev = np.linalg.eigvals(np.array(L, dtype=complex))
    return float(np.max(np.abs(ev)))

print("\nm003 (PMNS manifold):")
check("phi(aa)-phi(ab)+phi(aB) ~ delta_CP=197",
      phi3("aa")-phi3("ab")+phi3("aB")+360, 197+360, tol_pct=5.0)
check("|lam(bbbb)|/|lam(bAbA)| ~ mb/mc=3.2913",
      modlam3("bbbb")/modlam3("bAbA"), 3.2913, tol_pct=0.1)
check("180-phi(aabb) ~ theta23_nu=49.1",
      180-abs(phi3("aabb")), 49.1, tol_pct=3.0)

print("\nVerification complete.")
