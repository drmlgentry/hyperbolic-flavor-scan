import snappy
import numpy as np
from scipy.linalg import logm as scipy_logm
import csv

def axis_imlogm(A):
    L = scipy_logm(A)
    imL = L.imag
    nx = float((imL[0,1] + imL[1,0]).real)
    ny = float((1j*(imL[1,0] - imL[0,1])).real)
    nz = float((imL[0,0] - imL[1,1]).real)
    n = np.array([nx, ny, nz])
    norm = np.linalg.norm(n)
    if norm < 1e-10:
        return None
    return n / norm

def axis_eigenvector(A):
    """Alternative: fixed points of Mobius action on CP^1 -> axis on S^2."""
    # Loxodromic A has two fixed points z+, z- on CP^1
    # Map CP^1 -> S^2 via stereographic projection
    # Axis = midpoint direction between the two stereographic images
    a, b, c, d = A[0,0], A[0,1], A[1,0], A[1,1]
    # Fixed points: az+b = z(cz+d) => cz^2 + (d-a)z - b = 0
    disc = (d - a)**2 + 4*b*c
    if abs(disc) < 1e-12:
        return None
    sq = np.sqrt(disc)
    z1 = ((a - d) + sq) / (2*c) if abs(c) > 1e-12 else None
    z2 = ((a - d) - sq) / (2*c) if abs(c) > 1e-12 else None
    if z1 is None or z2 is None:
        return None
    # Stereographic projection CP^1 -> S^2: z -> (2Re(z), 2Im(z), |z|^2-1)/(|z|^2+1)
    def stereo(z):
        r2 = abs(z)**2
        return np.array([2*z.real, 2*z.imag, r2-1]) / (r2+1)
    p1 = stereo(z1)
    p2 = stereo(z2)
    # Axis = direction of p1 - p2 (the line connecting fixed points on S^2)
    n = p1 - p2
    norm = np.linalg.norm(n)
    if norm < 1e-10:
        return None
    return n / norm

def get_axes(M, w1, w2):
    G = M.fundamental_group()
    results = {}
    for w, method_fn in [(w1, None), (w2, None)]:
        try:
            mat = G.SL2C(w)
            A = np.array([[complex(mat[i][j]) for j in range(2)] for i in range(2)])
            A /= np.sqrt(np.linalg.det(A))
            results[w] = A
        except:
            results[w] = None
    return results

def dot_product(M, w1, w2, method='imlogm'):
    G = M.fundamental_group()
    try:
        A1 = np.array([[complex(G.SL2C(w1)[i][j]) for j in range(2)] for i in range(2)])
        A2 = np.array([[complex(G.SL2C(w2)[i][j]) for j in range(2)] for i in range(2)])
        A1 /= np.sqrt(np.linalg.det(A1))
        A2 /= np.sqrt(np.linalg.det(A2))
        if method == 'imlogm':
            n1, n2 = axis_imlogm(A1), axis_imlogm(A2)
        else:
            n1, n2 = axis_eigenvector(A1), axis_eigenvector(A2)
        if n1 is None or n2 is None:
            return None, n1, n2
        return abs(float(np.dot(n1, n2))), n1, n2
    except Exception as e:
        return None, None, None

COS_TW = 80.377 / 91.1876
print(f"cos(theta_W) = {COS_TW:.8f}", flush=True)
print()

m006 = snappy.OrientableClosedCensus[43]
m003 = snappy.OrientableClosedCensus[1]

# === LINE 1: METHOD COMPARISON ===
print("=== LINE 1: METHOD COMPARISON ===", flush=True)
for label, M, w1, w2 in [("m006", m006, "AAb","BaaB"),
                           ("m003", m003, "BAB","BBAb")]:
    d1, n1a, n1b = dot_product(M, w1, w2, 'imlogm')
    d2, n2a, n2b = dot_product(M, w1, w2, 'eigen')
    s1 = f"{d1:.8f}  diff={d1-COS_TW:+.8f}" if d1 is not None else "FAILED"
    s2 = f"{d2:.8f}  diff={d2-COS_TW:+.8f}" if d2 is not None else "FAILED"
    print(f"  {label} {w1}/{w2}:")
    print(f"    imlogm:      |dot|={s1}")
    print(f"    eigenvector: |dot|={s2}")
    if d1 and d2:
        print(f"    Methods agree to: {abs(d1-d2):.2e}")
print()

# === XZ-PLANE INVARIANCE TEST ===
print("=== XZ-PLANE INVARIANCE: IS y=0 FOR ALL pi_1(m006) GENERATORS? ===", flush=True)
print("Testing axis_imlogm(g) for all generators and short words in m006", flush=True)
G6 = m006.fundamental_group()
test_words = ['a','b','A','B','aa','ab','aB','Ab','AB','ba','bA','ba',
              'aaa','aab','aaB','AbA','AAb','BaaB','Baa','bAAb']
max_y = 0.0
for w in test_words:
    try:
        mat = G6.SL2C(w)
        A = np.array([[complex(mat[i][j]) for j in range(2)] for i in range(2)])
        A /= np.sqrt(np.linalg.det(A))
        n = axis_imlogm(A)
        if n is not None:
            y_component = abs(n[1])
            max_y = max(max_y, y_component)
            flag = " *** NON-ZERO Y ***" if y_component > 1e-6 else ""
            print(f"  {w:<8} y={n[1]:.2e}  axis=[{n[0]:.4f}, {n[1]:.6f}, {n[2]:.4f}]{flag}")
    except Exception as e:
        print(f"  {w:<8} ERROR: {e}")

print(f"\n  Maximum |y| across all tested words: {max_y:.2e}")
if max_y < 1e-6:
    print("  CONFIRMED: xz-plane is invariant for all tested words")
    print("  This suggests the holonomy of m006 preserves the xz-plane")
    print("  => Holonomy representation lies in a subgroup preserving y=0")
else:
    print("  xz-plane is NOT universally invariant")
print()

# === LINE 2: CONJUGACY CLASS TEST ===
print("=== LINE 2: CONJUGACY INVARIANCE ===", flush=True)
conj_words = []
for g in ['a','b','A','B','aa','ab','aB']:
    gi = g.swapcase() if len(g)==1 else g[::-1].swapcase()
    conj_words.append((g, g + "AAb" + gi))

print(f"  {'Conjugator':<8} {'Conjugated word':<16} {'|dot w/ BaaB|':>14}  {'diff':>10}")
for g, w_conj in conj_words:
    d, _, _ = dot_product(m006, w_conj, "BaaB", 'imlogm')
    if d is not None:
        print(f"  {g:<8} {w_conj:<16} {d:>14.8f}  {d-COS_TW:>+10.8f}")
    else:
        print(f"  {g:<8} {w_conj:<16} {'FAILED':>14}")
print()

# === LINE 5: DEHN SURGERY NULL TEST ===
print("=== LINE 5: DEHN SURGERY NULL TEST ===", flush=True)
from math import gcd
cusped = snappy.Manifold("m006")
results_surgery = []
count = 0
for p in range(-12, 13):
    for q in range(1, 7):
        if p == 0: continue
        if gcd(abs(p), q) != 1: continue
        try:
            M = cusped.copy()
            M.dehn_fill((p, q))
            vol = float(M.volume())
            if vol < 0.4: continue
            d, _, _ = dot_product(M, "AAb", "BaaB", 'imlogm')
            if d is not None:
                results_surgery.append((p, q, vol, d))
                count += 1
        except:
            pass

print(f"  Computed {count} valid fillings", flush=True)
results_surgery.sort(key=lambda x: abs(x[3]-COS_TW))

print(f"\n  Fillings CLOSEST to cos(theta_W)={COS_TW:.6f}:")
print(f"  {'(p,q)':<10} {'vol':>10}  {'|dot|':>10}  {'diff':>12}  {'angle':>8}")
print("  " + "-"*55)
for p,q,vol,d in results_surgery[:15]:
    angle = np.degrees(np.arccos(min(d,1.0)))
    marker = " <-- CKM" if p==-5 and q==2 else ""
    print(f"  ({p:3d},{q:2d})   {vol:>10.6f}  {d:>10.6f}  {d-COS_TW:>+12.8f}  {angle:>8.4f}{marker}")

# Rank of CKM filling
ckm_rank = next((i+1 for i,(p,q,v,d) in enumerate(results_surgery) if p==-5 and q==2), None)
total = len(results_surgery)
print(f"\n  CKM filling (-5,2) rank: {ckm_rank} out of {total} fillings")
print(f"  (rank 1 = closest to cos(theta_W))")

# Save CSV
with open('data/weinberg_surgery_scan.csv','w',newline='') as f:
    w = csv.writer(f)
    w.writerow(['p','q','vol','dot','diff_from_cosTW','angle_deg'])
    for p,q,vol,d in sorted(results_surgery, key=lambda x: (x[0],x[1])):
        w.writerow([p,q,f"{vol:.8f}",f"{d:.8f}",
                    f"{d-COS_TW:.8f}",f"{np.degrees(np.arccos(min(d,1.0))):.6f}"])
print(f"  Saved to data/weinberg_surgery_scan.csv", flush=True)
