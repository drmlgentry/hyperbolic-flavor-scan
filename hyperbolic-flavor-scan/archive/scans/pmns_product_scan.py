import numpy as np
from scipy.linalg import qr, logm
import snappy
import pandas as pd
import time

PMNS = np.array([[0.821,0.550,0.148],[0.357,0.339,0.871],[0.442,0.762,0.471]])

def to_numpy(m):
    return np.array([[complex(m[0,0]), complex(m[0,1])],
                     [complex(m[1,0]), complex(m[1,1])]])

def axis_from_matrix(M):
    try:
        L = logm(M)
    except:
        return None
    nx = L[0,1].real + L[1,0].real
    ny = (L[1,0].imag - L[0,1].imag)
    nz = (L[0,0] - L[1,1]).real
    n = np.array([nx, ny, nz])
    norm = np.linalg.norm(n)
    if norm < 1e-10: return None
    return n / norm

def symmetric_overlap(axes, sigma):
    n = len(axes)
    O = np.eye(n, dtype=float)
    for i in range(n):
        for j in range(i+1, n):
            cos_ang = np.clip(np.dot(axes[i], axes[j]), -1, 1)
            theta = np.arccos(abs(cos_ang))
            val = np.exp(-theta**2 / (2*sigma**2))
            O[i,j] = O[j,i] = val
    return O

def U_qr(O):
    if not np.all(np.linalg.eigvalsh(O) > 1e-10): return None
    U, _ = qr(O)
    for i in range(3):
        if U[i,i] < 0: U[:,i] *= -1
    return np.abs(U)

def fitness(U):
    return float(np.linalg.norm(U - PMNS, 'fro'))

def generate_words(gens, length):
    from itertools import product as iproduct
    letters = list(gens) + [g.upper() for g in gens]
    for combo in iproduct(letters, repeat=length):
        valid = True
        for k in range(len(combo)-1):
            if combo[k].lower() == combo[k+1].lower() and combo[k] != combo[k+1]:
                valid = False
                break
        if valid:
            yield ''.join(combo)

def best_U_for_manifold(mi, sigmas=[0.3,0.5,0.8,1.0,1.5,2.0,2.5,3.0]):
    mfld = snappy.OrientableClosedCensus[mi]
    mname = mfld.name()
    try:
        G = mfld.fundamental_group()
        gens = G.generators()
    except:
        return mname, None, 1e9, None

    axes, words = [], []
    for wlen in [2, 3]:
        for word in generate_words(gens, wlen):
            try:
                mat = to_numpy(G.SL2C(word))
            except:
                continue
            ax = axis_from_matrix(mat)
            if ax is None: continue
            axes.append(ax)
            words.append(word)
            if len(axes) >= 40: break
        if len(axes) >= 40: break

    if len(axes) < 3: return mname, None, 1e9, None

    best_f, best_U, best_meta = 1e9, None, None
    na = len(axes)
    for i in range(na):
        for j in range(i+1, na):
            for k in range(j+1, na):
                ax3 = [axes[i], axes[j], axes[k]]
                for sigma in sigmas:
                    O = symmetric_overlap(ax3, sigma)
                    U = U_qr(O)
                    if U is None: continue
                    f = fitness(U)
                    if f < best_f:
                        best_f = f
                        best_U = U.copy()
                        best_meta = (words[i], words[j], words[k], sigma)

    return mname, best_U, best_f, best_meta

print('PMNS Product Construction:  PMNS ~ U_e† @ U_nu')
print('='*60)
print('For each manifold pair (i,j): test U_i† @ U_j vs PMNS')
print('Benchmark: 0.018968')
print('='*60)

N = 60  # manifolds to precompute
print(f'Step 1: Computing best U for first {N} manifolds...')
t0 = time.time()

manifold_Us = []
for mi in range(N):
    mname, U, f, meta = best_U_for_manifold(mi)
    manifold_Us.append((mname, U, f, meta))
    if (mi+1) % 10 == 0:
        elapsed = time.time() - t0
        valid = sum(1 for _,U,_,_ in manifold_Us if U is not None)
        print(f'  [{mi+1:3d}/{N} | {elapsed:.0f}s] {valid} valid manifolds')

valid_Us = [(name, U, f, meta) for name, U, f, meta in manifold_Us if U is not None]
print(f'  Done. {len(valid_Us)} valid manifolds with U matrices.')
print()

print('Step 2: Testing all pairs U_i† @ U_j...')
t1 = time.time()
results = []
best_f = 1e9

for i in range(len(valid_Us)):
    for j in range(len(valid_Us)):
        if i == j: continue
        name_i, Ui, fi, meta_i = valid_Us[i]
        name_j, Uj, fj, meta_j = valid_Us[j]
        if Ui is None or Uj is None: continue

        # Test Ui.T @ Uj  (transpose = dagger for real orthogonal)
        P = Ui.T @ Uj
        f = fitness(np.abs(P))
        if f < best_f:
            best_f = f
            print(f'  NEW BEST: {f:.6f}  {name_i}† x {name_j}')
            print(f'    |P|: [{np.abs(P[0,0]):.3f} {np.abs(P[0,1]):.3f} {np.abs(P[0,2]):.3f}]  target: [0.821 0.550 0.148]')
            print(f'         [{np.abs(P[1,0]):.3f} {np.abs(P[1,1]):.3f} {np.abs(P[1,2]):.3f}]  target: [0.357 0.339 0.871]')
            print(f'         [{np.abs(P[2,0]):.3f} {np.abs(P[2,1]):.3f} {np.abs(P[2,2]):.3f}]  target: [0.442 0.762 0.471]')
            print(f'    U_e from {name_i}: words={meta_i[:3]} sigma={meta_i[3]} fit={fi:.4f}')
            print(f'    U_nu from {name_j}: words={meta_j[:3]} sigma={meta_j[3]} fit={fj:.4f}')
            print()

        results.append({
            'manifold_e': name_i, 'manifold_nu': name_j,
            'fitness': round(f, 6),
            'fit_e': round(fi, 4), 'fit_nu': round(fj, 4)
        })

elapsed = time.time() - t1
print(f'Pair scan done in {elapsed:.0f}s')
print()
print('='*60)
print('TOP 20 PAIRS')
print('='*60)
top20 = sorted(results, key=lambda x: x['fitness'])[:20]
for r in top20:
    print(f'  {r["manifold_e"]:8s}† x {r["manifold_nu"]:8s}  fit={r["fitness"]:.6f}  (fit_e={r["fit_e"]} fit_nu={r["fit_nu"]})')

df = pd.DataFrame(results)
df.to_csv('pmns_product_scan_results.csv', index=False)
print()
print('Saved: pmns_product_scan_results.csv')