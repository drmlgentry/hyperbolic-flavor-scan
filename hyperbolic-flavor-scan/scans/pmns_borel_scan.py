import numpy as np
from scipy.linalg import qr, logm
from itertools import permutations
import snappy
import pandas as pd
import time

PMNS = np.array([[0.821,0.550,0.148],[0.357,0.339,0.871],[0.442,0.762,0.471]])

PERMS = list(permutations([0,1,2]))

def fitness(U):
    # Best over all column permutations and sign flips
    best = 1e9
    for perm in PERMS:
        Up = U[:, perm]
        best = min(best, float(np.linalg.norm(Up - PMNS, 'fro')))
    return best

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

def U_from_L(L):
    try:
        if not np.all(np.isfinite(L)): return None
        U, _ = qr(L)
        for i in range(3):
            if U[i,i] < 0: U[:,i] *= -1
        return np.abs(U)
    except:
        return None

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

print('PMNS Borel/Triangular Construction Scan (permutation-aware)')
print('='*60)
print('Target L lower triangle: l21=0.443, l31=0.530, l32=-0.381')
print('Benchmark: 0.018968')
print()

# Verify L_TARGET with permutation-aware fitness
L_TARGET = np.array([[1.0,   0.0,   0.0],
                     [0.443, 1.0,   0.0],
                     [0.530,-0.381, 1.0]])

print('Verification: U_QR(L_TARGET) vs PMNS (permutation-aware)')
U_check = U_from_L(L_TARGET)
if U_check is not None:
    f_raw = float(np.linalg.norm(U_check - PMNS, 'fro'))
    f_perm = fitness(U_check)
    print(f'  raw fitness = {f_raw:.6f}')
    print(f'  perm fitness = {f_perm:.6f}')
    # Find best permutation
    for perm in PERMS:
        Up = U_check[:, perm]
        f = float(np.linalg.norm(Up - PMNS, 'fro'))
        if abs(f - f_perm) < 1e-6:
            print(f'  best perm: {perm}')
            for row, trow in zip(Up, PMNS):
                print(f'  [{row[0]:.3f} {row[1]:.3f} {row[2]:.3f}]  target [{trow[0]:.3f} {trow[1]:.3f} {trow[2]:.3f}]')
            break
print()

scales = np.linspace(0.1, 5.0, 50)
signs = [(1,1,1),(1,1,-1),(1,-1,1),(-1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1),(-1,-1,-1)]

census = snappy.OrientableClosedCensus
N_manifolds = 100
results = []
t0 = time.time()

for mi in range(N_manifolds):
    mfld = census[mi]
    mname = mfld.name()
    try:
        G = mfld.fundamental_group()
        gens = G.generators()
    except:
        continue

    axes, mats, words = [], [], []
    for wlen in [2, 3]:
        for word in generate_words(gens, wlen):
            try:
                mat = to_numpy(G.SL2C(word))
            except:
                continue
            ax = axis_from_matrix(mat)
            if ax is None: continue
            axes.append(ax)
            mats.append(mat)
            words.append(word)
            if len(axes) >= 40: break
        if len(axes) >= 40: break

    if len(axes) < 3: continue

    best_f = 1e9
    best_row = None
    na = len(axes)

    for i in range(na):
        for j in range(i+1, na):
            for k in range(j+1, na):
                ax3 = [axes[i], axes[j], axes[k]]

                for scale in scales:
                    for s1,s2,s3 in signs:
                        L = np.array([[1.0, 0.0, 0.0],
                                      [s1*scale*np.dot(ax3[1], ax3[0]), 1.0, 0.0],
                                      [s2*scale*np.dot(ax3[2], ax3[0]),
                                       s3*scale*np.dot(ax3[2], ax3[1]), 1.0]])
                        U = U_from_L(L)
                        if U is None: continue
                        f = fitness(U)
                        if f < best_f:
                            best_f = f
                            best_row = {
                                'manifold': mname,
                                'words': f'{words[i]}/{words[j]}/{words[k]}',
                                'scale': round(scale,3),
                                'signs': f'{s1},{s2},{s3}',
                                'fitness': round(f,6),
                                'l21': round(L[1,0],4),
                                'l31': round(L[2,0],4),
                                'l32': round(L[2,1],4),
                                'U': U.copy()
                            }

    if best_row:
        results.append(best_row)
        if best_f < 0.08:
            elapsed = time.time() - t0
            U = best_row['U']
            print(f'[{mi:3d} | {elapsed:5.1f}s] {mname}  fit={best_f:.6f}  scale={best_row["scale"]}  {best_row["words"]}')
            print(f'  L: l21={best_row["l21"]}  l31={best_row["l31"]}  l32={best_row["l32"]}')
            print(f'  tgt:   0.443         0.530         -0.381')
            # show best permutation
            for perm in PERMS:
                Up = U[:, perm]
                f = float(np.linalg.norm(Up - PMNS, 'fro'))
                if abs(f - best_f) < 1e-6:
                    print(f'  best perm {perm}:')
                    for row, trow in zip(Up, PMNS):
                        print(f'    [{row[0]:.3f} {row[1]:.3f} {row[2]:.3f}]  target [{trow[0]:.3f} {trow[1]:.3f} {trow[2]:.3f}]')
                    break
            print()

    if (mi+1) % 10 == 0:
        elapsed = time.time() - t0
        top = sorted(results, key=lambda x: x['fitness'])[:3]
        print(f'--- [{mi+1:3d}/{N_manifolds} | {elapsed:.0f}s] top3: ' +
              ' | '.join(f'{r["manifold"]}={r["fitness"]:.4f}' for r in top) + ' ---')

print()
print('='*60)
print('TOP 20')
print('='*60)
top20 = sorted(results, key=lambda x: x['fitness'])[:20]
for r in top20:
    U = r['U']
    print(f'{r["manifold"]:10s}  fit={r["fitness"]:.6f}  scale={r["scale"]}  signs={r["signs"]}  {r["words"]}')
    print(f'  L: l21={r["l21"]}  l31={r["l31"]}  l32={r["l32"]}  (target: 0.443 / 0.530 / -0.381)')
    for perm in PERMS:
        Up = U[:, perm]
        f = float(np.linalg.norm(Up - PMNS, 'fro'))
        if abs(f - r["fitness"]) < 1e-6:
            for row, trow in zip(Up, PMNS):
                print(f'  [{row[0]:.3f} {row[1]:.3f} {row[2]:.3f}]  target [{trow[0]:.3f} {trow[1]:.3f} {trow[2]:.3f}]')
            break

df = pd.DataFrame([{k:v for k,v in r.items() if k != 'U'} for r in results])
df.to_csv('pmns_borel_scan_results.csv', index=False)
print()
print('Saved: pmns_borel_scan_results.csv')