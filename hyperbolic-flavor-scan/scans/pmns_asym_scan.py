import numpy as np
from scipy.linalg import qr, logm
from itertools import product as iproduct
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

def twist_angle(M):
    tr = M[0,0] + M[1,1]
    disc = np.sqrt(tr**2 - 4.0 + 0j)
    lam = (tr + disc) / 2
    return float(np.angle(lam))

def asymmetric_overlap(axes, twists, sigma):
    n = len(axes)
    O = np.eye(n, dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j: continue
            cos_ang = np.clip(np.dot(axes[i], axes[j]), -1, 1)
            theta = np.arccos(abs(cos_ang))
            gauss = np.exp(-theta**2 / (2*sigma**2))
            phase = np.cos(twists[i] - twists[j])
            O[i,j] = gauss * phase
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
    letters = list(gens) + [g.upper() for g in gens]
    for combo in iproduct(letters, repeat=length):
        valid = True
        for k in range(len(combo)-1):
            if combo[k].lower() == combo[k+1].lower() and combo[k] != combo[k+1]:
                valid = False
                break
        if valid:
            yield ''.join(combo)

print('PMNS Asymmetric Overlap Scan')
print('='*60)
print('Benchmark: 0.018968 (theoretical min, asymmetric QR)')
print('='*60)

census = snappy.OrientableClosedCensus
N_manifolds = 150
sigmas = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0]
word_lengths = [2, 3]

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

    axes, twists, words = [], [], []

    for wlen in word_lengths:
        for word in generate_words(gens, wlen):
            try:
                mat = to_numpy(G.SL2C(word))
            except:
                continue
            ax = axis_from_matrix(mat)
            if ax is None: continue
            tw = twist_angle(mat)
            axes.append(ax)
            twists.append(tw)
            words.append(word)
            if len(axes) >= 60: break
        if len(axes) >= 60: break

    if len(axes) < 3: continue

    best_f = 1e9
    best_row = None

    na = len(axes)
    for i in range(na):
        for j in range(i+1, na):
            for k in range(j+1, na):
                ax3 = [axes[i], axes[j], axes[k]]
                tw3 = [twists[i], twists[j], twists[k]]
                for sigma in sigmas:
                    O = asymmetric_overlap(ax3, tw3, sigma)
                    U = U_qr(O)
                    if U is None: continue
                    f = fitness(U)
                    if f < best_f:
                        best_f = f
                        best_row = {
                            'manifold': mname,
                            'words': f'{words[i]}/{words[j]}/{words[k]}',
                            'sigma': sigma,
                            'fitness': round(f, 6),
                            'U00': round(U[0,0],4), 'U01': round(U[0,1],4), 'U02': round(U[0,2],4),
                            'U10': round(U[1,0],4), 'U11': round(U[1,1],4), 'U12': round(U[1,2],4),
                            'U20': round(U[2,0],4), 'U21': round(U[2,1],4), 'U22': round(U[2,2],4),
                        }

    if best_row:
        results.append(best_row)
        if best_f < 0.12:
            elapsed = time.time() - t0
            U = np.array([[best_row[f'U{r}{c}'] for c in range(3)] for r in range(3)])
            print(f'[{mi:3d} | {elapsed:5.1f}s] {mname}  fit={best_f:.6f}  sigma={best_row["sigma"]}  {best_row["words"]}')
            print(f'  [{U[0,0]:.3f} {U[0,1]:.3f} {U[0,2]:.3f}]  target: [0.821 0.550 0.148]')
            print(f'  [{U[1,0]:.3f} {U[1,1]:.3f} {U[1,2]:.3f}]  target: [0.357 0.339 0.871]')
            print(f'  [{U[2,0]:.3f} {U[2,1]:.3f} {U[2,2]:.3f}]  target: [0.442 0.762 0.471]')
            print()

    if (mi+1) % 15 == 0:
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
    U = np.array([[r[f'U{row}{col}'] for col in range(3)] for row in range(3)])
    print(f'{r["manifold"]:10s}  fit={r["fitness"]:.6f}  s={r["sigma"]}  {r["words"]}')
    print(f'  [{U[0,0]:.3f} {U[0,1]:.3f} {U[0,2]:.3f}] [{U[1,0]:.3f} {U[1,1]:.3f} {U[1,2]:.3f}] [{U[2,0]:.3f} {U[2,1]:.3f} {U[2,2]:.3f}]')

df = pd.DataFrame(results)
df.to_csv('pmns_asym_scan_results.csv', index=False)
print()
print('Saved: pmns_asym_scan_results.csv')