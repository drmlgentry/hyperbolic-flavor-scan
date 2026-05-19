import numpy as np
from scipy.linalg import qr, logm
from itertools import permutations, product as iproduct
from scipy.optimize import minimize
import snappy
import time

PMNS = np.array([[0.821,0.550,0.148],[0.357,0.339,0.871],[0.442,0.762,0.471]])
PERMS = list(permutations([0,1,2]))

def fitness_perm(U):
    best = 1e9
    for perm in PERMS:
        best = min(best, float(np.linalg.norm(U[:,perm] - PMNS, 'fro')))
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

def U_from_L(l21, l31, l32):
    L = np.array([[1.0, 0.0, 0.0],
                  [l21, 1.0, 0.0],
                  [l31, l32, 1.0]])
    try:
        U, _ = qr(L)
        for i in range(3):
            if U[i,i] < 0: U[:,i] *= -1
        return np.abs(U)
    except:
        return None

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

def deep_scan_manifold(mi, label):
    mfld = snappy.OrientableClosedCensus[mi]
    mname = mfld.name()
    G = mfld.fundamental_group()
    gens = G.generators()

    axes, words = [], []
    for wlen in [2, 3, 4]:
        for word in generate_words(gens, wlen):
            try:
                mat = to_numpy(G.SL2C(word))
            except:
                continue
            ax = axis_from_matrix(mat)
            if ax is None: continue
            axes.append(ax)
            words.append(word)
            if len(axes) >= 80: break
        if len(axes) >= 80: break

    print(f'{label} ({mname}): {len(axes)} word axes extracted')

    best_f = 1e9
    best_info = None
    na = len(axes)
    t0 = time.time()
    count = 0

    for i in range(na):
        for j in range(i+1, na):
            for k in range(j+1, na):
                d21 = np.dot(axes[j], axes[i])
                d31 = np.dot(axes[k], axes[i])
                d32 = np.dot(axes[k], axes[j])

                # Vectorized: scan scale in [-6,6] x sign simultaneously
                # by just scanning l21,l31,l32 in fine grid around dot products
                def fn(p):
                    U = U_from_L(p[0], p[1], p[2])
                    if U is None: return 1e9
                    return fitness_perm(U)

                # Coarse grid first
                for scale in [0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0]:
                    for s1,s2,s3 in [(1,1,1),(1,1,-1),(1,-1,1),(-1,1,1),
                                     (1,-1,-1),(-1,1,-1),(-1,-1,1),(-1,-1,-1)]:
                        l21 = s1*scale*d21
                        l31 = s2*scale*d31
                        l32 = s3*scale*d32
                        U = U_from_L(l21, l31, l32)
                        if U is None: continue
                        f = fitness_perm(U)
                        if f < best_f:
                            best_f = f
                            best_info = (words[i], words[j], words[k], l21, l31, l32, U.copy())
                            # Refine with Nelder-Mead from this point
                            res = minimize(fn, [l21, l31, l32], method='Nelder-Mead',
                                          options={'xatol':1e-9,'fatol':1e-9,'maxiter':2000})
                            if res.fun < best_f:
                                best_f = res.fun
                                U2 = U_from_L(*res.x)
                                best_info = (words[i], words[j], words[k],
                                             res.x[0], res.x[1], res.x[2], U2.copy())
                            elapsed = time.time() - t0
                            print(f'  [{count:5d} | {elapsed:6.1f}s] NEW BEST {best_f:.6f}  '
                                  f'words={words[i]}/{words[j]}/{words[k]}  '
                                  f'l=({best_info[3]:.3f},{best_info[4]:.3f},{best_info[5]:.3f})')
                count += 1

    return best_f, best_info

print('Deep Borel Scan: m003 and m007')
print('='*60)
print('Permutation-aware fitness, Nelder-Mead refinement on every new best')
print()

for mi, label in [(1,'m003'), (2,'m007')]:
    print(f'Scanning {label}...')
    t0 = time.time()
    best_f, info = deep_scan_manifold(mi, label)
    elapsed = time.time() - t0
    w1,w2,w3,l21,l31,l32,U = info
    print()
    print(f'FINAL {label}: fitness={best_f:.6f}  ({elapsed:.0f}s)')
    print(f'  Words: {w1}/{w2}/{w3}')
    print(f'  L entries: l21={l21:.6f}  l31={l31:.6f}  l32={l32:.6f}')
    for perm in PERMS:
        Up = U[:,perm]
        f = float(np.linalg.norm(Up - PMNS, 'fro'))
        if abs(f - best_f) < 1e-6:
            print(f'  Best perm: {perm}')
            for row, trow in zip(Up, PMNS):
                print(f'    [{row[0]:.4f} {row[1]:.4f} {row[2]:.4f}]  target [{trow[0]:.3f} {trow[1]:.3f} {trow[2]:.3f}]')
            break
    print()

print('Done.')