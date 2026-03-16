import numpy as np
from scipy.linalg import qr, logm
import snappy
import time

PMNS = np.array([[0.821,0.550,0.148],[0.357,0.339,0.871],[0.442,0.762,0.471]])

def to_numpy(m):
    return np.array([[complex(m[0,0]), complex(m[0,1])],
                     [complex(m[1,0]), complex(m[1,1])]])

def fixed_points(M):
    # Returns (z+, z-) attracting/repelling fixed points on CP1
    a,b,c,d = M[0,0], M[0,1], M[1,0], M[1,1]
    if abs(c) < 1e-10:
        return (b/(a-d) if abs(a-d)>1e-10 else 1e10+0j), np.inf
    disc = np.sqrt((d-a)**2 + 4*b*c + 0j)
    z1 = (-(d-a) + disc) / (2*c)
    z2 = (-(d-a) - disc) / (2*c)
    # Attracting = smaller |eigenvalue|
    lam1 = a + c*z1
    lam2 = a + c*z2
    if abs(lam1) < abs(lam2):
        return z1, z2
    return z2, z1

def stereo(z):
    if abs(z) > 1e8: return np.array([0.,0.,1.])
    x,y = z.real, z.imag
    r2 = x**2 + y**2
    return np.array([2*x/(1+r2), 2*y/(1+r2), (r2-1)/(r2+1)])

def cp1_inner(z1, z2):
    # Fubini-Study inner product on CP1
    # |<z1|z2>|^2 = |1 + conj(z1)*z2|^2 / ((1+|z1|^2)(1+|z2|^2))
    if abs(z1) > 1e8 or abs(z2) > 1e8:
        return 0.0 + 0j
    num = 1 + np.conj(z1)*z2
    den = np.sqrt((1+abs(z1)**2)*(1+abs(z2)**2))
    return num / den  # complex — keeps phase

def hermitian_overlap(zplus_list, zminus_list, sigma):
    # O_ij = <z+_i | z+_j> * exp(-theta_ij^2/2s^2)
    # Hermitian: O_ij = conj(O_ji)  =>  |U| from complex QR is genuinely asymmetric
    n = len(zplus_list)
    O = np.eye(n, dtype=complex)
    vecs = [stereo(z) for z in zplus_list]
    for i in range(n):
        for j in range(i+1, n):
            cos_ang = np.clip(np.dot(vecs[i], vecs[j]), -1, 1)
            theta = np.arccos(abs(cos_ang))
            gauss = np.exp(-theta**2 / (2*sigma**2))
            phase = cp1_inner(zplus_list[i], zplus_list[j])
            val = gauss * phase
            O[i,j] = val
            O[j,i] = np.conj(val)
    return O

def hermitian_overlap_mixed(zplus_list, zminus_list, sigma):
    # O_ij = <z+_i | z-_j>  -- cross overlap between attracting/repelling
    # Fully asymmetric: O_ij != conj(O_ji) in general
    n = len(zplus_list)
    O = np.eye(n, dtype=complex)
    for i in range(n):
        for j in range(n):
            if i == j: continue
            val = cp1_inner(zplus_list[i], zminus_list[j])
            vecs_i = stereo(zplus_list[i])
            vecs_j = stereo(zminus_list[j])
            cos_ang = np.clip(np.dot(vecs_i, vecs_j), -1, 1)
            theta = np.arccos(abs(cos_ang))
            gauss = np.exp(-theta**2 / (2*sigma**2))
            O[i,j] = gauss * val
    return O

def U_from_complex_O(O):
    # QR of complex matrix -- |U| is generically asymmetric
    try:
        if not np.all(np.linalg.eigvalsh(O) > 1e-10):
            return None
    except:
        return None
    U, _ = qr(O)
    for i in range(3):
        if U[i,i].real < 0: U[:,i] *= -1
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

print('PMNS Complex QR / Fixed-Point Scan — Prototype')
print('='*60)
print('Testing Hermitian and mixed fixed-point overlap matrices')
print('Benchmark: 0.018968')
print()

sigmas = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0]
census = snappy.OrientableClosedCensus
N = 50
results = []
t0 = time.time()

for mi in range(N):
    mfld = census[mi]
    mname = mfld.name()
    try:
        G = mfld.fundamental_group()
        gens = G.generators()
    except:
        continue

    zplus_list, zminus_list, words = [], [], []

    for wlen in [2, 3]:
        for word in generate_words(gens, wlen):
            try:
                mat = to_numpy(G.SL2C(word))
            except:
                continue
            zp, zm = fixed_points(mat)
            if abs(zp) > 1e7 or abs(zm) > 1e7: continue
            zplus_list.append(zp)
            zminus_list.append(zm)
            words.append(word)
            if len(words) >= 50: break
        if len(words) >= 50: break

    if len(words) < 3: continue

    best_f = 1e9
    best_row = None
    nw = len(words)

    for i in range(nw):
        for j in range(i+1, nw):
            for k in range(j+1, nw):
                zp3 = [zplus_list[i], zplus_list[j], zplus_list[k]]
                zm3 = [zminus_list[i], zminus_list[j], zminus_list[k]]
                for sigma in sigmas:
                    for name, fn in [('herm', hermitian_overlap),
                                     ('mixed', hermitian_overlap_mixed)]:
                        O = fn(zp3, zm3, sigma)
                        U = U_from_complex_O(O)
                        if U is None: continue
                        f = fitness(U)
                        if f < best_f:
                            best_f = f
                            best_row = {
                                'manifold': mname,
                                'method': name,
                                'words': f'{words[i]}/{words[j]}/{words[k]}',
                                'sigma': sigma,
                                'fitness': round(f,6),
                                'U': U.copy()
                            }

    if best_row:
        results.append(best_row)
        if best_f < 0.10:
            elapsed = time.time() - t0
            U = best_row['U']
            print(f'[{mi:3d} | {elapsed:5.1f}s] {mname} {best_row["method"]}  fit={best_f:.6f}  s={best_row["sigma"]}  {best_row["words"]}')
            print(f'  [{U[0,0]:.3f} {U[0,1]:.3f} {U[0,2]:.3f}]  target: [0.821 0.550 0.148]')
            print(f'  [{U[1,0]:.3f} {U[1,1]:.3f} {U[1,2]:.3f}]  target: [0.357 0.339 0.871]')
            print(f'  [{U[2,0]:.3f} {U[2,1]:.3f} {U[2,2]:.3f}]  target: [0.442 0.762 0.471]')
            print()

    if (mi+1) % 10 == 0:
        elapsed = time.time() - t0
        top = sorted(results, key=lambda x: x['fitness'])[:3]
        print(f'--- [{mi+1:3d}/{N} | {elapsed:.0f}s] top3: ' +
              ' | '.join(f'{r["manifold"]}({r["method"]})={r["fitness"]:.4f}' for r in top) + ' ---')

print()
print('='*60)
print('TOP 15')
print('='*60)
for r in sorted(results, key=lambda x: x['fitness'])[:15]:
    U = r['U']
    print(f'{r["manifold"]:10s} {r["method"]:6s}  fit={r["fitness"]:.6f}  s={r["sigma"]}  {r["words"]}')
    print(f'  [{U[0,0]:.3f} {U[0,1]:.3f} {U[0,2]:.3f}] [{U[1,0]:.3f} {U[1,1]:.3f} {U[1,2]:.3f}] [{U[2,0]:.3f} {U[2,1]:.3f} {U[2,2]:.3f}]')