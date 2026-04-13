import snappy
import numpy as np
from scipy.linalg import logm, qr
from itertools import product as iproduct, combinations

def get_axis(rho, word):
    try:
        mat = np.array(rho(word), dtype=complex)
        mat = mat / np.sqrt(np.linalg.det(mat))
        # Check condition before logm
        eigs = np.linalg.eigvals(mat)
        if any(abs(e.real) > 1e6 for e in eigs):
            return None
        L = logm(mat)
        if not np.all(np.isfinite(L)):
            return None
        x = float(np.real(L[0,1]+L[1,0]))/2
        y = float(np.imag(L[1,0]-L[0,1]))/2
        z = float(np.real(L[0,0]-L[1,1]))/2
        v = np.array([x,y,z])
        n = np.linalg.norm(v)
        return v/n if n>1e-10 else None
    except:
        return None

def borel_fitness(n1, n2, n3, target):
    best = 999.
    best_params = None
    for lam in np.linspace(0.1, 5.0, 50):
        for s21 in [1,-1]:
            for s31 in [1,-1]:
                for s32 in [1,-1]:
                    l21 = lam*s21*np.dot(n1,n2)
                    l31 = lam*s31*np.dot(n1,n3)
                    l32 = lam*s32*np.dot(n2,n3)
                    Lm = np.array([[1.,0.,0.],
                                   [l21,1.,0.],
                                   [l31,l32,1.]])
                    Q,_ = qr(Lm)
                    Qabs = np.abs(Q)
                    for perm in [[0,1,2],[0,2,1],[1,0,2],
                                 [1,2,0],[2,0,1],[2,1,0]]:
                        f = np.linalg.norm(Qabs[:,perm]-target,'fro')
                        if f < best:
                            best = f
                            best_params = (lam,(s21,s31,s32),
                                          l21,l31,l32,perm)
    return best, best_params

PMNS = np.array([[0.821,0.550,0.148],
                 [0.357,0.339,0.871],
                 [0.442,0.762,0.471]])

# ── First establish theoretical minimum ──────────────────────
print("=== Theoretical minimum of Borel construction ===")
from scipy.optimize import minimize

def neg_fitness(params):
    l21,l31,l32 = params
    Lm = np.array([[1.,0.,0.],[l21,1.,0.],[l31,l32,1.]])
    Q,_ = qr(Lm)
    Qabs = np.abs(Q)
    return min(np.linalg.norm(Qabs[:,p]-PMNS,'fro')
               for p in [[0,1,2],[0,2,1],[1,0,2],
                         [1,2,0],[2,0,1],[2,1,0]])

best_th = 999.
best_th_params = None
np.random.seed(42)
for _ in range(200):
    x0 = np.random.uniform(-3,3,3)
    res = minimize(neg_fitness, x0, method='Nelder-Mead',
                   options={'xatol':1e-9,'fatol':1e-9,'maxiter':10000})
    if res.fun < best_th:
        best_th = res.fun
        best_th_params = res.x

print(f"Theoretical minimum: {best_th:.5f}")
print(f"Optimal (l21,l31,l32) = {best_th_params}")
print(f"Paper reports: 0.01897")
print()

# ── Scan closed OrientableClosedCensus[1] ────────────────────
print("=== Scanning closed m003 (OrientableClosedCensus[1]) ===")
M = snappy.OrientableClosedCensus[1]
rho = M.polished_holonomy()
G = M.fundamental_group()
print(f"vol={float(M.volume()):.4f} H1={M.homology()}")
print(f"Generators: {G.generators()}")
print()

letters = ['a','b','A','B']
hyp = []
for L in range(1,5):  # up to length 4
    for w in [''.join(p) for p in iproduct(letters,repeat=L)]:
        try:
            mat = np.array(rho(w), dtype=complex)
            tr = abs(float(abs(np.trace(mat))))
            if tr > 2.01:
                n = get_axis(rho, w)
                if n is not None:
                    hyp.append((w,n))
        except: pass
    if L == 2:
        print(f"  Hyperbolic words len<={L}: {len(hyp)}")

print(f"  Hyperbolic words len<=4: {len(hyp)}")
print()

# Deduplicate by axis direction
unique = []
seen_axes = []
for w, n in hyp:
    is_dup = False
    for n2 in seen_axes:
        if abs(np.dot(n,n2)) > 0.9999:
            is_dup = True
            break
    if not is_dup:
        unique.append((w,n))
        seen_axes.append(n)

print(f"  Unique axis directions: {len(unique)}")
print()

# Scan triples
print("Scanning triples (this may take a minute)...")
best_f = 999.
best_triple = None
best_U = None
count = 0
for (w1,n1),(w2,n2),(w3,n3) in combinations(unique[:50], 3):
    f, params = borel_fitness(n1, n2, n3, PMNS)
    count += 1
    if f < best_f:
        best_f = f
        best_triple = (w1,w2,w3)
        best_params_found = params
    if count % 500 == 0:
        print(f"  {count} triples... best so far: {best_f:.5f}")

lam,s,l21,l31,l32,perm = best_params_found
print()
print(f"Best triple: {best_triple}")
print(f"Fitness: {best_f:.5f}")
print(f"lambda={lam:.3f} signs={s}")
print(f"l21={l21:.4f} l31={l31:.4f} l32={l32:.4f}")
print()

# Show matrix
Lm = np.array([[1.,0.,0.],[l21,1.,0.],[l31,l32,1.]])
Q,_ = qr(Lm)
print(f"|U_geom| (col perm {perm}):")
for row in np.abs(Q)[:,perm]:
    print('  '+'  '.join(f'{x:.4f}' for x in row))
print()
print("PDG PMNS:")
for row in PMNS:
    print('  '+'  '.join(f'{x:.4f}' for x in row))
print()
print(f"Theoretical min: {best_th:.5f}")
print(f"Achieved:        {best_f:.5f}")
print(f"Ratio:           {best_f/best_th:.4f}")
