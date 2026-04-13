import snappy
import numpy as np
from scipy.linalg import logm, qr
from itertools import product as iproduct, combinations
from scipy.optimize import minimize

def get_axis(rho, word):
    try:
        mat = np.array(rho(word), dtype=complex)
        mat = mat / np.sqrt(np.linalg.det(mat))
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

PMNS = np.array([[0.821,0.550,0.148],
                 [0.357,0.339,0.871],
                 [0.442,0.762,0.471]])

def borel_fit_continuous(n1, n2, n3):
    """Continuous optimization over lambda and signs."""
    d12 = np.dot(n1,n2)
    d13 = np.dot(n1,n3)
    d23 = np.dot(n2,n3)

    def f(params):
        l21,l31,l32 = params
        Lm = np.array([[1.,0.,0.],[l21,1.,0.],[l31,l32,1.]])
        Q,_ = qr(Lm)
        Qabs = np.abs(Q)
        return min(np.linalg.norm(Qabs[:,p]-PMNS,'fro')
                   for p in [[0,1,2],[0,2,1],[1,0,2],
                              [1,2,0],[2,0,1],[2,1,0]])

    best = 999.
    best_x = None
    # Try multiple starting points based on dot products
    for lam in [0.5, 1.0, 1.5, 2.0]:
        for s21 in [1,-1]:
            for s31 in [1,-1]:
                for s32 in [1,-1]:
                    x0 = [lam*s21*d12, lam*s31*d13, lam*s32*d23]
                    res = minimize(f, x0, method='Nelder-Mead',
                                   options={'xatol':1e-8,'fatol':1e-8,
                                            'maxiter':5000})
                    if res.fun < best:
                        best = res.fun
                        best_x = res.x
    return best, best_x

M = snappy.OrientableClosedCensus[1]
rho = M.polished_holonomy()
letters = ['a','b','A','B']

# Collect hyperbolic words up to length 4
hyp = []
for L in range(1,5):
    for w in [''.join(p) for p in iproduct(letters,repeat=L)]:
        try:
            mat = np.array(rho(w), dtype=complex)
            if abs(float(abs(np.trace(mat)))) > 2.01:
                n = get_axis(rho, w)
                if n is not None:
                    hyp.append((w,n))
        except: pass

# Deduplicate
unique = []
seen = []
for w,n in hyp:
    if not any(abs(np.dot(n,n2))>0.9999 for n2 in seen):
        unique.append((w,n))
        seen.append(n)

print(f"Unique axes: {len(unique)}")
print()

# Refine best triple with continuous optimization
print("=== Refining best triple with continuous optimization ===")
best_triple_words = ('aBAb', 'baba', 'bABa')
axes = [get_axis(rho, w) for w in best_triple_words]
fit, x = borel_fit_continuous(*axes)
l21,l31,l32 = x
print(f"Triple: {best_triple_words}")
print(f"Fitness (continuous): {fit:.5f}")
print(f"l21={l21:.4f} l31={l31:.4f} l32={l32:.4f}")
Lm = np.array([[1.,0.,0.],[l21,1.,0.],[l31,l32,1.]])
Q,_ = qr(Lm)
Qabs = np.abs(Q)
best_perm = min([[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]],
                key=lambda p: np.linalg.norm(Qabs[:,p]-PMNS,'fro'))
print(f"|U_geom|:")
for row in Qabs[:,best_perm]:
    print('  '+'  '.join(f'{x:.4f}' for x in row))
print()

# Scan with continuous optimization on top 50 unique axes
print("=== Full scan with continuous optimization ===")
print("(scanning top 50 unique axes, ~19k triples)")
best_f = 999.
best_t = None
best_x = None
count = 0
for (w1,n1),(w2,n2),(w3,n3) in combinations(unique[:50],3):
    f, x = borel_fit_continuous(n1,n2,n3)
    count += 1
    if f < best_f:
        best_f = f
        best_t = (w1,w2,w3)
        best_x = x
    if count % 1000 == 0:
        print(f"  {count} triples... best: {best_f:.5f}")

l21,l31,l32 = best_x
print()
print(f"Best triple: {best_t}")
print(f"Fitness: {best_f:.5f}  (theoretical min: 0.01897)")
print(f"l21={l21:.4f} l31={l31:.4f} l32={l32:.4f}")
print(f"Paper: l21=0.443 l31=-0.530 l32=0.432")
Lm = np.array([[1.,0.,0.],[l21,1.,0.],[l31,l32,1.]])
Q,_ = qr(Lm)
Qabs = np.abs(Q)
best_perm = min([[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]],
                key=lambda p: np.linalg.norm(Qabs[:,p]-PMNS,'fro'))
print(f"|U_geom| (perm {best_perm}):")
for row in Qabs[:,best_perm]:
    print('  '+'  '.join(f'{x:.4f}' for x in row))
print()
print("PDG PMNS:")
for row in PMNS:
    print('  '+'  '.join(f'{x:.4f}' for x in row))
