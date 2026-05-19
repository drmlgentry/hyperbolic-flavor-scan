import numpy as np
from scipy.linalg import qr, svd
from scipy.optimize import minimize

PMNS = np.array([[0.821,0.550,0.148],[0.357,0.339,0.871],[0.442,0.762,0.471]])

def make_O(o12, o13, o23):
    return np.array([[1.0,o12,o13],[o12,1.0,o23],[o13,o23,1.0]])

def is_pd(O):
    return np.all(np.linalg.eigvalsh(O) > 1e-10)

def U_qr(O):
    U, _ = qr(O)
    for i in range(3):
        if U[i,i] < 0: U[:,i] *= -1
    return np.abs(U)

def U_lq(O):
    U, _ = qr(O.T)
    for i in range(3):
        if U[i,i] < 0: U[:,i] *= -1
    return np.abs(U.T)

def U_svd(O):
    Umat, s, Vt = svd(O)
    W = Umat @ Vt
    for i in range(3):
        if W[i,i] < 0: W[:,i] *= -1
    return np.abs(W)

def U_chol(O):
    try:
        L = np.linalg.cholesky(O)
        norms = np.linalg.norm(L, axis=1, keepdims=True)
        U = L / norms
        for i in range(3):
            if U[i,i] < 0: U[i,:] *= -1
        return np.abs(U)
    except:
        return None

def fit(U, target):
    return float(np.linalg.norm(U - target, 'fro'))

print('='*60)
print('PART 1: Four decompositions vs PMNS')
print('='*60)

np.random.seed(42)
methods = ['QR','LQ','SVD','Chol']
best = {m: 1e9 for m in methods}
best_p = {m: None for m in methods}
fns = {'QR':U_qr,'LQ':U_lq,'SVD':U_svd,'Chol':U_chol}

def make_fitness(name):
    def f(p):
        o12,o13,o23 = p
        if any(x<0 or x>0.9999 for x in p): return 1e9
        O = make_O(o12,o13,o23)
        if not is_pd(O): return 1e9
        U = fns[name](O)
        if U is None: return 1e9
        return fit(U, PMNS)
    return f

for name in methods:
    fn = make_fitness(name)
    b, bp = 1e9, None
    for i in range(1000):
        p0 = np.random.uniform(0.01, 0.99, 3)
        res = minimize(fn, p0, method='Nelder-Mead',
                       options={'xatol':1e-8,'fatol':1e-8,'maxiter':2000})
        if res.fun < b:
            b, bp = res.fun, res.x
            print(name, i, round(b,6), flush=True)
    best[name] = b
    best_p[name] = bp
    O = make_O(*bp)
    U = fns[name](O)
    print(name, 'FINAL', round(b,6), 'overlaps', [round(x,3) for x in bp])
    if U is not None:
        for row in U:
            print('  ', [round(x,3) for x in row])

print()
print('='*60)
print('PART 2: Two-matrix constructions')
print('='*60)

def fitness_avg(p):
    if any(x<0 or x>0.9999 for x in p): return 1e9
    O1 = make_O(p[0],p[1],p[2])
    O2 = make_O(p[3],p[4],p[5])
    if not is_pd(O1) or not is_pd(O2): return 1e9
    O = (O1 + O2) / 2
    if not is_pd(O): return 1e9
    return fit(U_qr(O), PMNS)

def fitness_prod(p):
    if any(x<0 or x>0.9999 for x in p): return 1e9
    O1 = make_O(p[0],p[1],p[2])
    O2 = make_O(p[3],p[4],p[5])
    if not is_pd(O1) or not is_pd(O2): return 1e9
    O = O1 @ O2
    return fit(U_qr(O), PMNS)

def fitness_asym(p):
    # Fully asymmetric 3x3 with positive entries, QR
    O = np.array([[1.0,   p[0], p[1]],
                  [p[2],  1.0,  p[3]],
                  [p[4],  p[5], 1.0 ]])
    if not is_pd(O): return 1e9
    return fit(U_qr(O), PMNS)

for name, fn, ndim in [
    ('Average(O1,O2)', fitness_avg, 6),
    ('Product(O1,O2)', fitness_prod, 6),
    ('Asymmetric O',   fitness_asym, 6),
]:
    b, bp = 1e9, None
    for i in range(500):
        p0 = np.random.uniform(0.01, 0.99, ndim)
        res = minimize(fn, p0, method='Nelder-Mead',
                       options={'xatol':1e-7,'fatol':1e-7,'maxiter':3000})
        if res.fun < b:
            b, bp = res.fun, res.x
            print(name, i, round(b,6), flush=True)
    print(name, 'FINAL', round(b,6))

print()
print('='*60)
print('SUMMARY')
print('='*60)
print('Gaussian floor:  0.300379')
print('CKM best:        0.017290')
for name in methods:
    print(name, 'PMNS floor:', round(best[name],6))