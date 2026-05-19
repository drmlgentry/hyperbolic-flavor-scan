import snappy
import numpy as np
from scipy.linalg import logm, qr
from itertools import product as iproduct, combinations

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

def borel_fitness(axes, target):
    from itertools import product as ip
    best = 999.
    best_params = None
    for lam in np.linspace(0.1, 5.0, 30):
        for s in ip([1,-1], repeat=3):
            s21,s31,s32 = s
            l21 = lam*s21*np.dot(axes[0],axes[1])
            l31 = lam*s31*np.dot(axes[0],axes[2])
            l32 = lam*s32*np.dot(axes[1],axes[2])
            Lm = np.array([[1.,0.,0.],[l21,1.,0.],[l31,l32,1.]])
            Q,_ = qr(Lm)
            Qabs = np.abs(Q)
            for perm in [[0,1,2],[0,2,1],[1,0,2],
                         [1,2,0],[2,0,1],[2,1,0]]:
                f = np.linalg.norm(Qabs[:,perm]-target,'fro')
                if f < best:
                    best = f
                    best_params = (lam,s,l21,l31,l32,perm)
    return best, best_params

from scipy.optimize import minimize

def borel_fitness_continuous(axes, target):
    """Continuous optimization over l21, l31, l32."""
    d12 = np.dot(axes[0],axes[1])
    d13 = np.dot(axes[0],axes[2])
    d23 = np.dot(axes[1],axes[2])

    def f(params):
        l21,l31,l32 = params
        Lm = np.array([[1.,0.,0.],[l21,1.,0.],[l31,l32,1.]])
        Q,_ = qr(Lm)
        Qabs = np.abs(Q)
        return min(np.linalg.norm(Qabs[:,p]-target,'fro')
                   for p in [[0,1,2],[0,2,1],[1,0,2],
                              [1,2,0],[2,0,1],[2,1,0]])

    best = 999.
    best_x = None
    for lam in [0.5,1.0,1.5,2.0,3.0]:
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

PMNS = np.array([[0.821,0.550,0.148],
                 [0.357,0.339,0.871],
                 [0.442,0.762,0.471]])

letters = ['a','b','A','B']

manifolds = [
    (1,   "m003 base (PMNS)",      "Z/5"),
    (39,  "deg-2 cover idx39",     "Z/55"),
    (238, "deg-3 cover idx238",    "Z/80"),
]

print("=== Borel construction on PMNS covering tower ===")
print("Theoretical minimum: 0.01897")
print()

for idx, label, h1 in manifolds:
    M = snappy.OrientableClosedCensus[idx]
    rho = M.polished_holonomy()
    vol = float(M.volume())
    print(f"--- {label} (idx={idx}, vol={vol:.4f}, H1={h1}) ---")

    # Collect hyperbolic words up to length 4
    hyp = []
    for Lw in range(1,5):
        for w in [''.join(p) for p in iproduct(letters,repeat=Lw)]:
            try:
                mat = np.array(rho(w), dtype=complex)
                if abs(float(abs(np.trace(mat)))) > 2.01:
                    n = get_axis(rho, w)
                    if n is not None:
                        hyp.append((w,n))
            except: pass

    # Deduplicate by axis
    unique = []
    seen = []
    for w,n in hyp:
        if not any(abs(np.dot(n,n2))>0.9999 for n2 in seen):
            unique.append((w,n))
            seen.append(n)

    print(f"  Unique axes: {len(unique)}")

    # Scan triples
    best_f = 999.
    best_t = None
    best_x = None

    for (w1,n1),(w2,n2),(w3,n3) in combinations(unique[:40],3):
        axes = [n1,n2,n3]
        f, x = borel_fitness_continuous(axes, PMNS)
        if f < best_f:
            best_f = f
            best_t = (w1,w2,w3)
            best_x = x

    print(f"  Best triple: {best_t}")
    print(f"  Fitness: {best_f:.5f}  (base m003 target: 0.01897)")

    if best_x is not None and best_t is not None:
        l21,l31,l32 = best_x
        axes = [dict(unique)[w] for w in best_t]
        Lm = np.array([[1.,0.,0.],[l21,1.,0.],[l31,l32,1.]])
        Q,_ = qr(Lm)
        Qabs = np.abs(Q)
        best_perm = min([[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]],
                        key=lambda p: np.linalg.norm(Qabs[:,p]-PMNS,'fro'))
        print(f"  |U_geom|:")
        for row in Qabs[:,best_perm]:
            print("    " + "  ".join(f"{x:.4f}" for x in row))
        print(f"  l21={l21:.4f}  l31={l31:.4f}  l32={l32:.4f}")
    print()

print("PDG PMNS:")
for row in PMNS:
    print("  " + "  ".join(f"{x:.4f}" for x in row))
