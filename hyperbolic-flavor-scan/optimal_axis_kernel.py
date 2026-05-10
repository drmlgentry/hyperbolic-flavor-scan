import snappy
import numpy as np
import itertools

phi = (1+np.sqrt(5))/2
log_phi = np.log(phi)

def axis_from_SL2C(A):
    a,b = complex(A[0,0]),complex(A[0,1])
    c,d = complex(A[1,0]),complex(A[1,1])
    if abs(c) < 1e-10:
        if abs(d-a) < 1e-10: return None
        z = -b/(d-a)
    else:
        disc = (d-a)**2 + 4*b*c
        sq = np.sqrt(disc+0j)
        z1 = (-(d-a)+sq)/(2*c)
        z2 = (-(d-a)-sq)/(2*c)
        z = z1 if abs(z1) >= abs(z2) else z2
    if abs(z) > 1e8: return np.array([0.,0.,1.])
    x,y = z.real,z.imag
    r2 = x*x+y*y
    d2 = 1.0+r2
    n = np.array([2*x/d2, 2*y/d2, (r2-1)/d2])
    nrm = np.linalg.norm(n)
    return n/nrm if nrm > 1e-10 else None

def angle_between(n1,n2):
    return np.arccos(np.clip(abs(np.dot(n1,n2)),0,1))

def build_kernel(axes, sigma):
    N = len(axes)
    K = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            th = angle_between(axes[i],axes[j])
            K[i,j] = np.exp(-th*th/(2*sigma*sigma))
    return K

def eff_rank(K):
    ev = np.maximum(np.linalg.eigvalsh(K),1e-12)
    p = ev/ev.sum()
    return np.exp(-np.sum(p*np.log(p)))

def energy_frac(K,n):
    ev = np.sort(np.linalg.eigvalsh(K))[::-1]
    return ev[:n].sum()/ev.sum()

PMNS_exp = np.array([
    [0.821, 0.551, 0.150],
    [0.358, 0.694, 0.624],
    [0.444, 0.461, 0.767]
])

CKM_exp = np.array([
    [0.974, 0.225, 0.004],
    [0.225, 0.973, 0.041],
    [0.009, 0.040, 0.999]
])

print("="*65, flush=True)
print("OPTIMAL 3-AXIS KERNEL: finding best word triple", flush=True)
print("="*65, flush=True)

for label, mname, filling, exp_matrix, sigma_opt in [
    ('PMNS', 'm003', (-2,3), PMNS_exp, log_phi),
    ('CKM',  'm006', (-5,2), CKM_exp,  1.5*log_phi),
]:
    print(f"\n--- {label}: {mname}{filling} ---", flush=True)
    M = snappy.Manifold(mname)
    M.dehn_fill(filling)
    G = M.fundamental_group(simplify_presentation=True)
    gens = G.generators()
    PH = M.polished_holonomy()

    all_axes = []
    all_words = []
    letters = gens + [g.upper() for g in gens]

    for length in [1, 2, 3]:
        for wt in itertools.product(letters, repeat=length):
            word = ''.join(wt)
            try:
                A = np.array(PH.SL2C(word), dtype=complex)
                tr = abs(np.trace(A))
                if tr <= 2.001: continue
                n = axis_from_SL2C(A)
                if n is None or np.any(np.isnan(n)): continue
                if not any(abs(np.dot(n,m)) > 0.9999 for m in all_axes):
                    all_axes.append(n)
                    all_words.append(word)
            except: pass

    N = len(all_axes)
    print(f"  Total unique axes: {N}", flush=True)

    best_fitness = 1e9
    best_triple = None
    best_axes = None

    for idx in itertools.combinations(range(N), 3):
        axes3 = np.array([all_axes[i] for i in idx])
        K3 = build_kernel(axes3, sigma_opt)
        Q, R = np.linalg.qr(K3)
        for i in range(3):
            if R[i,i] < 0: Q[:,i] = -Q[:,i]
        U = np.abs(Q)
        fitness = np.sqrt(np.mean((U - exp_matrix)**2))
        if fitness < best_fitness:
            best_fitness = fitness
            best_triple = tuple(all_words[i] for i in idx)
            best_axes = axes3.copy()

    print(f"  Best triple: {best_triple}  fitness={best_fitness:.6f}",
          flush=True)
    for i,(w,n) in enumerate(zip(best_triple, best_axes)):
        print(f"    [{i}] '{w}': n={np.round(n,4)}", flush=True)

    print(f"  Pairwise angles:", flush=True)
    for i in range(3):
        for j in range(i+1,3):
            th = angle_between(best_axes[i],best_axes[j])
            print(f"    axis_{i} vs axis_{j}: {np.degrees(th):.2f} deg",
                  flush=True)

    print(f"\n  Kernel rank sweep:", flush=True)
    print(f"  {'sigma':>8} {'R_eff':>8} {'E3':>8} {'E2':>8}  note", flush=True)
    print("  "+"-"*48, flush=True)
    for sigma in np.arange(0.20, 1.31, 0.10):
        K = build_kernel(best_axes, sigma)
        R = eff_rank(K)
        E3 = energy_frac(K,3)
        E2 = energy_frac(K,2)
        note = ''
        if abs(sigma-log_phi) < 0.06:       note = '<-- log(phi)'
        elif abs(sigma-1.5*log_phi) < 0.06: note = '<-- (3/2)log(phi)'
        print(f"  {sigma:>8.4f} {R:>8.4f} {E3:>8.4f} {E2:>8.4f}  {note}",
              flush=True)

    K_opt = build_kernel(best_axes, sigma_opt)
    ev = np.sort(np.linalg.eigvalsh(K_opt))[::-1]
    ev_n = ev/ev[0]
    print(f"\n  Eigenvalues at sigma={sigma_opt:.4f}:", flush=True)
    for i in range(3):
        print(f"    lambda_{i}={ev_n[i]:.5f}  gap={ev_n[i]/ev_n[i+1]:.4f}",
              flush=True)

    Q2, R2 = np.linalg.qr(K_opt)
    for i in range(3):
        if R2[i,i] < 0: Q2[:,i] = -Q2[:,i]
    U = np.abs(Q2)
    print(f"\n  Reconstructed |U|:", flush=True)
    print(f"  {np.round(U,4)}", flush=True)
    print(f"  Experimental:", flush=True)
    print(f"  {exp_matrix}", flush=True)
    print(f"  RMS: {np.sqrt(np.mean((U-exp_matrix)**2)):.6f}", flush=True)

print("\nDone.", flush=True)
