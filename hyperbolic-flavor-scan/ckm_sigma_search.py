import snappy
import numpy as np
import itertools
from scipy.linalg import logm as scipy_logm

phi = (1+np.sqrt(5))/2
log_phi = np.log(phi)

def axis_imlogm(A):
    try:
        L = scipy_logm(np.array(A, dtype=complex))
        imL = L.imag
        nx = float((imL[0,1] + imL[1,0]).real)
        ny = float((1j*(imL[1,0] - imL[0,1])).real)
        nz = float((imL[0,0] - imL[1,1]).real)
        n = np.array([nx, ny, nz])
        nrm = np.linalg.norm(n)
        return n/nrm if nrm > 1e-10 else None
    except:
        return None

def angle_between(n1, n2):
    return np.arccos(np.clip(abs(np.dot(n1,n2)), 0, 1))

def build_kernel(axes, sigma):
    N = len(axes)
    K = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            th = angle_between(axes[i], axes[j])
            K[i,j] = np.exp(-th*th/(2*sigma*sigma))
    return K

def eff_rank(K):
    ev = np.maximum(np.linalg.eigvalsh(K), 1e-12)
    p = ev/ev.sum()
    return np.exp(-np.sum(p*np.log(p)))

def best_fit(axes, exp_matrix):
    best = 1e9
    best_s = 0
    best_U = None
    for sigma100 in range(10, 151):
        sigma = sigma100/100.0
        K = build_kernel(axes, sigma)
        Q, R = np.linalg.qr(K)
        for i in range(3):
            if R[i,i] < 0:
                Q[:,i] = -Q[:,i]
        U = np.abs(Q)
        f = np.sqrt(np.mean((U - exp_matrix)**2))
        if f < best:
            best = f
            best_s = sigma
            best_U = U.copy()
    return best, best_s, best_U

CKM_exp = np.array([
    [0.974, 0.225, 0.004],
    [0.225, 0.973, 0.041],
    [0.009, 0.040, 0.999]
])

PMNS_exp = np.array([
    [0.821, 0.551, 0.150],
    [0.358, 0.694, 0.624],
    [0.444, 0.461, 0.767]
])

print("="*60, flush=True)
print("CKM: searching all word triples up to length 3", flush=True)
print("="*60, flush=True)

M = snappy.Manifold('m006')
M.dehn_fill((-5,2))
G = M.fundamental_group(simplify_presentation=True)
gens = G.generators()
PH = M.polished_holonomy()
letters = gens + [g.upper() for g in gens]

all_axes = []
all_words = []
for length in [1,2,3]:
    for wt in itertools.product(letters, repeat=length):
        word = ''.join(wt)
        try:
            A = np.array(PH.SL2C(word), dtype=complex)
            tr = abs(np.trace(A))
            if tr <= 2.001:
                continue
            n = axis_imlogm(A)
            if n is None or np.any(np.isnan(n)):
                continue
            if not any(abs(np.dot(n,m)) > 0.999 for m in all_axes):
                all_axes.append(n)
                all_words.append(word)
        except:
            pass

N = len(all_axes)
print("Total unique axes: %d" % N, flush=True)

best_global = 1e9
best_info = None

for idx in itertools.combinations(range(N), 3):
    axes3 = np.array([all_axes[i] for i in idx])
    words3 = tuple(all_words[i] for i in idx)
    f, s, U = best_fit(axes3, CKM_exp)
    if f < best_global:
        best_global = f
        best_info = (words3, s, U, axes3.copy())
    if f < 0.05:
        print("  GOOD: %s fit=%.6f sigma=%.4f" % (str(words3), f, s),
              flush=True)

words3, sigma_best, U_best, axes_best = best_info
print("", flush=True)
print("Best CKM triple: %s" % str(words3), flush=True)
print("Best fitness: %.6f at sigma=%.4f" % (best_global, sigma_best),
      flush=True)
print("log(phi)=%.4f  (3/2)log(phi)=%.4f" % (log_phi, 1.5*log_phi),
      flush=True)
print("Axes:", flush=True)
for i,(w,n) in enumerate(zip(words3, axes_best)):
    print("  [%d] '%s': %s" % (i, w, str(np.round(n,4))), flush=True)
print("Pairwise angles:", flush=True)
for i in range(3):
    for j in range(i+1,3):
        th = angle_between(axes_best[i], axes_best[j])
        print("  %d vs %d: %.2f deg" % (i,j,np.degrees(th)), flush=True)
print("Reconstructed:", flush=True)
for row in U_best:
    print("  %s" % str(np.round(row,4)), flush=True)
K_best = build_kernel(axes_best, sigma_best)
print("R_eff=%.5f" % eff_rank(K_best), flush=True)

print("", flush=True)
print("="*60, flush=True)
print("PMNS: same search", flush=True)
print("="*60, flush=True)

M2 = snappy.Manifold('m003')
M2.dehn_fill((-2,3))
G2 = M2.fundamental_group(simplify_presentation=True)
gens2 = G2.generators()
PH2 = M2.polished_holonomy()
letters2 = gens2 + [g.upper() for g in gens2]

all_axes2 = []
all_words2 = []
for length in [1,2,3]:
    for wt in itertools.product(letters2, repeat=length):
        word = ''.join(wt)
        try:
            A = np.array(PH2.SL2C(word), dtype=complex)
            tr = abs(np.trace(A))
            if tr <= 2.001:
                continue
            n = axis_imlogm(A)
            if n is None or np.any(np.isnan(n)):
                continue
            if not any(abs(np.dot(n,m)) > 0.999 for m in all_axes2):
                all_axes2.append(n)
                all_words2.append(word)
        except:
            pass

N2 = len(all_axes2)
print("Total unique axes: %d" % N2, flush=True)

best_global2 = 1e9
best_info2 = None
for idx in itertools.combinations(range(N2), 3):
    axes3 = np.array([all_axes2[i] for i in idx])
    words3 = tuple(all_words2[i] for i in idx)
    f, s, U = best_fit(axes3, PMNS_exp)
    if f < best_global2:
        best_global2 = f
        best_info2 = (words3, s, U, axes3.copy())
    if f < 0.05:
        print("  GOOD: %s fit=%.6f sigma=%.4f" % (str(words3), f, s),
              flush=True)

words3, sigma_best, U_best, axes_best = best_info2
print("Best PMNS triple: %s" % str(words3), flush=True)
print("Best fitness: %.6f at sigma=%.4f" % (best_global2, sigma_best),
      flush=True)
print("Axes:", flush=True)
for i,(w,n) in enumerate(zip(words3, axes_best)):
    print("  [%d] '%s': %s" % (i, w, str(np.round(n,4))), flush=True)
print("Reconstructed:", flush=True)
for row in U_best:
    print("  %s" % str(np.round(row,4)), flush=True)
K_best = build_kernel(axes_best, sigma_best)
print("R_eff=%.5f" % eff_rank(K_best), flush=True)
ev = np.sort(np.linalg.eigvalsh(K_best))[::-1]
print("Eigenvalues: %s" % str(np.round(ev/ev[0],5)), flush=True)
