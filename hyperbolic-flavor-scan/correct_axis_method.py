import snappy, numpy as np, itertools
from scipy.linalg import logm as scipy_logm

phi = (1+np.sqrt(5))/2
log_phi = np.log(phi)

def axis_from_holonomy_logm(A):
    """
    Correct method: axis from imaginary part of log of SL2C matrix.
    Decomposes into rotation + translation; axis is from rotation part.
    Maps to S^2 via Bloch sphere / Pauli decomposition.
    """
    try:
        L = scipy_logm(np.array(A, dtype=complex))
        # Imaginary part gives the rotation generator
        # Pauli decomposition of im(L):
        imL = L.imag
        # sigma_x, sigma_y, sigma_z components
        nx = imL[0,1] + imL[1,0]  # Tr(sigma_x * imL)
        ny = 1j*(imL[0,1] - imL[1,0])  # Tr(sigma_y * imL)
        nz = imL[0,0] - imL[1,1]  # Tr(sigma_z * imL)
        # ny is real for SL2C
        n = np.array([nx.real, ny.real if hasattr(ny,'real') else float(ny), nz.real])
        nrm = np.linalg.norm(n)
        if nrm < 1e-10: return None
        return n / nrm
    except:
        return None

def axis_from_real_logm(A):
    """Alternative: real part of log."""
    try:
        L = scipy_logm(np.array(A, dtype=complex))
        reL = L.real
        nx = reL[0,1] + reL[1,0]
        ny = reL[0,1] - reL[1,0]
        nz = reL[0,0] - reL[1,1]
        n = np.array([nx, ny, nz])
        nrm = np.linalg.norm(n)
        if nrm < 1e-10: return None
        return n / nrm
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

print("Testing axis extraction methods on known CKM triple", flush=True)
print("Known best: aaB, AbA, AAb  fitness~0.017", flush=True)
print()

M = snappy.Manifold('m006')
M.dehn_fill((-5,2))
PH = M.polished_holonomy()

known_words = ['aaB','AbA','AAb']

for method_name, method_fn in [
    ('imaginary logm', axis_from_holonomy_logm),
    ('real logm',      axis_from_real_logm),
]:
    print(f"--- Method: {method_name} ---", flush=True)
    axes = []
    for w in known_words:
        try:
            A = np.array(PH.SL2C(w), dtype=complex)
            n = method_fn(A)
            if n is not None:
                axes.append(n)
                tr = abs(np.trace(A))
                print(f"  '{w}': tr={tr:.4f}  n={np.round(n,4)}", flush=True)
        except Exception as e:
            print(f"  '{w}': error {e}", flush=True)

    if len(axes) == 3:
        axes = np.array(axes)
        K = build_kernel(axes, 1.5*log_phi)
        Q, R = np.linalg.qr(K)
        for i in range(3):
            if R[i,i] < 0: Q[:,i] = -Q[:,i]
        U = np.abs(Q)
        fitness = np.sqrt(np.mean((U - CKM_exp)**2))
        print(f"  Fitness vs CKM: {fitness:.6f}", flush=True)
        print(f"  R_eff at (3/2)log(phi): {eff_rank(K):.5f}", flush=True)
        print(f"  Reconstructed:", flush=True)
        for row in U: print(f"    {np.round(row,4)}", flush=True)
    print()

# Also try the known PMNS best triple
print("--- PMNS: known words (if any from earlier scan) ---", flush=True)
M2 = snappy.Manifold('m003')
M2.dehn_fill((-2,3))
PH2 = M2.polished_holonomy()
# From session memory: PMNS best was found with word_triple_scan_corrected.py
# Try some candidates
for words in [['a','b','ab'],['b','ab','ba'],['a','ab','aba']]:
    axes = []
    for w in words:
        try:
            A = np.array(PH2.SL2C(w), dtype=complex)
            n = axis_from_holonomy_logm(A)
            if n is not None: axes.append(n)
        except: pass
    if len(axes)==3:
        axes = np.array(axes)
        K = build_kernel(axes, log_phi)
        Q,R = np.linalg.qr(K)
        for i in range(3):
            if R[i,i]<0: Q[:,i]=-Q[:,i]
        U = np.abs(Q)
        fit = np.sqrt(np.mean((U-PMNS_exp)**2))
        print(f"  {words}: fitness={fit:.6f}", flush=True)
