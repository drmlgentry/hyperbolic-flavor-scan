from snappy import *
import numpy as np

# ── Helper functions ──────────────────────────────────────────────────
def get_SU2(M, word):
    """Get SU(2) matrix from holonomy, restricted to compact factor K."""
    G = M.fundamental_group()
    mat = G.SL2C(word)
    # Convert to numpy
    A = np.array([[complex(mat[0][0]), complex(mat[0][1])],
                  [complex(mat[1][0]), complex(mat[1][1])]])
    # Project to SU(2): A -> A / sqrt(det(A))
    d = np.linalg.det(A)
    return A / np.sqrt(d)

def bloch_vector(U):
    """Extract Bloch sphere axis from SU(2) matrix U = exp(i*theta*n.sigma)."""
    # U = cos(theta)*I + i*sin(theta)*(n.sigma)
    # Tr(U) = 2*cos(theta)
    tr = np.trace(U)
    cos_theta = tr.real / 2
    cos_theta = np.clip(cos_theta, -1, 1)
    theta = np.arccos(cos_theta)
    if abs(np.sin(theta)) < 1e-10:
        return np.array([0., 0., 1.]), 0.0
    # Extract axis from off-diagonal
    nx = U[0,1].imag / np.sin(theta)
    ny = U[0,1].real / np.sin(theta)
    nz = U[0,0].imag / np.sin(theta)
    n = np.array([nx, ny, nz])
    # Normalize
    n = n / (np.linalg.norm(n) + 1e-15)
    return n, float(theta)

def qubit_state(n, theta):
    """Convert Bloch vector to qubit state |psi> = cos(t/2)|0> + e^{i*phi}*sin(t/2)|1>."""
    # n = (sin(polar)*cos(az), sin(polar)*sin(az), cos(polar))
    polar = np.arccos(np.clip(n[2], -1, 1))
    az    = np.arctan2(n[1], n[0])
    psi   = np.array([np.cos(polar/2),
                      np.exp(1j*az)*np.sin(polar/2)])
    return psi

def inner_product(psi1, psi2):
    """Inner product <psi1|psi2>."""
    return np.dot(psi1.conj(), psi2)

def fidelity(psi1, psi2):
    """Fidelity |<psi1|psi2>|^2."""
    return abs(inner_product(psi1, psi2))**2

# ── CKM manifold: m006 = OrientableClosedCensus[43] ──────────────────
print("=== CKM sector: m006 (OrientableClosedCensus[43]) ===")
M_CKM = OrientableClosedCensus[43]
print(f"Vol={float(M_CKM.volume()):.4f}, H1={M_CKM.homology()}")

ckm_words = ['aaB', 'AbA', 'AAb']
ckm_axes  = []
ckm_states = []

for w in ckm_words:
    try:
        U = get_SU2(M_CKM, w)
        n, theta = bloch_vector(U)
        psi = qubit_state(n, theta)
        ckm_axes.append(n)
        ckm_states.append(psi)
        print(f"\nWord: {w}")
        print(f"  SU(2) rotation angle: {np.degrees(theta):.4f} deg")
        print(f"  Bloch axis: ({n[0]:.4f}, {n[1]:.4f}, {n[2]:.4f})")
        print(f"  Qubit state: ({psi[0]:.4f}, {psi[1]:.4f})")
    except Exception as e:
        print(f"  {w}: {e}")

print("\n--- CKM Fidelity matrix |<psi_i|psi_j>|^2 ---")
print("  (= cos^2(half-angle between Bloch vectors)")
print("  (compare to |U_CKM|^2 entries)")
for i in range(3):
    row = []
    for j in range(3):
        row.append(f"{fidelity(ckm_states[i], ckm_states[j]):.4f}")
    print(f"  {' '.join(row)}")

# Z/5 phase gate
print("\n--- Z/5 phase structure ---")
phase = np.exp(2j * np.pi / 5)
print(f"e^(2pi*i/5) = {phase:.6f}")
print(f"arg/pi      = {np.angle(phase)/np.pi:.6f}")
print(f"This is the T-gate analog for Z/5 topology")

# ── PMNS manifold: m003 = OrientableClosedCensus[1] ──────────────────
print("\n=== PMNS sector: m003 (OrientableClosedCensus[1]) ===")
M_PMNS = OrientableClosedCensus[1]
print(f"Vol={float(M_PMNS.volume()):.4f}, H1={M_PMNS.homology()}")

pmns_words = ['aa', 'ab', 'aB']
pmns_axes  = []
pmns_states = []

for w in pmns_words:
    try:
        U = get_SU2(M_PMNS, w)
        n, theta = bloch_vector(U)
        psi = qubit_state(n, theta)
        pmns_axes.append(n)
        pmns_states.append(psi)
        print(f"\nWord: {w}")
        print(f"  SU(2) rotation angle: {np.degrees(theta):.4f} deg")
        print(f"  Bloch axis: ({n[0]:.4f}, {n[1]:.4f}, {n[2]:.4f})")
        print(f"  Qubit state: ({psi[0]:.4f}, {psi[1]:.4f})")
    except Exception as e:
        print(f"  {w}: {e}")

print("\n--- PMNS Fidelity matrix |<psi_i|psi_j>|^2 ---")
for i in range(3):
    row = []
    for j in range(3):
        row.append(f"{fidelity(pmns_states[i], pmns_states[j]):.4f}")
    print(f"  {' '.join(row)}")

# ── Cross-sector entanglement ─────────────────────────────────────────
print("\n=== Cross-sector fidelities (CKM vs PMNS axes) ===")
print("Measures how 'orthogonal' the two sectors are as qubit systems")
for i, w1 in enumerate(ckm_words):
    for j, w2 in enumerate(pmns_words):
        f = fidelity(ckm_states[i], pmns_states[j])
        print(f"  CKM[{w1}] vs PMNS[{w2}]: {f:.4f}")

# ── Entanglement entropy ──────────────────────────────────────────────
print("\n=== Entanglement structure ===")
print("Von Neumann entropy of each qubit state (should be 0 for pure states)")
print("Concurrence between pairs:")
def concurrence(psi1, psi2):
    """Concurrence for two-qubit pure state |psi1>|psi2>."""
    # For product states this is 0; for entangled it's > 0
    # Here we treat the pair as a 2-qubit system
    rho = np.outer(np.kron(psi1, psi2), np.kron(psi1, psi2).conj())
    # Partial trace over second qubit
    rho_A = rho[:2,:2] + rho[2:,2:]  # simplified
    evals = np.linalg.eigvalsh(rho_A)
    evals = np.clip(evals, 0, 1)
    S = -sum(e*np.log2(e+1e-15) for e in evals)
    return float(S)

print("\nCKM pairs (entanglement entropy of reduced state):")
for i in range(3):
    for j in range(i+1, 3):
        S = concurrence(ckm_states[i], ckm_states[j])
        print(f"  {ckm_words[i]}-{ckm_words[j]}: S={S:.4f}")

print("\nPMNS pairs:")
for i in range(3):
    for j in range(i+1, 3):
        S = concurrence(pmns_states[i], pmns_states[j])
        print(f"  {pmns_words[i]}-{pmns_words[j]}: S={S:.4f}")
