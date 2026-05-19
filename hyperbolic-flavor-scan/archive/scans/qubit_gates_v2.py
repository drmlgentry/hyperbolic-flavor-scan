import snappy
import numpy as np
from scipy.linalg import logm, qr

def get_axis_and_angle(rho, word):
    mat = np.array(rho(word), dtype=complex)
    mat = mat / np.sqrt(np.linalg.det(mat))
    L = logm(mat)
    x = float(np.real(L[0,1]+L[1,0]))/2
    y = float(np.imag(L[1,0]-L[0,1]))/2
    z = float(np.real(L[0,0]-L[1,1]))/2
    v = np.array([x,y,z])
    norm = np.linalg.norm(v)
    # Reduce rotation angle to [0, 2pi)
    angle_full = norm * 2
    angle_reduced = angle_full % (2*np.pi)
    axis = v / norm if norm > 1e-10 else v
    return axis, angle_reduced, angle_full

PMNS_PDG = np.array([[0.821,0.550,0.148],
                     [0.357,0.339,0.871],
                     [0.442,0.762,0.471]])

PERMS = [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]]

def borel_matrix(axes, l21, l31, l32):
    Lm = np.array([[1.,0.,0.],[l21,1.,0.],[l31,l32,1.]])
    Q,_ = qr(Lm)
    return np.abs(Q)

print("=== Qubit gate analysis of m003 PMNS word triple ===")
print("(Rotation angles reduced mod 2pi)")
print()

M = snappy.OrientableClosedCensus[1]
rho = M.polished_holonomy()

words_base = ['aBaB', 'baab', 'baaB']
axes_base = []

for w in words_base:
    ax, ang_red, ang_full = get_axis_and_angle(rho, w)
    axes_base.append(ax)
    print(f"Word {w}:")
    print(f"  Bloch axis:       ({ax[0]:.5f}, {ax[1]:.5f}, {ax[2]:.5f})")
    print(f"  Rotation (full):  {np.degrees(ang_full):.2f} deg")
    print(f"  Rotation (mod 2π): {np.degrees(ang_red):.4f} deg = {ang_red:.6f} rad")

print()
print("=== Pairwise axis angles (no abs — signed dot product) ===")
pairs = [(0,1,'θ₁₂'),(0,2,'θ₁₃'),(1,2,'θ₂₃')]
for i,j,name in pairs:
    dot = np.dot(axes_base[i], axes_base[j])
    angle = np.degrees(np.arccos(np.clip(dot,-1,1)))
    print(f"  {name} = {angle:.4f} deg  (dot={dot:.6f})")

print()
print("=== What IS invariant: the mixing matrix ===")
for idx, label, triple in [
    (1,   "m003 base",   ['aBaB','baab','baaB']),
    (39,  "deg-2 cover", ['bbb','aBBB','bAAb']),
    (238, "deg-3 cover", ['AAb','aaaa','aBAb']),
]:
    Mi = snappy.OrientableClosedCensus[idx]
    ri = Mi.polished_holonomy()
    axs = []
    for w in triple:
        ax,_,_ = get_axis_and_angle(ri, w)
        axs.append(ax)
    # Borel construction with optimal parameters from earlier scan
    l21,l31,l32 = -6.4381, 3.4996, -0.8368
    Lm = np.array([[1.,0.,0.],[l21,1.,0.],[l31,l32,1.]])
    Q,_ = qr(Lm)
    Qabs = np.abs(Q)
    best = min(PERMS, key=lambda p: np.linalg.norm(Qabs[:,p]-PMNS_PDG,'fro'))
    fitness = np.linalg.norm(Qabs[:,best]-PMNS_PDG,'fro')
    print(f"\n{label} (idx={idx}):")
    print(f"  Fitness: {fitness:.5f}")
    row = Qabs[:,best]
    for r in row:
        print(f"    {r[0]:.4f}  {r[1]:.4f}  {r[2]:.4f}")
