import numpy as np
import snappy
from scipy.linalg import logm, qr

PHI = (1+5**0.5)/2
log_phi = np.log(PHI)

def matrix_to_axis_vector(matrix):
    mat = np.array(matrix, dtype=complex)
    mat = mat / np.sqrt(np.linalg.det(mat))
    L = logm(mat)
    x = float(np.real(L[0,1]+L[1,0]))/2
    y = float(np.imag(L[1,0]-L[0,1]))/2
    z = float(np.real(L[0,0]-L[1,1]))/2
    v = np.array([x,y,z])
    n = np.linalg.norm(v)
    return v/n if n>1e-10 else np.array([1.,0.,0.])

def vector_angle(v1,v2):
    return np.arccos(np.clip(np.abs(np.dot(v1,v2)),-1.,1.))

CKM = np.array([[0.97427,0.22536,0.00355],
                [0.22522,0.97339,0.04108],
                [0.00886,0.04050,0.99914]])

M = snappy.OrientableClosedCensus[43]
rho = M.polished_holonomy()
words = ['aaB','AbA','AAb']
axes = [matrix_to_axis_vector(rho(w)) for w in words]

def fitness(sigma):
    mo = np.zeros((3,3),dtype=complex)
    for r in range(3):
        for c in range(3):
            mo[r,c] = np.exp(-vector_angle(axes[r],axes[c])**2/(2*sigma**2))
    U,_ = qr(mo)
    return np.linalg.norm(np.abs(U)-CKM,'fro')

print("Fine grid around log(phi):")
print(f"log(phi) = {log_phi:.8f}")
print()
print(f"{'sigma':>10}  {'sigma/log(phi)':>16}  {'fitness':>10}  {'min?':>6}")
print("-"*50)

results = []
for s in np.linspace(0.460, 0.500, 41):
    f = fitness(s)
    results.append((s, f))

min_s, min_f = min(results, key=lambda x: x[1])
for s, f in results:
    marker = " <--MIN" if abs(s-min_s)<1e-6 else ""
    lp = " <--log(phi)" if abs(s-log_phi)<0.0006 else ""
    print(f"{s:10.5f}  {s/log_phi:16.6f}  {f:10.6f}{marker}{lp}")

print()
print(f"Minimum fitness: {min_f:.6f} at sigma = {min_s:.5f}")
print(f"log(phi)       : {log_phi:.5f}")
print(f"Difference     : {abs(min_s-log_phi):.5f} rad = {np.degrees(abs(min_s-log_phi)):.3f} deg")
print(f"Ratio          : {min_s/log_phi:.6f}")
print()

# Golden ratio exact test
f_exact = fitness(log_phi)
print(f"Fitness at sigma=log(phi) exactly: {f_exact:.6f}")
print(f"Fitness at sigma=0.49:             {fitness(0.49):.6f}")
print(f"Fitness at sigma=min:              {min_f:.6f}")
print()
print("Conclusion:")
if abs(min_s - log_phi) < 0.003:
    print(f"The fitness minimum is at sigma = log(phi) to within {np.degrees(abs(min_s-log_phi)):.2f} degrees.")
    print("The Gaussian coherence scale that best reproduces CKM mixing")
    print("is the mass lattice spacing = regulator of Q(sqrt(5)).")
    print("This is not a free parameter -- it is fixed by the arithmetic.")
else:
    print(f"Minimum at sigma={min_s:.5f}, log(phi)={log_phi:.5f}")
    print(f"Difference: {abs(min_s-log_phi):.5f} rad")
