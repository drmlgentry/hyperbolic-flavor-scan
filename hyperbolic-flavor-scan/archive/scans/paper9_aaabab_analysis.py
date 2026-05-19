"""
paper9_aaabab_analysis.py
Compute full geometric data for the near-identity geodesic AAABAB on m006.
"""
import snappy, numpy as np

M = snappy.OrientableClosedCensus[43]
G = M.fundamental_group()

def to_numpy(m):
    return np.array([[complex(m[i,j]) for j in range(2)] for i in range(2)])

word = "AAABAB"
L = to_numpy(G.SL2C(word))

# Eigenvalues
ev = np.linalg.eigvals(L)
lam = ev[np.argmax(np.abs(ev))]
lam_small = ev[np.argmin(np.abs(ev))]

print(f"Word: {word}")
print(f"Eigenvalues: {ev}")
print(f"lambda (large): {lam:.8f}")
print(f"lambda (small): {lam_small:.8f}")
print(f"|lambda|: {abs(lam):.8f}")
print(f"phi = Im(log lambda): {np.degrees(np.angle(lam)):.6f} deg")
print(f"phi_fold: {min(abs(np.degrees(np.angle(lam)))%180, 180-abs(np.degrees(np.angle(lam)))%180):.6f} deg")

# Complex geodesic length: lambda = exp((l + i*phi)/2)
# so l + i*phi = 2*log(lambda)
complex_length = 2 * np.log(lam)
real_length = complex_length.real
twist = complex_length.imag
print(f"\nComplex length l + i*phi = 2*log(lambda):")
print(f"  Real length l = {real_length:.6f}")
print(f"  Twist phi = {np.degrees(twist):.6f} deg")

# Trace
trace = np.trace(L)
print(f"\nTrace Tr(AAABAB) = {trace:.8f}")
print(f"|Tr| = {abs(trace):.8f}")
print(f"For identity: Tr = ±2, so deviation = {abs(abs(trace)-2):.8f}")

# Homology class
h1 = sum(1 if c.islower() else -1 for c in word) % 5
print(f"\nHomology class: sum = {sum(1 if c.islower() else -1 for c in word)}, mod 5 = {h1}")

# Compare with generator traces
print(f"\nGenerator comparison:")
for gen in ["a","b","A","B","aa","bb","AAAAA","BBBBB"]:
    try:
        Lg = to_numpy(G.SL2C(gen))
        tr = np.trace(Lg)
        evg = np.linalg.eigvals(Lg)
        lamg = evg[np.argmax(np.abs(evg))]
        phi_g = np.degrees(np.angle(lamg))
        phi_fold = min(abs(phi_g)%180, 180-abs(phi_g)%180)
        print(f"  {gen:8s}: Tr={tr:.4f}, phi_fold={phi_fold:.3f} deg, |lam|={abs(lamg):.4f}")
    except: pass
