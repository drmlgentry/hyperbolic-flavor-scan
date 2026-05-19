import snappy, numpy as np

M = snappy.OrientableClosedCensus[43]
G = M.fundamental_group()

def to_numpy(m):
    return np.array([[complex(m[i,j]) for j in range(2)] for i in range(2)])

# Analyze the pure-translation geodesic
word = "aaBBABABB"
L = to_numpy(G.SL2C(word))
ev = np.linalg.eigvals(L)
lam = ev[np.argmax(np.abs(ev))]

print(f"Word: {word}")
print(f"Eigenvalue lambda: {lam:.8f}")
print(f"|lambda|: {abs(lam):.8f}")
print(f"phi = Im(log lambda): {np.degrees(np.angle(lam)):.8f} deg")
print(f"phi_fold: {min(abs(np.degrees(np.angle(lam)))%180, 180-abs(np.degrees(np.angle(lam)))%180):.8f} deg")
print(f"Real length l: {2*np.log(abs(lam)):.8f}")
print(f"Trace: {np.trace(L):.8f}")
h1 = sum(1 if c.islower() else -1 for c in word) % 5
print(f"Homology class: {h1}")

# Is lambda real positive? That would mean pure translation
print(f"\nlambda real? {abs(lam.imag) < 1e-6}")
print(f"lambda positive? {lam.real > 0}")
print(f"phi exactly 0? {abs(np.degrees(np.angle(lam))) < 1e-4}")

# Summary of first near-zero per class
print("\n\nSummary: first length achieving phi_fold < 0.01 per class:")
import pandas as pd
df = pd.read_csv(r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len10_m006.csv")
df["h1_class"] = df["word"].apply(lambda w: sum(1 if c.islower() else -1 for c in w) % 5)
for k in range(5):
    sub = df[(df.h1_class == k) & (df.phi_fold < 0.01)].sort_values("length")
    if len(sub):
        row = sub.iloc[0]
        print(f"  Class {k}: first at length {int(row.length)}, word={row.word}, phi_fold={row.phi_fold:.6f}")
    else:
        print(f"  Class {k}: no phi_fold < 0.01 found at length <= 10")
