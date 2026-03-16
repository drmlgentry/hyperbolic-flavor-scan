import snappy, numpy as np

M = snappy.OrientableClosedCensus[1]
G = M.fundamental_group()

# Correct conversion from SimpleMatrix to numpy array
mat = G.SL2C('a')
print(f"Type: {type(mat)}")
print(f"mat: {mat}")

# Convert via list of lists
L = np.array([[complex(mat[i,j]) for j in range(2)] for i in range(2)])
print(f"numpy array:\n{L}")

ev = np.linalg.eigvals(L)
lam = ev[np.argmax(np.abs(ev))]
phi = float(np.degrees(np.angle(lam)))
print(f"phi(a) = {phi:.4f} deg -- SUCCESS")
