import snappy
import numpy as np
from scipy.linalg import logm as scipy_logm
import itertools

phi = (1 + 5**0.5) / 2
log_phi = np.log(phi)
THETA_W_DEG = 28.1826
TOL = 1.0

def axis_imlogm(A):
    L = scipy_logm(A)
    imL = L.imag
    nx = float((imL[0,1] + imL[1,0]).real)
    ny = float((1j*(imL[1,0] - imL[0,1])).real)
    nz = float((imL[0,0] - imL[1,1]).real)
    n = np.array([nx, ny, nz])
    norm = np.linalg.norm(n)
    if norm < 1e-10:
        return None
    return n / norm

def pairwise_angle_deg(n1, n2):
    cos_a = np.clip(np.dot(n1, n2), -1.0, 1.0)
    return np.degrees(np.arccos(abs(cos_a)))

def gen_words(gens, max_len):
    words = []
    for length in range(1, max_len+1):
        for w in itertools.product(gens, repeat=length):
            word = "".join(w)
            skip = False
            for i in range(len(word)-1):
                if (word[i]=="a" and word[i+1]=="A") or \
                   (word[i]=="A" and word[i+1]=="a") or \
                   (word[i]=="b" and word[i+1]=="B") or \
                   (word[i]=="B" and word[i+1]=="b"):
                    skip = True
                    break
            if not skip:
                words.append(word)
    return words

gens = ["a","b","A","B"]
words = gen_words(gens, 4)
print(f"Words generated: {len(words)}", flush=True)

for label, idx in [("m006 (CKM)", 43), ("m003 (PMNS)", 1)]:
    print(f"\n=== {label} ===", flush=True)
    M = snappy.OrientableClosedCensus[idx]
    G = M.fundamental_group()
    axes = {}
    for word in words:
        try:
            mat = G.SL2C(word)
            A = np.array([[complex(mat[i][j]) for j in range(2)] for i in range(2)])
            A = A / np.sqrt(np.linalg.det(A))
            n = axis_imlogm(A)
            if n is not None:
                axes[word] = n
        except Exception:
            pass
    print(f"Axes extracted: {len(axes)}", flush=True)

    hits = []
    pairs = list(itertools.combinations(axes.keys(), 2))
    print(f"Scanning {len(pairs)} pairs...", flush=True)
    for w1, w2 in pairs:
        angle = pairwise_angle_deg(axes[w1], axes[w2])
        if abs(angle - THETA_W_DEG) < TOL:
            hits.append((angle, w1, w2))
    hits.sort(key=lambda x: abs(x[0] - THETA_W_DEG))

    print(f"Hits within {TOL} deg of theta_W={THETA_W_DEG}:")
    if hits:
        print(f"{'Angle':>10}  {'Diff':>8}  {'Word1':<10}  Word2")
        for angle, w1, w2 in hits[:20]:
            diff = angle - THETA_W_DEG
            print(f"{angle:>10.4f}  {diff:>8.4f}  {w1:<10}  {w2}")
    else:
        print("  No hits.")
        all_p = []
        for w1, w2 in list(itertools.combinations(axes.keys(), 2))[:500]:
            angle = pairwise_angle_deg(axes[w1], axes[w2])
            all_p.append((abs(angle - THETA_W_DEG), angle, w1, w2))
        all_p.sort()
        for diff, angle, w1, w2 in all_p[:5]:
            print(f"  {angle:.4f} deg (diff {diff:.4f})  {w1} vs {w2}")
