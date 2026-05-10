import snappy
import numpy as np

phi = (1+np.sqrt(5))/2
log_phi = np.log(phi)

def matrix_from_spec_item(item):
    """Extract 2x2 complex matrix from length spectrum item."""
    raw = item['matrix']
    # raw is a 4x4 real matrix encoding the complex 2x2
    # SnapPy stores loxodromic as O31 matrix or as SL2C
    # Try to get SL2C directly
    return raw

def axis_from_SL2C(A):
    """
    A is a 2x2 complex matrix in SL(2,C).
    Fixed points of Mobius z->(az+b)/(cz+d):
    c*z^2 + (d-a)*z - b = 0
    Returns unit vector on S^2 via stereographic projection.
    """
    a, b = complex(A[0,0]), complex(A[0,1])
    c, d = complex(A[1,0]), complex(A[1,1])

    if abs(c) < 1e-10:
        if abs(d-a) < 1e-10:
            return None
        z = -b/(d-a)
    else:
        disc = (d-a)**2 + 4*b*c
        sq = np.sqrt(disc + 0j)
        z1 = (-(d-a) + sq)/(2*c)
        z2 = (-(d-a) - sq)/(2*c)
        # Use fixed point with |z|>1 (repelling) or just z1
        z = z1 if abs(z1) >= abs(z2) else z2

    # Stereographic: z=x+iy -> S^2
    if abs(z) > 1e8:
        return np.array([0.,0.,1.])
    x, y = z.real, z.imag
    r2 = x*x + y*y
    d2 = 1.0 + r2
    n = np.array([2*x/d2, 2*y/d2, (r2-1)/d2])
    nrm = np.linalg.norm(n)
    return n/nrm if nrm > 1e-10 else None

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

def energy_frac(K, n):
    ev = np.sort(np.linalg.eigvalsh(K))[::-1]
    return ev[:n].sum()/ev.sum()

print("="*60, flush=True)
print("HOLONOMY AXIS KERNEL ON S^2", flush=True)
print("Using G.SL2C() and length_spectrum matrices", flush=True)
print("="*60, flush=True)

all_axes = {}

for label, mname, filling in [
    ('PMNS','m003',(-2,3)),
    ('CKM', 'm006',(-5,2)),
]:
    print(f"\n--- {label}: {mname}{filling} ---", flush=True)
    M = snappy.Manifold(mname)
    M.dehn_fill(filling)

    # Method: use G.SL2C() to get holonomy matrices for generators
    G = M.fundamental_group(simplify_presentation=True)
    gens = G.generators()
    print(f"  Generators: {gens}", flush=True)

    # G.SL2C(word) returns SL(2,C) matrix for a word
    axes = []
    words_tried = []

    # Build words up to length 3
    import itertools
    all_words = list(gens)
    # Add length-2 words
    for g1,g2 in itertools.product(gens+[g.upper() for g in gens], repeat=2):
        all_words.append(g1+g2)
    # Add length-3 words
    for g1,g2,g3 in itertools.product(gens, repeat=3):
        all_words.append(g1+g2+g3)

    seen_axes = []
    for word in all_words[:200]:
        try:
            A = np.array(G.SL2C(word), dtype=complex)
            tr = abs(np.trace(A))
            if tr <= 2.001:
                continue  # not loxodromic
            n = axis_from_SL2C(A)
            if n is None or np.any(np.isnan(n)):
                continue
            # Check not duplicate
            is_dup = any(abs(np.dot(n, m)) > 0.9999 for m in seen_axes)
            if not is_dup:
                seen_axes.append(n)
                axes.append((word, tr, n))
        except:
            pass

    print(f"  Found {len(axes)} unique axis directions", flush=True)

    # Also try from length spectrum matrices
    print(f"  Adding axes from length_spectrum...", flush=True)
    try:
        spec = M.length_spectrum(4.0)
        for item in spec:
            try:
                raw = item['matrix']
                # raw is 4x4 real O31 matrix - convert via polished_holonomy
                pass
            except: pass
    except: pass

    # Use polished_holonomy for better precision
    print(f"  Trying polished_holonomy...", flush=True)
    try:
        PH = M.polished_holonomy()
        print(f"  polished_holonomy type: {type(PH)}", flush=True)
        print(f"  dir: {[x for x in dir(PH) if not x.startswith('_')][:15]}", flush=True)
        for g in gens:
            try:
                A = np.array(PH(g), dtype=complex)
                tr = abs(np.trace(A))
                if tr > 2.001:
                    n = axis_from_SL2C(A)
                    if n is not None:
                        is_dup = any(abs(np.dot(n,m)) > 0.9999 for m in seen_axes)
                        if not is_dup:
                            seen_axes.append(n)
                            axes.append((g+'_polished', tr, n))
                            print(f"    Added polished axis for '{g}': tr={tr:.4f}", flush=True)
            except Exception as e:
                print(f"    polished {g}: {e}", flush=True)
    except Exception as e:
        print(f"  polished_holonomy error: {e}", flush=True)

    if len(axes) < 3:
        print(f"  WARNING: only {len(axes)} axes found", flush=True)
        all_axes[label] = None
        continue

    ax_array = np.array([a[2] for a in axes])
    all_axes[label] = ax_array

    print(f"\n  First 5 axes and pairwise angles:", flush=True)
    for i,(w,tr,n) in enumerate(axes[:5]):
        print(f"    [{i}] '{w}' tr={tr:.4f}  n={np.round(n,4)}", flush=True)

    print(f"\n  Kernel rank sweep:", flush=True)
    print(f"  {'sigma':>8} {'R_eff':>8} {'E3':>8} {'E2':>8}  note", flush=True)
    print("  "+"-"*48, flush=True)

    for sigma in np.arange(0.20, 1.21, 0.10):
        K = build_kernel(ax_array, sigma)
        R = eff_rank(K)
        E3 = energy_frac(K, 3)
        E2 = energy_frac(K, 2)
        note = ''
        if abs(sigma-log_phi) < 0.06:       note = '<-- log(phi)'
        elif abs(sigma-1.5*log_phi) < 0.06: note = '<-- (3/2)log(phi)'
        print(f"  {sigma:>8.4f} {R:>8.4f} {E3:>8.4f} {E2:>8.4f}  {note}", flush=True)

    # Eigenvalues at log(phi)
    K_lp = build_kernel(ax_array, log_phi)
    ev = np.sort(np.linalg.eigvalsh(K_lp))[::-1]
    ev_n = ev/ev[0]
    print(f"\n  Eigenvalues at sigma=log(phi):", flush=True)
    for i in range(min(6,len(ev_n)-1)):
        print(f"    lambda_{i}={ev_n[i]:.5f}  gap={ev_n[i]/ev_n[i+1]:.4f}", flush=True)

# Final comparison
print("\n"+"="*60, flush=True)
print("COMPARISON", flush=True)
for sigma,name in [(log_phi,'log(phi)'),(1.5*log_phi,'(3/2)log(phi)')]:
    print(f"\nsigma={sigma:.5f} ({name}):", flush=True)
    for label in ['PMNS','CKM']:
        ax = all_axes.get(label)
        if ax is not None:
            K = build_kernel(ax, sigma)
            print(f"  {label}: R_eff={eff_rank(K):.4f}  "
                  f"E3={energy_frac(K,3):.4f}  N={len(ax)}", flush=True)
