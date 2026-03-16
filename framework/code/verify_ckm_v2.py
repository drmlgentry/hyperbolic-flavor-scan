#!/usr/bin/env python3
"""
Arithmetic Flavor Geometry: Verification Computation v2
Fixed sector geometry: Gram matrix construction for CKM.
"""
import numpy as np
from numpy.linalg import norm, eigh
import itertools

np.set_printoptions(precision=6, suppress=True)

PDG = {
    'J': 3.00e-5,
    'theta12': 0.22739, 'theta13': 0.003692, 'theta23': 0.04085,
    'Vud': 0.97435, 'Vus': 0.22500, 'Vub': 0.003690,
    'Vcb': 0.04182, 'Vtb': 0.99918
}

def find_real_roots(coeffs):
    roots = np.roots(coeffs)
    return sorted([r.real for r in roots if abs(r.imag) < 1e-10])

def compute_log_vectors(roots):
    r = roots
    def embed(c):
        a0,a1,a2,a3 = c
        return np.array([a0+a1*ri+a2*ri**2+a3*ri**3 for ri in r])
    def log_vec(c):
        vals = embed(c)
        if any(v <= 0 for v in vals): return None
        lv = np.log(np.abs(vals)); lv -= lv.mean()
        return lv[:3]
    units = []
    for c in itertools.product(range(-4,5), repeat=4):
        if list(c) in [[0,0,0,0],[1,0,0,0],[-1,0,0,0]]: continue
        n = np.prod(embed(c))
        if abs(abs(n)-1.0) < 1e-6:
            lv = log_vec(c)
            if lv is not None and norm(lv) > 0.05:
                units.append((list(c), lv, norm(lv)))
    units.sort(key=lambda x: x[2])
    return units

def select_independent_units(units, n=3):
    sel, vecs = [], []
    for c, lv, ln in units:
        test = np.array(vecs + [lv])
        if len(sel)==0 or np.linalg.matrix_rank(test, tol=1e-8)==len(sel)+1:
            sel.append((c,lv,ln)); vecs.append(lv)
        if len(sel)==n: break
    return sel

def sector_gram(u_dom, u_shared, n_vals, phi=0.0):
    """
    Build 3x3 Gram matrix for generational vectors.
    v_f = n_f * u_dom * exp(i*phi) + u_shared
    All three vectors lie in span(u_dom, u_shared).
    Gram[i,j] = <v_i | v_j>
    """
    ud = u_dom.astype(complex) * np.exp(1j*phi)
    us = u_shared.astype(complex)
    vecs = [n*ud + us for n in sorted(n_vals)]
    G = np.zeros((3,3), dtype=complex)
    for i,vi in enumerate(vecs):
        for j,vj in enumerate(vecs):
            G[i,j] = np.vdot(vi, vj)
    return G

def ckm_from_grams(G_up, G_dn):
    """CKM = U_up^dag U_dn where U diagonalizes each Gram matrix."""
    def diag(G):
        evals, evecs = eigh(G)
        idx = np.argsort(evals)
        return evals[idx], evecs[:,idx]
    _, U_up = diag(G_up)
    _, U_dn = diag(G_dn)
    return U_up.conj().T @ U_dn

def pdg_angles(V):
    Vub=abs(V[0,2]); Vus=abs(V[0,1]); Vcb=abs(V[1,2])
    s13=min(Vub,1.); c13=np.sqrt(max(0,1-s13**2))
    s12=min(Vus/c13,1.) if c13>1e-10 else 0.
    s23=min(Vcb/c13,1.) if c13>1e-10 else 0.
    return np.arcsin(s12), np.arcsin(s13), np.arcsin(s23)

def jarlskog(V):
    return np.imag(V[0,0]*V[1,1]*np.conj(V[0,1])*np.conj(V[1,0]))

def angle_deg(a,b):
    return np.degrees(np.arccos(np.clip(np.dot(a,b)/(norm(a)*norm(b)),-1,1)))

# ---- MAIN ANALYSIS ----

def run_field(name, poly, phi=0.3, up_c=(1,5,9), dn_c=(3,7,12)):
    print(f"\n{'='*65}")
    print(f"FIELD: {name}  phi={phi:.3f}")
    roots = find_real_roots(poly)
    n_real = sum(1 for r in np.roots(poly) if abs(r.imag)<1e-8)
    print(f"Roots: {[f'{r:.5f}' for r in roots]}")
    print(f"Totally real: {n_real==4}")
    if n_real != 4: return None

    units = compute_log_vectors(roots)
    sel = select_independent_units(units, 3)
    if len(sel)<3: print("ERROR: <3 independent units"); return None

    u1,u2,u3 = sel[0][1],sel[1][1],sel[2][1]
    print(f"u1={[f'{x:.4f}' for x in u1]}  |u1|={norm(u1):.4f}")
    print(f"u2={[f'{x:.4f}' for x in u2]}  |u2|={norm(u2):.4f}")
    print(f"u3={[f'{x:.4f}' for x in u3]}  |u3|={norm(u3):.4f}")
    print(f"angles: u1-u2={angle_deg(u1,u2):.1f}  "
          f"u1-u3={angle_deg(u1,u3):.1f}  u2-u3={angle_deg(u2,u3):.1f} deg")

    # Mass check
    lam = [norm(u1), norm(u2), norm(u3)]
    me = 0.511
    print(f"\nLepton mass check (coords 0,2,3 on u1, lam1={lam[0]:.4f}):")
    print(f"  m_mu pred={me*np.exp(2*lam[0]):.1f} MeV  actual=105.7  "
          f"ratio={me*np.exp(2*lam[0])/105.7:.3f}")
    print(f"  m_tau pred={me*np.exp(3*lam[0]):.1f} MeV  actual=1776.9  "
          f"ratio={me*np.exp(3*lam[0])/1776.9:.3f}")
    print(f"Up check (1,5,9 on u2, lam2={lam[1]:.4f}):")
    for n,mact in zip(sorted(up_c),[2.16,1270,172760]):
        mpred=me*np.exp(n*lam[1])
        print(f"  n={n}: pred={mpred:.1f} MeV  actual={mact}  ratio={mpred/mact:.3f}")
    print(f"Down check (3,7,12 on u3, lam3={lam[2]:.4f}):")
    for n,mact in zip(sorted(dn_c),[4.67,93.4,4180]):
        mpred=me*np.exp(n*lam[2])
        print(f"  n={n}: pred={mpred:.1f} MeV  actual={mact}  ratio={mpred/mact:.3f}")

    # CKM
    G_up = sector_gram(u2, u1, up_c, phi)
    G_dn = sector_gram(u3, u1, dn_c, 0.0)
    V = ckm_from_grams(G_up, G_dn)
    J = jarlskog(V)
    t12,t13,t23 = pdg_angles(V)

    print(f"\nCKM |V|:")
    for row in np.abs(V): print(f"  {[f'{x:.5f}' for x in row]}")
    print(f"\ntheta_12: {np.degrees(t12):.3f} deg  (PDG {np.degrees(PDG['theta12']):.3f})  "
          f"ratio={t12/PDG['theta12']:.3f}")
    print(f"theta_13: {np.degrees(t13):.4f} deg  (PDG {np.degrees(PDG['theta13']):.4f})  "
          f"ratio={t13/PDG['theta13']:.3f}")
    print(f"theta_23: {np.degrees(t23):.3f} deg  (PDG {np.degrees(PDG['theta23']):.3f})  "
          f"ratio={t23/PDG['theta23']:.3f}")
    print(f"J = {J:.4e}  (PDG {PDG['J']:.2e})  ratio={J/PDG['J']:.4f}")

    return V, (t12,t13,t23), J, u1, u2, u3

def phase_scan(name, poly, up_c=(1,5,9), dn_c=(3,7,12)):
    print(f"\n--- PHASE SCAN: {name} ---")
    print(f"{'phi':>7} {'J':>12} {'J/PDG':>8} {'t12 deg':>9} {'t23 deg':>9}")
    roots = find_real_roots(poly)
    sel = select_independent_units(compute_log_vectors(roots), 3)
    if len(sel)<3: return
    u1,u2,u3 = sel[0][1],sel[1][1],sel[2][1]
    Js, t12s = [], []
    for phi in np.linspace(0.05,1.5,20):
        try:
            Gu=sector_gram(u2,u1,up_c,phi); Gd=sector_gram(u3,u1,dn_c,0.)
            V=ckm_from_grams(Gu,Gd); J=jarlskog(V)
            t12,t13,t23=pdg_angles(V)
            print(f"  {phi:7.3f} {J:12.4e} {J/PDG['J']:8.3f} "
                  f"{np.degrees(t12):9.3f} {np.degrees(t23):9.3f}")
            Js.append(J); t12s.append(np.degrees(t12))
        except Exception as e:
            print(f"  {phi:7.3f}  {e}")
    if Js:
        print(f"  J range: {min(Js):.3e}--{max(Js):.3e}  "
              f"t12 range: {min(t12s):.2f}--{max(t12s):.2f} deg")

def integer_scan(name, poly, phi=0.3, base_up=(1,5,9), base_dn=(3,7,12)):
    print(f"\n--- INTEGER SCAN: {name}  phi={phi:.2f} ---")
    roots = find_real_roots(poly)
    sel = select_independent_units(compute_log_vectors(roots), 3)
    if len(sel)<3: return
    u1,u2,u3 = sel[0][1],sel[1][1],sel[2][1]
    Gu0=sector_gram(u2,u1,base_up,phi); Gd0=sector_gram(u3,u1,base_dn,0.)
    V0=ckm_from_grams(Gu0,Gd0); J0=jarlskog(V0); t0,_,_=pdg_angles(V0)
    print(f"Base up={base_up} dn={base_dn}: J={J0:.4e} t12={np.degrees(t0):.3f} deg")
    print(f"{'Perturb':<22} {'J':>12} {'J/J0':>8} {'t12':>8}")
    for sec,base,fn in [('up',base_up,lambda c:sector_gram(u2,u1,c,phi)),
                         ('dn',base_dn,lambda c:sector_gram(u3,u1,c,0.))]:
        for idx in range(3):
            for d in [-1,+1]:
                nc=list(base); nc[idx]+=d
                if any(x<=0 for x in nc) or len(set(nc))<3: continue
                try:
                    Gn=fn(tuple(nc))
                    V=ckm_from_grams(Gn if sec=='up' else Gu0,
                                     Gd0 if sec=='up' else Gn)
                    J=jarlskog(V); t,_,_=pdg_angles(V)
                    print(f"  {sec}[{idx}] {base[idx]}->{nc[idx]}: "
                          f"J={J:.4e} J/J0={J/J0:.3f} t12={np.degrees(t):.3f}")
                except Exception as e:
                    print(f"  {sec}[{idx}]{d:+d}: {e}")

def null_tests(name, poly, phi=0.3, up_c=(1,5,9), dn_c=(3,7,12)):
    print(f"\n--- NULL TESTS: {name} ---")
    roots = find_real_roots(poly)
    sel = select_independent_units(compute_log_vectors(roots), 3)
    if len(sel)<3: return
    u1,u2,u3 = sel[0][1],sel[1][1],sel[2][1]
    Gu=sector_gram(u2,u1,up_c,phi); Gd=sector_gram(u3,u1,dn_c,0.)
    V_ref=ckm_from_grams(Gu,Gd); J_ref=jarlskog(V_ref); t_ref,_,_=pdg_angles(V_ref)
    print(f"Reference: J={J_ref:.4e}  t12={np.degrees(t_ref):.3f} deg")

    # Null A: no shared direction
    print("Null A: no shared u1 (sectors use unrelated spanning vectors)")
    rng=np.random.RandomState(42)
    for trial in range(5):
        v1=rng.randn(3); v2=rng.randn(3)
        v1/=norm(v1); v2/=norm(v2)
        try:
            Gu_A=sector_gram(u2,v1,up_c,phi); Gd_A=sector_gram(u3,v2,dn_c,0.)
            V_A=ckm_from_grams(Gu_A,Gd_A); J_A=jarlskog(V_A); t_A,_,_=pdg_angles(V_A)
            print(f"  trial {trial+1}: J={J_A:.4e}  t12={np.degrees(t_A):.3f}")
        except: pass

    # Null B: orthogonalized units
    print("Null B: orthogonalized unit vectors")
    u_orth=[]
    for v in [u1,u2,u3]:
        w=v.astype(float).copy()
        for b in u_orth: w-=np.dot(b,w)*b
        u_orth.append(w/norm(w))
    try:
        Gu_B=sector_gram(u_orth[1],u_orth[0],up_c,phi)
        Gd_B=sector_gram(u_orth[2],u_orth[0],dn_c,0.)
        V_B=ckm_from_grams(Gu_B,Gd_B); J_B=jarlskog(V_B); t_B,_,_=pdg_angles(V_B)
        print(f"  J={J_B:.4e}  t12={np.degrees(t_B):.3f} deg  "
              f"J/J_ref={J_B/J_ref:.3f}")
    except Exception as e: print(f"  Error: {e}")

    # Null C: 20 random unit triples
    print("Null C: 20 random unit triples")
    rng2=np.random.RandomState(7)
    Jlist,t12list=[],[]
    for _ in range(20):
        vs=[rng2.randn(3) for _ in range(3)]
        vs=[v/norm(v) for v in vs]
        try:
            Gu_C=sector_gram(vs[1],vs[0],up_c,phi)
            Gd_C=sector_gram(vs[2],vs[0],dn_c,0.)
            V_C=ckm_from_grams(Gu_C,Gd_C)
            Jlist.append(abs(jarlskog(V_C)))
            t12list.append(np.degrees(pdg_angles(V_C)[0]))
        except: pass
    if Jlist:
        print(f"  J range: {min(Jlist):.3e}--{max(Jlist):.3e}  "
              f"t12 range: {min(t12list):.2f}--{max(t12list):.2f} deg")
        print(f"  (Reference: J={J_ref:.3e}  t12={np.degrees(t_ref):.2f} deg)")

if __name__=="__main__":
    print("="*70)
    print("ARITHMETIC FLAVOR GEOMETRY: VERIFICATION v2")
    print("="*70)

    poly1=[1,-1,-4,4,1]
    poly2=[1,-2,-1,2,1]

    results={}
    for poly,name,key in [(poly1,"Field 1: x^4-x^3-4x^2+4x+1","f1"),
                           (poly2,"Field 2: x^4-2x^3-x^2+2x+1","f2")]:
        r=run_field(name, poly)
        if r: results[key]=(r[2],r[1])  # J, angles
        phase_scan(name, poly)
        integer_scan(name, poly)
        null_tests(name, poly)

    print(f"\n{'='*70}")
    print("CROSS-FIELD SUMMARY")
    print(f"{'='*70}")
    if 'f1' in results and 'f2' in results:
        print(f"\n{'Quantity':<22} {'Field 1':>12} {'Field 2':>12} {'PDG':>12}")
        print("-"*60)
        for key,label in [('f1','Field 1'),('f2','Field 2')]:
            J,a = results[key]
        J1,a1=results['f1']; J2,a2=results['f2']
        print(f"{'J':<22} {J1:12.3e} {J2:12.3e} {PDG['J']:12.3e}")
        print(f"{'theta_12 (deg)':<22} {np.degrees(a1[0]):12.3f} "
              f"{np.degrees(a2[0]):12.3f} {np.degrees(PDG['theta12']):12.3f}")
        print(f"{'theta_13 (deg)':<22} {np.degrees(a1[1]):12.4f} "
              f"{np.degrees(a2[1]):12.4f} {np.degrees(PDG['theta13']):12.4f}")
        print(f"{'theta_23 (deg)':<22} {np.degrees(a1[2]):12.3f} "
              f"{np.degrees(a2[2]):12.3f} {np.degrees(PDG['theta23']):12.3f}")
        print(f"\nJ_F1/J_PDG = {J1/PDG['J']:.4f}")
        print(f"J_F2/J_PDG = {J2/PDG['J']:.4f}")