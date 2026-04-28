import snappy
import numpy as np

phi = (1 + np.sqrt(5)) / 2
log_phi = np.log(phi)

def get_lengths(M, cutoff=5.0):
    return np.array(sorted([float(x.length.real())
                             for x in M.length_spectrum(cutoff)]))

def build_kernel(L, s):
    d = L[:,None] - L[None,:]
    return np.exp(-d*d/(2*s*s))

def eff_rank(K):
    ev = np.linalg.eigvalsh(K)
    ev = np.maximum(ev, 1e-12)
    p = ev/ev.sum()
    return np.exp(-np.sum(p*np.log(p)))

def energy_frac(K, n):
    ev = np.sort(np.linalg.eigvalsh(K))[::-1]
    return ev[:n].sum()/ev.sum()

M1 = snappy.Manifold('m003'); M1.dehn_fill((-2,3))
M2 = snappy.Manifold('m006'); M2.dehn_fill((-5,2))
L1 = get_lengths(M1, 5.0)
L2 = get_lengths(M2, 5.0)
M3 = snappy.Manifold('m003')
M4 = snappy.Manifold('m006')
L3 = get_lengths(M3, 5.0)
L4 = get_lengths(M4, 5.0)

print('PMNS (m003(-2,3)): %d geodesics' % len(L1), flush=True)
print('CKM  (m006(-5,2)): %d geodesics' % len(L2), flush=True)
print('m003 cusped:       %d geodesics' % len(L3), flush=True)
print('m006 cusped:       %d geodesics' % len(L4), flush=True)
print(flush=True)

print('%8s  %9s  %9s  %9s  %9s  note' % (
    'sigma','R_PMNS','R_CKM','R_m003c','R_m006c'), flush=True)
print('-'*70, flush=True)

for sigma in [round(0.05*i,2) for i in range(4,29)]:
    K1 = build_kernel(L1, sigma)
    K2 = build_kernel(L2, sigma)
    K3 = build_kernel(L3, sigma)
    K4 = build_kernel(L4, sigma)
    R1 = eff_rank(K1)
    R2 = eff_rank(K2)
    R3 = eff_rank(K3)
    R4 = eff_rank(K4)
    note = ''
    if abs(sigma - log_phi) < 0.03:       note = '<-- log(phi)'
    elif abs(sigma - 1.5*log_phi) < 0.03: note = '<-- (3/2)log(phi)'
    if 0.30 <= sigma <= 1.10 or note:
        print('%8.4f  %9.4f  %9.4f  %9.4f  %9.4f  %s' % (
            sigma, R1, R2, R3, R4, note), flush=True)

print(flush=True)
print('ENERGY FRACTIONS:', flush=True)
print('%8s  %10s  %10s  %10s  %10s' % (
    'sigma','E3(PMNS)','E2(PMNS)','E3(CKM)','E2(CKM)'), flush=True)
print('-'*55, flush=True)
for sigma in [log_phi, 1.5*log_phi]:
    K1 = build_kernel(L1, sigma)
    K2 = build_kernel(L2, sigma)
    print('%8.4f  %10.4f  %10.4f  %10.4f  %10.4f' % (
        sigma,
        energy_frac(K1,3), energy_frac(K1,2),
        energy_frac(K2,3), energy_frac(K2,2)), flush=True)

print(flush=True)
print('EIGENVALUE GAPS AT log(phi):', flush=True)
for label, L in [('PMNS', L1), ('CKM', L2)]:
    K = build_kernel(L, log_phi)
    ev = np.sort(np.linalg.eigvalsh(K))[::-1]
    ev_n = ev/ev[0]
    gaps = ev_n[:-1]/ev_n[1:]
    print('  %s: lambda_0..5 = %s' % (
        label, str([round(float(x),4) for x in ev_n[:6]])), flush=True)
    print('       gaps 0/1 1/2 2/3 3/4 = %s' % (
        str([round(float(x),4) for x in gaps[:4]])), flush=True)

print(flush=True)
print('SPACING ANALYSIS:', flush=True)
for label, L in [('PMNS',L1),('CKM',L2),('m003c',L3),('m006c',L4)]:
    sp = np.diff(L)
    rt = np.exp(sp)
    n_lp = int(np.sum(np.abs(sp - log_phi) < 0.025))
    n_ph = int(np.sum(np.abs(rt - phi) < 0.05))
    print('  %6s: mean_sp=%.5f  mean_ratio=%.5f  '
          'near_logphi=%d/%d  near_phi=%d/%d' % (
        label, sp.mean(), rt.mean(),
        n_lp, len(sp), n_ph, len(rt)), flush=True)
