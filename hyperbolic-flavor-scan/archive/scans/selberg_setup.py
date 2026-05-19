"""
selberg_setup.py -- Length spectrum + Selberg zeta zero search
Estimates Laplace eigenvalues of m003 and m006
Fixed: float() conversion on all SnapPy Number objects
"""
import snappy, numpy as np
from scipy.signal import argrelmin
import pandas as pd

for manifold_idx, manifold_name in [(1, 'm003'), (43, 'm006')]:
    M = snappy.OrientableClosedCensus[manifold_idx]
    print(f'\n{"="*60}\n{M.name()} vol={float(M.volume()):.6f}\n{"="*60}')
    vol = float(M.volume())
    try:
        spec = M.length_spectrum(5.0)
        lengths = np.array([float(s.length.real) for s in spec])
        twists  = np.array([float(s.length.imag) for s in spec])
        mult    = np.array([int(s.multiplicity)   for s in spec])
        print(f'Geodesics up to length 5: {len(lengths)}')
        print(f'Shortest geodesic: {lengths[0]:.6f}')

        pd.DataFrame({
            'length':   lengths,
            'twist':    twists,
            'mult':     mult,
        }).to_csv(
            f'C:\\dev\\hyperbolic-flavor-scan\\length_spectrum_{manifold_name}.csv',
            index=False)

        print('\nFirst 15 geodesics (ell, twist_angle_deg, mult):')
        for i in range(min(15, len(lengths))):
            print(f'  ell={lengths[i]:.6f}  phi={np.degrees(twists[i]):.4f} deg  mult={mult[i]}')

        # Selberg zeta |Z(1/2+ir)| zero search
        print('\nSelberg zeta zero search (lambda=1/4+r^2):')
        r_vals   = np.linspace(0.05, np.sqrt(4.75), 600)
        zeta_log = np.zeros(len(r_vals))
        for ell_p in lengths[:200]:
            for k in range(8):
                s_re  = 0.5
                s_im  = r_vals
                log_term = np.log(np.abs(1 - np.exp(-(s_re + k + 1j*s_im)*ell_p)))
                zeta_log += log_term
        zeta_abs = np.exp(zeta_log)
        mins = argrelmin(zeta_abs, order=5)[0]

        print('Candidate eigenvalues lambda_i = 1/4 + r_i^2:')
        for idx in mins[:15]:
            r   = r_vals[idx]
            lam = 0.25 + r**2
            print(f'  r={r:.5f}  lambda={lam:.6f}  |Z|={zeta_abs[idx]:.6f}')

        # Selberg eigenvalue conjecture check: is lambda_1 >= 1/4?
        if len(mins) > 0:
            r1   = r_vals[mins[0]]
            lam1 = 0.25 + r1**2
            print(f'\nSmallest candidate: lambda_1 = {lam1:.6f}')
            if lam1 >= 0.25:
                print('  -> Consistent with Selberg eigenvalue conjecture (lambda_1 >= 1/4)')
            else:
                print('  -> WARNING: lambda_1 < 1/4, exceptional spectrum')

        # Weyl law sanity check
        n_found = len(mins)
        if n_found > 0:
            lam_max = 0.25 + r_vals[mins[-1]]**2
            weyl    = vol * lam_max / (4 * np.pi)
            print(f'\nWeyl law check:')
            print(f'  Found {n_found} candidates up to lambda={lam_max:.3f}')
            print(f'  Weyl prediction N(lambda) ~ {weyl:.1f}  (ratio={n_found/weyl:.2f})')

        # Save zeta values
        pd.DataFrame({
            'r':        r_vals,
            'lambda':   0.25 + r_vals**2,
            'zeta_abs': zeta_abs,
        }).to_csv(
            f'C:\\dev\\hyperbolic-flavor-scan\\selberg_zeta_{manifold_name}.csv',
            index=False)
        print(f'Saved zeta scan to selberg_zeta_{manifold_name}.csv')

    except Exception as e:
        import traceback
        print(f'Error: {e}')
        traceback.print_exc()

print('\nSelberg setup complete.')
