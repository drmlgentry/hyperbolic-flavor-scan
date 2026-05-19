import snappy
import mpmath
mpmath.mp.dps = 50

PHI = (1 + mpmath.sqrt(5)) / 2
LOG_PHI = mpmath.log(PHI)

M = snappy.OrientableClosedCensus[1]
rho = M.polished_holonomy(bits_prec=200)

print("High-precision muon geodesic check")
print(f"PHI    = {mpmath.nstr(PHI, 30)}")
print(f"log(phi) = {mpmath.nstr(LOG_PHI, 30)}")
print()

for word in ['BAAAAAB', 'baaaaab', 'aabbaaa']:
    mat = rho(word)
    # Extract matrix entries as strings then convert
    # snappy polished holonomy entries are sage complex numbers
    # Convert via string representation
    def to_mpc(z):
        s = str(z)
        # format: "a + b*I" or "a - b*I" or just "a"
        s = s.replace('*I','j').replace(' ','')
        try:
            return mpmath.mpc(complex(s))
        except:
            # try splitting manually
            import re
            parts = re.split(r'(?=[+-])', s)
            real = imag = 0.0
            for p in parts:
                if 'j' in p:
                    imag += float(p.replace('j',''))
                elif p:
                    real += float(p)
            return mpmath.mpc(real, imag)
    
    try:
        a = to_mpc(mat[0,0])
        b = to_mpc(mat[0,1])
        c = to_mpc(mat[1,0])
        d = to_mpc(mat[1,1])
    except:
        # Alternative: use numpy via float conversion
        import numpy as np
        mat_np = np.array(rho(word), dtype=complex)
        det = np.linalg.det(mat_np)
        mat_np = mat_np / np.sqrt(det)
        tr_np = np.trace(mat_np)
        disc = tr_np**2 - 4
        lam = (tr_np + np.sqrt(disc)) / 2
        if abs(lam) < 1: lam = 1/lam
        ell = 2 * abs(np.log(lam).real)
        ratio = ell / float(LOG_PHI)
        print(f"{word}: ell={ell:.10f}  ell/log(phi)={ratio:.10f}  "
              f"dist_to_11={abs(ratio-11):.2e}  (numpy, lower precision)")
        continue

    det = a*d - b*c
    sq_det = mpmath.sqrt(det)
    a /= sq_det; b /= sq_det; c /= sq_det; d /= sq_det

    tr = a + d
    disc = tr**2 - 4
    lam1 = (tr + mpmath.sqrt(disc)) / 2
    lam2 = (tr - mpmath.sqrt(disc)) / 2
    lam = lam1 if abs(lam1) > 1 else lam2

    log_lam = mpmath.log(lam)
    ell = 2 * abs(mpmath.re(log_lam))
    ratio = ell / LOG_PHI

    print(f"Word: {word}")
    print(f"  ell         = {mpmath.nstr(ell, 25)}")
    print(f"  11*log(phi) = {mpmath.nstr(11*LOG_PHI, 25)}")
    print(f"  ell/log(phi)= {mpmath.nstr(ratio, 25)}")
    print(f"  dist to 11  = {mpmath.nstr(abs(ratio-11), 5)}")
    print(f"  rel error   = {mpmath.nstr(abs(ratio-11)/11*100, 5)}%")
    print()
