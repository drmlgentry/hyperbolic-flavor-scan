import snappy

for idx, name in [(1, "m003"), (43, "m006")]:
    M = snappy.OrientableClosedCensus[idx]
    print(f"\n=== {name} (idx {idx}) ===")
    print(f"  Vol={float(M.volume()):.6f}, H1={M.homology()}")

    # ── Hyperbolic torsion (Reidemeister-Ray-Singer) ───────────────
    # This is the analytic torsion of the hyperbolic representation
    # Directly related to the one-loop quantum gravity partition function
    for method, desc in [
        ("hyperbolic_torsion",
         "Reidemeister torsion (std rep)"),
        ("hyperbolic_adjoint_torsion",
         "Adjoint torsion (sl2C rep)"),
    ]:
        try:
            tor = getattr(M, method)()
            print(f"\n  {desc}:")
            print(f"    {tor}")
            print(f"    type: {type(tor)}")
            try:
                print(f"    float: {complex(tor):.6f}")
            except:
                try:
                    print(f"    str: {str(tor)[:80]}")
                except: pass
        except Exception as e:
            print(f"  {method}: {e}")

    # ── SL(N) torsion for N=2,3 ───────────────────────────────────
    for N in [2, 3, 4]:
        try:
            tor_N = M.hyperbolic_SLN_torsion(N)
            print(f"\n  SL({N}) torsion:")
            print(f"    {str(tor_N)[:100]}")
        except Exception as e:
            print(f"  SL({N}) torsion: {e}")

    # ── Length spectrum analysis ───────────────────────────────────
    # Torsion angles from length spectrum = twist angles from CKM/PMNS paper
    # These are Im(complex_length) = the geodesic holonomy angles
    print(f"\n  Length spectrum (first 10 geodesics):")
    print(f"  {'Re(len)':>10} {'Im(len)=torsion':>16} {'mult':>5} "
          f"{'|torsion/pi|':>14}")
    ls = M.length_spectrum(4.0)
    for geo in list(ls)[:10]:
        re = float(geo.length.real)
        im = float(geo.length.imag)
        mult = geo.multiplicity
        print(f"  {re:>10.6f} {im:>16.6f} {mult:>5d} "
              f"{abs(im/3.14159265):>14.6f}")

    # ── CS invariant context ──────────────────────────────────────
    name_str = M.name()
    M_cusp = snappy.Manifold(name_str)
    cs = float(M_cusp.chern_simons())
    print(f"\n  CS(cusped) = {cs:.6f}")
    print(f"  CS * 4     = {cs*4:.6f}  (should be integer if arithmetic)")
    print(f"  CS * 5     = {cs*5:.6f}")
    print(f"  CS * 20    = {cs*20:.6f}")

    # ── Connection to quantum gravity ────────────────────────────
    # The one-loop quantum gravity partition function on M is:
    # Z_1loop = |Tor(M, Ad rho)|^(1/2) * exp(i * pi * CS(M) * k)
    # where k = Chern-Simons level
    # For k=5 (matching H1=Z/5):
    import cmath
    for k in [1, 4, 5, 10, 20]:
        phase = cmath.exp(2j * 3.14159265 * cs * k)
        print(f"  exp(2pi*i*CS*k={k:2d}) = {phase.real:+.6f} {phase.imag:+.6f}i")
