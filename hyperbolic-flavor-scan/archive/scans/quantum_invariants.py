import snappy

for idx, name in [(1, "m003 (PMNS)"), (43, "m006 (CKM)")]:
    M = snappy.OrientableClosedCensus[idx]
    print(f"\n=== {name} (idx {idx}) ===")
    print(f"  Volume:  {float(M.volume()):.6f}")
    print(f"  H1:      {M.homology()}")
    print(f"  Isom:    {M.symmetry_group()}")

    # CS via cusped ancestor
    try:
        name_str = M.name()
        M_cusp = snappy.Manifold(name_str)
        cs = float(M_cusp.chern_simons())
        print(f"  CS (cusped ancestor {name_str}): {cs:.6f}")
    except Exception as e:
        print(f"  CS: {e}")

    # Arithmetic invariants
    for attr in ["invariant_trace_field", "quaternion_algebra",
                 "is_arithmetic", "alexander_polynomial"]:
        try:
            val = getattr(M, attr)()
            print(f"  {attr}: {val}")
        except Exception as e:
            print(f"  {attr}: {e}")

    # All available methods with quantum/torsion/eta content
    quantum_methods = [m for m in dir(M)
                      if any(x in m.lower() for x in
                             ["quantum","jones","witten","turaev",
                              "torsion","eta","reidemeister","reshetikhin",
                              "colored","homfly","kauffman"])]
    print(f"  Quantum methods available: {quantum_methods if quantum_methods else 'none'}")

    # Reidemeister torsion if available
    try:
        rt = M.reidemeister_torsion()
        print(f"  Reidemeister torsion: {rt}")
    except Exception as e:
        pass

    # Eta invariant
    try:
        eta = M.eta_invariant()
        print(f"  Eta invariant: {eta}")
    except:
        pass

    # Length spectrum (first 5 geodesics)
    try:
        ls = M.length_spectrum(3.0)
        print(f"  Length spectrum (len<3):")
        for geo in list(ls)[:5]:
            print(f"    length={float(geo.length.real):.6f} "
                  f"torsion={float(geo.length.imag):.6f} "
                  f"mult={geo.multiplicity}")
    except Exception as e:
        print(f"  Length spectrum: {e}")
