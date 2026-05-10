import snappy

for idx, name, order in [(39,"idx39",55), (238,"idx238",80)]:
    print(f"\n=== {name} ===")
    M = snappy.OrientableClosedCensus[idx]
    print(f"H1={M.homology()}, Isom={M.symmetry_group()}")
    
    try:
        isoms = M.symmetry_group().isometries()
        print(f"Number of isometries: {len(isoms)}")
        for i, f in enumerate(isoms):
            try:
                preserves = f.preserves_orientation()
                mat = f.homology_matrix()
                # mat is a matrix over Z; we want its action on H1=Z/order
                # For cyclic H1, mat is 1x1
                entry = int(mat[0][0]) % order
                print(f"  Isom {i}: orientation_preserving={preserves}, "
                      f"H1 action = x -> {entry}x mod {order}")
                if not preserves:
                    print(f"    *** ORIENTATION-REVERSING ***")
                    if entry == order - 1:
                        print(f"    Acts as full inversion (-1) on Z/{order}")
                        print(f"    -> J=0 on ALL torsion factors")
                    else:
                        print(f"    Acts as x -> {entry}x (partial action)")
            except Exception as e:
                print(f"  Isom {i}: error {e}")
    except Exception as e:
        print(f"isometries() failed: {e}")
        # Fallback: check if amphicheiral directly
        try:
            amp = M.symmetry_group().is_amphicheiral()
            print(f"is_amphicheiral(): {amp}")
        except Exception as e2:
            print(f"is_amphicheiral() also failed: {e2}")
