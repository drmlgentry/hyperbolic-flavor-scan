import snappy
import numpy as np

# We need to compute the ACTION of each isometry on H1
# SnapPy's Isometry object only exposes cusp data for closed manifolds
# Alternative: use the peripheral structure or try the permutation representation

# Method: for each generator word in pi_1(M), apply the isometry
# and track how homology class changes

for idx, name, p in [(1,"m003",5),(43,"m006",5),(39,"idx39",55),(238,"idx238",80)]:
    print(f"\n=== {name} (idx {idx}, H1=Z/{p}) ===")
    M = snappy.OrientableClosedCensus[idx]
    G = M.symmetry_group()
    print(f"  Isom = {G}, order = {G.order()}")
    
    # Try to get the induced action via the abelianization
    try:
        # SnapPy can compute the abelianization map
        # Try to find a word that generates H1
        fg = M.fundamental_group()
        print(f"  Num generators: {fg.num_generators()}")
        print(f"  Relators: {len(fg.relators())}")
        
        # Try homological action via symmetry group
        try:
            ab = G.abelianization()
            print(f"  G abelianization: {ab}")
        except Exception as e:
            print(f"  abelianization error: {e}")
            
        # Try: for each element of G, compute its action on H1
        # via the induced map on abelianization of pi_1
        isoms = G.isometries()
        print(f"  Num isometries: {len(isoms)}")
        
        # Each isometry permutes the generators of pi_1
        # We can track this via the permutation of relators
        for i, f in enumerate(isoms):
            try:
                # Try all available attributes
                attrs = [a for a in dir(f) if not a.startswith('_')]
                if i == 0:
                    print(f"  Isometry attrs: {attrs}")
                # cusp_images gives permutation of cusps (not useful for closed)
                # Try: does the isometry act on the Dehn filling?
                ci = f.cusp_images()
                cm = f.cusp_maps()
                print(f"  Isom {i}: cusp_images={ci}, cusp_maps={cm}")
            except Exception as e:
                print(f"  Isom {i}: {e}")
    except Exception as e:
        print(f"  Error: {e}")
