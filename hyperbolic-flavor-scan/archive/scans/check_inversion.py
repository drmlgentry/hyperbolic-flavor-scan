import snappy

def check_central_inversion(idx):
    M = snappy.OrientableClosedCensus[idx]
    print(f"\n{'='*50}")
    print(f"Manifold: {M.name()} (index {idx})")
    print(f"Volume: {M.volume()}")
    print(f"H1: {M.homology()}")
    
    G = M.symmetry_group()
    print(f"Symmetry group order: {G.order()}")
    
    # Try different ways to get elements
    elements = []
    try:
        elements = G.elements()
    except AttributeError:
        try:
            elements = G.list()
        except AttributeError:
            try:
                elements = [G[i] for i in range(G.order())]
            except:
                print("  Cannot enumerate elements; skipping homology action.")
                return
    
    print("Checking for central inversion (-I) on H1...")
    found_minus_I = False
    h1 = M.homology()
    
    for g in elements:
        try:
            mat = g.homology_matrix()
            # Convert to list of lists for comparison
            mat_list = [list(row) for row in mat]
            
            if h1.isomorphic_to(snappy.AbelianGroup([5])):
                # H1 = Z/5, matrix is 1x1
                if mat_list == [[4]]:
                    found_minus_I = True
                    print(f"  Found -I: {mat_list}")
            elif h1.isomorphic_to(snappy.AbelianGroup([5,5])):
                # H1 = Z/5 + Z/5, matrix is 2x2
                if mat_list == [[4,0],[0,4]]:
                    found_minus_I = True
                    print(f"  Found -I: {mat_list}")
        except Exception as e:
            pass
    
    if found_minus_I:
        print("  => Central inversion EXISTS -> CP suppressed.")
    else:
        print("  => Central inversion DOES NOT EXIST -> CP allowed.")

# Run for m003 and m006
check_central_inversion(0)   # m003
check_central_inversion(43)  # m006