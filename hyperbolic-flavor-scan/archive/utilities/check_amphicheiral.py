import snappy

for idx, name, order in [(39,"idx39",55),(238,"idx238",80),(43,"m006_ckm",5),(1,"m003_pmns",5)]:
    print(f"\n=== {name} (idx {idx}) ===")
    M = snappy.OrientableClosedCensus[idx]
    G = M.symmetry_group()
    
    # Inspect SymmetryGroup attributes
    attrs = [a for a in dir(G) if not a.startswith("__")]
    
    # Try orientation-related properties
    for attr in attrs:
        if any(x in attr.lower() for x in 
               ['amphicheiral','chiral','orient','reflect','mirror','direct','indirect']):
            try:
                val = getattr(G, attr)
                result = val() if callable(val) else val
                print(f"  {attr}: {result}")
            except Exception as e:
                print(f"  {attr}: error {e}")
    
    # Try all attributes to find anything useful
    print(f"  All attrs: {attrs}")
    
    # Try the group order and whether it's consistent with amphicheirality
    # For amphicheiral: Isom contains orientation-reversing elements
    # Index of orientation-preserving subgroup = 2 if amphicheiral
    try:
        print(f"  order: {G.order()}")
        # Try abelian invariants
        try:
            print(f"  abelian_invariants: {G.abelian_invariants()}")
        except: pass
        try:
            print(f"  is_abelian: {G.is_abelian()}")
        except: pass
        try:
            print(f"  is_amphicheiral: {G.is_amphicheiral()}")
        except: pass
        try:
            print(f"  amphicheiral: {G.amphicheiral()}")
        except: pass
    except Exception as e:
        print(f"  error: {e}")
