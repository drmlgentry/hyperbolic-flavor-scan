"""
arithmetic_invariants.py
Compute arithmetic invariants of m003 and m006 using SnapPy.
These determine whether the manifolds are arithmetic and what their
trace fields and quaternion algebras are.
"""
import snappy

for idx, name in [(1,"m003"),(43,"m006")]:
    M = snappy.OrientableClosedCensus[idx]
    print(f"\n{'='*55}")
    print(f"{name}  vol={float(M.volume()):.6f}  H1={M.homology()}")
    
    try:
        # Invariant trace field
        print(f"Invariant trace field: {M.invariant_trace_field()}")
    except Exception as e:
        print(f"Invariant trace field: {e}")
    
    try:
        # Trace field
        print(f"Trace field: {M.trace_field()}")
    except Exception as e:
        print(f"Trace field: {e}")
        
    try:
        # Is it arithmetic?
        print(f"Is arithmetic: {M.is_arithmetic()}")
    except Exception as e:
        print(f"Is arithmetic: {e}")

    try:
        # Chern-Simons invariant
        print(f"Chern-Simons: {M.chern_simons()}")
    except Exception as e:
        print(f"Chern-Simons: {e}")

    try:
        # Alexander polynomial (requires Sage)
        print(f"Eta invariant: {M.eta_invariant()}")
    except Exception as e:
        print(f"Eta invariant: {e}")
