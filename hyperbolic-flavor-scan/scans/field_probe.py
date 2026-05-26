import snappy

M = snappy.Manifold("m003(-2,3)")

print("=== Trace field generators ===")
try:
    print(M.trace_field_gens())
except Exception as e:
    print("ERROR:", e)

print()

print("=== Invariant trace field generators ===")
try:
    print(M.invariant_trace_field_gens())
except Exception as e:
    print("ERROR:", e)

print()

print("=== Tetrahedra field generators ===")
try:
    print(M.tetrahedra_field_gens())
except Exception as e:
    print("ERROR:", e)
