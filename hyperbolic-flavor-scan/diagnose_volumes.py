import snappy

def vol(M):
    return float(M.volume())

# ── Diagnose the control failure first ──────────────────────────────────────
# idx39 gave vol=1.9627, but m003 gave vol=2.0299 — ratio 0.967, not 2.0
# This means either:
#   (a) OrientableClosedCensus[39] is not the manifold we found before, or
#   (b) m003 is being loaded as a cusped manifold (different vol)

print("=== DIAGNOSING m003 ===")
m003_closed  = snappy.OrientableClosedCensus[1]   # should be m003
m003_cusped  = snappy.Manifold("m003")            # cusped version
print(f"OrientableClosedCensus[1]: vol={vol(m003_closed):.10f}, H1={m003_closed.homology()}, name={m003_closed.name()}")
print(f"snappy.Manifold('m003'):   vol={vol(m003_cusped):.10f}, H1={m003_cusped.homology()}")

print("\n=== DIAGNOSING idx39 ===")
# Search for the manifold with vol ~ 2*0.98136 = 1.96274 and H1=Z/55
target_vol = 2 * vol(m003_closed)
print(f"Looking for vol ~ {target_vol:.8f} (2 * vol(m003_closed))")
census = snappy.OrientableClosedCensus
for i in range(200):
    M = census[i]
    v = vol(M)
    if abs(v - target_vol) < 0.0001:
        print(f"  idx={i}, name={M.name()}, vol={v:.10f}, H1={M.homology()}")

print("\n=== DIAGNOSING m006 ===")
m006_closed = snappy.OrientableClosedCensus[43]
m006_cusped = snappy.Manifold("m006")
print(f"OrientableClosedCensus[43]: vol={vol(m006_closed):.10f}, H1={m006_closed.homology()}, name={m006_closed.name()}")
print(f"snappy.Manifold('m006'):    vol={vol(m006_cusped):.10f}, H1={m006_cusped.homology()}")

print("\n=== CORRECT VOLUMES TO USE ===")
print(f"vol(m003_closed) = {vol(m003_closed):.10f}")
print(f"vol(m006_closed) = {vol(m006_closed):.10f}")
print(f"2*vol(m003)      = {2*vol(m003_closed):.10f}")
print(f"2*vol(m006)      = {2*vol(m006_closed):.10f}")
