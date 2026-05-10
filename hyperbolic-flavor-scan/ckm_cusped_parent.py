import snappy

M2 = snappy.OrientableClosedCensus[43]
M3 = snappy.OrientableClosedCensus[1]

print("=== Drilling M2/CKM along its systole ===")
print(f"M2: vol={float(M2.volume()):.8f}, H1={M2.homology()}")

# Drill out the shortest geodesic to get the cusped parent
try:
    ls2 = M2.length_spectrum(2, full_rigor=True)
    systole_info = min(ls2, key=lambda x: float(abs(x.length.real)))
    print(f"Systole: length={float(systole_info.length.real):.8f}, "
          f"twist={float(systole_info.length.imag):.8f}")
    
    # Drill along the systole
    M2_drilled = M2.drill(systole_info)
    print(f"\nDrilled manifold:")
    print(f"  Cusps: {M2_drilled.num_cusps()}")
    print(f"  Volume: {float(M2_drilled.volume()):.8f}")
    print(f"  H1: {M2_drilled.homology()}")
    print(f"  Name (if known): {M2_drilled.identify()}")
except Exception as e:
    print(f"Drilling error: {e}")

print("\n=== Drilling M3/PMNS along its systole ===")
print(f"M3: vol={float(M3.volume()):.8f}, H1={M3.homology()}")

try:
    ls3 = M3.length_spectrum(2, full_rigor=True)
    systole_info3 = min(ls3, key=lambda x: float(abs(x.length.real)))
    print(f"Systole: length={float(systole_info3.length.real):.8f}, "
          f"twist={float(systole_info3.length.imag):.8f}")
    
    M3_drilled = M3.drill(systole_info3)
    print(f"\nDrilled manifold:")
    print(f"  Cusps: {M3_drilled.num_cusps()}")
    print(f"  Volume: {float(M3_drilled.volume()):.8f}")
    print(f"  H1: {M3_drilled.homology()}")
    print(f"  Name (if known): {M3_drilled.identify()}")
except Exception as e:
    print(f"Drilling error: {e}")

# Also: what Dehn filling of what cusped manifold gives M2?
# Try to identify M2 via its canonical retriangulation
print("\n=== Attempting to identify M2 via SnapPy database ===")
try:
    hits = M2.identify()
    print(f"M2 identifications: {hits}")
except Exception as e:
    print(f"Identification error: {e}")

print("\n=== Attempting to identify M3 via SnapPy database ===")
try:
    hits3 = M3.identify()
    print(f"M3 identifications: {hits3}")
except Exception as e:
    print(f"Identification error: {e}")
