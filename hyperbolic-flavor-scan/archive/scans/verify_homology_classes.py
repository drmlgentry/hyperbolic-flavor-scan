import snappy

print("=== HOMOLOGY CLASSES OF GEODESIC TRIPLE WORDS ===")
print()

# CKM manifold
print("m006 (CKM, idx 43, H1=Z/5):")
M6 = snappy.OrientableClosedCensus[43]
for word in ["aaB", "AbA", "AAb", "bAA", "aBa"]:
    try:
        h = M6.homology_of_word(word)
        print(f"  [{word}] = {h}")
    except Exception as e:
        print(f"  [{word}]: {e}")

print()
print("m003 (PMNS, idx 1, H1=Z/5):")
M3 = snappy.OrientableClosedCensus[1]
for word in ["aa","ab","aB","ba","bA","AB","aab","aaB","abA","aBa"]:
    try:
        h = M3.homology_of_word(word)
        print(f"  [{word}] = {h}")
    except Exception as e:
        print(f"  [{word}]: {e}")
