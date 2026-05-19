import snappy
import numpy as np

def homology_class(M, word):
    """
    Compute homology class of word in H1(M) = Z/p.
    Method: abelianize the word by counting signed generator exponents,
    then apply the abelianization map from the fundamental group presentation.
    """
    fg = M.fundamental_group()
    gens = fg.generators()  # e.g. ['a', 'b']
    n = len(gens)
    
    # Count signed exponents in the word
    exponents = {g: 0 for g in gens}
    for c in word:
        if c.islower():
            exponents[c] += 1
        else:
            exponents[c.lower()] -= 1
    
    vec = np.array([exponents[g] for g in gens], dtype=int)
    
    # Get abelianization matrix from relators
    # Each relator gives a row: sum of signed exponents
    relators = fg.relators()
    rel_matrix = []
    for rel in relators:
        row = {g: 0 for g in gens}
        for c in rel:
            if c.islower():
                row[c] += 1
            else:
                row[c.lower()] -= 1
        rel_matrix.append([row[g] for g in gens])
    
    return vec, np.array(rel_matrix, dtype=int)

print("=== m006 (CKM, idx 43) ===")
M6 = snappy.OrientableClosedCensus[43]
fg6 = M6.fundamental_group()
print(f"Generators: {fg6.generators()}")
print(f"Relators: {fg6.relators()}")
print(f"H1: {M6.homology()}")
print()

for word in ["aaB", "AbA", "AAb", "bAA", "aBa"]:
    vec, rel_mat = homology_class(M6, word)
    print(f"[{word}]: exponent vector = {vec}")
    # The homology class mod 5
    # For Z/5 with 2 generators a,b and relators,
    # we need the Smith normal form to find the class
    # Simple approach: compute the class mod 5 using the
    # known fact that H1=Z/5 has a single generator
    # Try: compute sum mod 5 (rough proxy)
    print(f"         sum mod 5 = {sum(vec) % 5}")
    print(f"         relators act as: {rel_mat}")

print()
print("=== m003 (PMNS, idx 1) ===")
M3 = snappy.OrientableClosedCensus[1]
fg3 = M3.fundamental_group()
print(f"Generators: {fg3.generators()}")
print(f"Relators: {fg3.relators()}")
print(f"H1: {M3.homology()}")
print()

for word in ["aa","ab","aB","ba","bA","aab","aaB","abA"]:
    vec, rel_mat = homology_class(M3, word)
    print(f"[{word}]: exponent vector = {vec}, sum mod 5 = {sum(vec) % 5}")
