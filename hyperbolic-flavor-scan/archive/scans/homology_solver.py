import snappy

M = snappy.OrientableClosedCensus[1]
G = M.fundamental_group()
relators = G.relators()
print("Relators:", relators)

# The relator is a string. We need to count net exponents of 'a' and 'b'.
# Example: "ababAbbAb" -> net a: count('a') - count('A'), net b: count('b') - count('B')
def net_exp(word):
    exp_a = word.count('a') - word.count('A')
    exp_b = word.count('b') - word.count('B')
    return exp_a, exp_b

# Use the first relator (there should be only one for a 1-cusped manifold)
relator = relators[0]
ra, rb = net_exp(relator)
print(f"Relator net exponents: a^{ra}, b^{rb}")

# H1 = Z/5, so ra*a + rb*b ≡ 0 (mod 5)
# We can set a = 1 (generator). Then b ≡ -ra * inv(rb) mod 5 if rb invertible mod 5.
# If rb ≡ 0 mod 5, then the equation gives ra*a ≡ 0 mod 5, so a may be 1 and ra must be multiple of 5.
from sympy import mod_inverse

# Handle cases:
if rb % 5 == 0:
    # Then b is free? Actually H1 = Z/5 means both a and b are in Z/5 and related.
    # If rb ≡ 0 mod 5, then ra*a ≡ 0 mod 5. Since a generates Z/5, ra must be multiple of 5.
    # Then b can be any element? But we need a specific map.
    # Usually the abelianization of a knot complement: a is meridian (1), b is longitude (0 mod 5 for a knot).
    # For m003, it's a knot complement, so one generator (meridian) maps to 1, the other (longitude) maps to 0.
    # Let's check: if rb ≡ 0, then b is in the commutator subgroup? We'll set a=1, b=0.
    print("rb ≡ 0 mod 5: assuming a=1, b=0 (knot complement convention).")
    b_class = 0
else:
    inv = mod_inverse(rb % 5, 5)
    b_class = (-ra * inv) % 5
    print(f"b_class = {b_class}")

# Now compute for our words
words = ['ab', 'bA', 'bba']
for w in words:
    ea, eb = net_exp(w)
    h = (ea * 1 + eb * b_class) % 5
    print(f"{w}: exponents ({ea},{eb}) -> homology class {h}")
