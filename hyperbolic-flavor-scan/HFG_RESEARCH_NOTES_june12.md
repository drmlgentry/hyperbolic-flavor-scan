# HFG Research Notes — June 12, 2026

## Census Scan Results: All Six H1=Z/5 Manifolds (first 300)

| idx | Name    | Volume  | vol/v0 | theta* | Symmetry       | ITF           |
|-----|---------|---------|--------|--------|----------------|---------------|
|   1 | m003    | 0.98137 | 1.000  | -180°  | none           | Q(sqrt(-3))   |
|  11 | m004    | 1.52948 | 1.559  | +120°  | none           | Q(sqrt(-3))   |
|  43 | m006    | 2.02885 | 2.067  |  +90°  | none           | cubic disc=-59|
| 117 | m038    | 2.50266 | 2.550  |  -45°  | tr(a)=tr(b)    | unknown       |
| 177 | m032    | 2.68223 | 2.733  |  -45°  | tr(a)=tr(ab)   | unknown       |
| 209 | m206    | 2.82812 | 2.882  |  -60°  | tr(a)=-tr(b)   | Q(sqrt(-3))   |

## Key Results

### D(w^3) = v0 CONFIRMED
- D(w^3) = 0.98136883 = v0 exactly (diff < 1e-16)
- Minimal polynomial of w^3: x^4 - 3x^3 + 3x^2 - x - 1
- Substitution x=1+t gives: t^4 + t^3 - 1 = 0
- So w^3 - 1 = -u1 where u1 = 1-w^3 is the fundamental unit
- Minimal polynomial of u1 = 1-w^3: x^4 - x^3 - 1
- N(u1) = -1 exactly (unit of norm -1)
- Q(w^3) = Q(w) = K (same field, different generator)

### Eisenstein Conjecture
Three of six H1=Z/5 manifolds have ITF = Q(sqrt(-3)):
  m003 (theta*=-180°), m004 (theta*=+120°), m206 (theta*=-60°)
All three have theta* as a multiple of 60° (Eisenstein angles).
Conjecture: theta* in 60*Z <=> ITF contains Q(sqrt(-3))

### m004 ITF Confirmation
- ITF gen satisfies z^2 + z + 7 = 0 (residual 8.26e-14 = EXACT)
- disc = -27 = -3^3, field disc = -3 = Q(sqrt(-3))
- Second gen satisfies z^2 + 8z + 28 = 0, disc = -48 = -3*16
- Both generators confirm Q(sqrt(-3)) ancestry

### Symmetry Classification (exact, machine precision)
- m038: tr(a) = tr(b)    [Z/2 swap symmetry]
- m032: tr(a) = tr(ab)   [absorption symmetry]  
- m206: tr(a) = -tr(b)   [Z/2 anti-conjugation]
- m003, m004, m006: no special trace symmetry

### Next Steps
1. Prove the Eisenstein Conjecture
2. Find ITF of m038 and m032 (higher degree search needed)
3. Compute mixing matrix for m004 (theta*=+120°, order 3)
4. Extend census scan to idx 300-1000
5. D(w^3)=v0 algebraic proof: connect N(u1)=-1 to Bloch group primitivity
