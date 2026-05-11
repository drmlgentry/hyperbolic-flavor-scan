## Priority 3: CP Phase — CLOSED (May 10, 2026)

### All approaches tested and failed

1. Word twist angles (logm method): hits 200.6 deg in default triangulation
   but std=36 deg across 20 trials -- triangulation artifact
   
2. CS formula 2*pi*CS*m: cycles through {90,180,270,0} deg only
   Cannot produce 197 deg. CS=1/4 is too constrained.
   
3. Character variety k=3: gives 216 deg (diff=19 deg). Not close enough.

4. Canonical eigenvalue twist: std=68 deg -- worse than logm method.

5. Dedekind sum D(-2,3)*360 = 20 deg. Wrong.

6. Slope phase arg(-2+3i) = 123.7 deg. Wrong.

### Conclusion
delta_CP = 197 deg is NOT derivable from m003(-2,3) by any known
triangulation-independent method. The J_PMNS != 0 result (homology
no-go theorem) correctly predicts CP violation EXISTS in the lepton
sector but does not give the specific phase value.

The specific phase may require:
- A non-Abelian flat connection on m003 beyond PSL(2,C)
- Coupling between the two flavor manifolds m003 and m006
- Physics input beyond pure topology (e.g., Yukawa texture)

### Interesting near-miss
Character variety angle at k=3: 2*pi*3/5 = 216 deg (diff=19 deg)
Closest of any triangulation-independent quantity. Worth noting
in a footnote but not claimable as a derivation.

### Retrodiction that DOES work  
arctan(|-2|/|3|) = arctan(2/3) = 33.69 deg ~ theta_12 = 33.65 deg
The theta_12 PMNS mixing angle matches the Dehn filling slope ratio
to 0.04 deg. This is a genuine (if unexplained) coincidence worth
noting in future papers.
