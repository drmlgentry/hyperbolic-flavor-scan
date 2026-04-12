=== Qubit/Bloch sphere analysis: summary of findings ===

1. CONFIRMED: logm/Pauli method reproduces CKM fitness=0.01729
   - Words: {aaB, AbA, AAb} on m006 (OrientableClosedCensus[43])
   - Method: matrix logarithm -> Pauli decomposition -> axis
   - Gaussian kernel sigma=0.49, QR orthogonalization
   - Axes confirmed: aaB=(0.237,-0.899,0.368) etc.

2. CONFIRMED: All three CKM words have identical SU(2) trace
   - Tr(U) = -0.1231 + 2.0108i for all three words
   - Same conjugacy class in SU(2)
   - J=0 numerically confirmed (1000 random trials)

3. CONFIRMED: sigma ~ log(phi) to 1.4%
   - Optimal sigma = 0.488
   - log(phi) = 0.481
   - Ratio = 1.014
   - NOT exact -- closest approximation, not a theorem

4. PMNS best fitness = 0.127 (words: ab, bA, bba on m003)
   - Worse than CKM (0.017) -- PMNS word search needs extension
   - The original PMNS paper used a different/longer word set

5. OPEN: exact sigma determination from geometry
   - Is sigma determined by the manifold invariants?
   - Candidate: sigma = log(phi) exactly (1.4% off)
   - Candidate: sigma from systole or spectral gap?
   - Requires independent derivation, not just numerical fit

Files:
  /tmp/ckm_correct.py    -- confirms fitness=0.017
  /tmp/sigma_fine.py     -- sigma scan, minimum at 0.488
  /tmp/sigma_logphi.py   -- log(phi) comparison
