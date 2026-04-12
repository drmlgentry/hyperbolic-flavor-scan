## Qubit analysis: final status

### CKM (K-factor, m006) -- COMPLETE
- Method: logm/Pauli axis + Gaussian kernel sigma=0.488 + QR
- Fitness: 0.01729 for words {aaB, AbA, AAb}
- All three words: same SU(2) conjugacy class, J=0 forced
- sigma_optimal = 0.488, log(phi) = 0.481, ratio = 1.014
- Qubit picture: three Bloch vectors at fixed latitude, sigma sets coherence length

### PMNS (N-factor, m003) -- DIFFERENT CONSTRUCTION
- The Gaussian/QR method gives fitness=0.127 at best (words: ab, bA, bba)
- J=0 regardless of Z/5 phase assignment (QR always gives real-like U)
- PMNS in the original paper uses Borel/Iwasawa N-factor, not Gaussian kernel
- These are genuinely different geometric constructions for the two sectors
- Qubit interpretation of N-factor requires separate analysis

### Key theorem confirmed
- Same SU(2) conjugacy class => J=0 (verified 1000 random trials)
- CKM: forced J=0 by topology (all words same conjugacy class)
- PMNS: J!=0 requires mixing of conjugacy classes (different construction)

### Open question
- Is sigma = log(phi) exact or approximate?
- Current answer: approximate (1.4% off), not a theorem
- Physical interpretation: coherence length ~ mass lattice spacing

## sigma = l2 result: CKM only

CKM: sigma_opt = 0.488, l2(m006) = 0.491, ratio = 1.005 -- CONFIRMED
PMNS: sigma_opt = 0.500, l2(m003) = 0.722, ratio = 1.44 -- NO MATCH

The l2 = sigma derivation is specific to CKM sector.
PMNS Gaussian method is wrong construction (word axes nearly degenerate:
bA-bba angle = 6.3 deg). PMNS requires N-factor/Borel construction.

PMNS axis angles: ab-bA=36.8, ab-bba=32.7, bA-bba=6.3 deg (degenerate)
CKM axis angles:  aaB-AbA=48.2, aaB-AAb=77.5, AbA-AAb=68.4 deg (well-separated)

## Covering tower mixing matrix scan

Gaussian/QR method, sigma=0.488 (CKM) or 0.500 (PMNS):
  m006 base (CKM, idx=43):   fitness=0.01673  -- EXCELLENT
  m003 base (PMNS, idx=1):   fitness=0.13299  -- POOR
  m003 deg-2 (idx=39):       fitness=0.16730  -- WORSE
  m003 deg-3 (idx=238):      fitness=0.13036  -- MARGINAL

Conclusions:
1. Covers do not improve PMNS fit -- no generation correspondence here
2. Gaussian/QR wrong construction for PMNS regardless of manifold
3. CKM result (0.017) is robust -- also found with {bAA,AbA,AAb}
4. The two sectors genuinely require different geometric constructions
