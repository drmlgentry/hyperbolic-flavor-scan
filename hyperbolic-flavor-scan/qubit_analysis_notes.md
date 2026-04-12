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
