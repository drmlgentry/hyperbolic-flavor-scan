# Length-7 Census Results Summary
Generated: 2026-03-16

## m003 (Meyerhoff manifold)
- 4,372 words evaluated
- Minimum phi_fold at length <= 7: **2.643 deg** (same as length-5 floor)
- No new small-angle geodesics found at lengths 6-7
- theta13_CKM (0.201 deg) not found at any length <= 7

## m006
- 4,372 words evaluated  
- **NEW: phi_fold = 0.005 deg at length 6** (word: aaabab, full family of 26)
  - phi = 179.9947 deg (near-identity: almost pure half-turn)
  - mod_lambda = 3.7297
- theta13_CKM = 0.201 deg falls in spectral gap [0.005, 1.611] at lengths <= 7
- No word with phi_fold in [0.05, 1.60] exists at any length <= 7

## Spectral gap structure on m006
The gap [0.005, 1.611] deg contains theta13_CKM = 0.201 deg.
This gap is likely arithmetic in origin (H1=Z/5 torsion structure).

## Implications for Paper 8
Section 5 (Null result) needs updating in referee response:
- Original claim: "theta13_CKM appears at word length ~12"  
- Revised: "theta13_CKM falls in a spectral gap [0.005, 1.61] at lengths <= 7;
  the gap origin is likely arithmetic and its width at longer lengths is unknown"
- The near-identity geodesic at 0.005 deg is a new Paper 9 result

## Near-identity geodesic: 5x torsion check
5 x 0.005298 = 0.0265 deg
5 x 1.687 = 8.435 deg ~ theta13_nu = 8.54 deg (1.2% error)
-- suggests torsion-5 structure may organize the spectral floor
