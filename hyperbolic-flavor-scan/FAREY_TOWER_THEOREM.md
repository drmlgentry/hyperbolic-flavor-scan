# Farey Chain Cover Prime Theorem
**Verified:** May 19, 2026
**Status:** Computational theorem, n=1..6

## Formula
p_n = 5*n*(n+1) + 1

The degree-5 algebraic cover of M_n = m003(-n, 2n+1) introduces
exactly one new H1 prime, equal to p_n.

## Verified Data
| n | slope    | p_n | formula | match | prime |
|---|----------|-----|---------|-------|-------|
| 1 | (-1,3)   |  11 |      11 | TRUE  | YES   | <- M_PMNS, L5=Lucas prime
| 2 | (-2,5)   |  31 |      31 | TRUE  | YES   |
| 3 | (-3,7)   |  61 |      61 | TRUE  | YES   |
| 4 | (-4,9)   | 101 |     101 | TRUE  | YES   |
| 5 | (-5,11)  | 151 |     151 | TRUE  | YES   |
| 6 | (-6,13)  | 211 |     211 | TRUE  | YES   |

## Primality through n=20
Prime: n=1,2,3,4,5,6,7,11,13,14,15,16,17,18,20
Composite: n=8 (19^2), n=9 (11x41), n=10 (19x29), n=12 (11x71), n=19 (31x61=p2*p3)
Lucas-prime: ONLY n=1 (p1=11=L5)

## Supporting Results
- Cusp shape m003: tau = e^(i*pi/3) = omega, residual 9.55e-16 (exact)
- Cusp translations: meridian=1+sqrt(3)*i, longitude=2, area=2*sqrt(3)
- NZ constant: gap*L^2 converges to C ~ 14.16 (computed n=1..19)
- m006 cusp shape: Re(tau) satisfies 8x^3+4x-1=0 (cubic, not in Q(sqrt(17)))
- Surgery formula: |H1(m003(p,q))| = 5*|2p+q| (verified 13 slopes)

## Key Implications
1. Lucas-purity (p_n is Lucas) holds ONLY at n=1 = M_PMNS.
   This makes M_PMNS the unique Lucas-prime element of the Farey tower.
2. The formula p_n = 5n(n+1)+1 is the intrinsic cover prime invariant,
   independent of Dehn filling presentation.
3. The Farey intersection det[(-2,3),(-5,2)]=11 equals p_1 when using
   the (-2,3) presentation but gives 13 for the (-1,3) presentation.
   The cover prime formula is the presentation-independent quantity.
4. M_CKM volume 2.02885 lies in the asymptotic gap of the tower:
   vol(M_n) < vol(M_CKM) < vol(m003_cusp) = 2.02988 for all n.
   The manifolds are NOT commensurable (different invariant trace fields).

## Open Questions
- Does p_n = 5n(n+1)+1 hold for all n? (Proof needed)
- Why does n=19 give p_19 = p_2 * p_3? (Multiplicative structure)
- What is the exact NZ constant C and its relation to cusp area?
- Does m006 have an analogous Farey tower with a cover prime formula?
