# Spectral-Arithmetic Bridge — Complete Results
## Date: 2026-04-25

## The Bridge Statement
A prime p from the Dehn filling slope arithmetic satisfying p ~ phi^q
corresponds to a geodesic of length ~ q*log(phi) in the flavor manifolds.

## Complete correspondence table

| Target | p | q | m003 delta | m006 delta | Source |
|---|---|---|---|---|---|
| phi^(3/2) ~ 2 | 2 | 1.5 | **0.03%** *** | 30.5% | |p_PMNS| (local) |
| phi^2 ~ 3 | 3 | 2.0 | 3.73% * | 2.16% * | |q_PMNS| |
| phi^3 ~ 5 | 5 | 3.0 | 11.5% | **1.47%** ** | torsion (local to CKM) |
| phi^4 ~ 7 | 7 | 4.0 | 5.33% | 3.77% | coord sum (weak) |
| phi^5 ~ 11 | 11 | 5.0 | **0.55%** ** | **1.07%** ** | det (inter-sector) |
| phi^7 ~ 29 | 29 | 7.0 | **0.30%** *** | **0.22%** *** | CKM norm^2 (inter-sector) |

## Asymmetry structure
- phi^(3/2) ~ 2: present in m003 (0.03%) but absent in m006 -- PMNS slope
  coordinate encoded locally in PMNS manifold only
- phi^3 ~ 5: present in m006 (1.47%) but absent in m003 -- torsion order
  encoded locally in CKM manifold only
- phi^5 ~ 11: present in BOTH (0.55%, 1.07%) -- inter-sector determinant
- phi^7 ~ 29: present in BOTH (0.30%, 0.22%) -- inter-sector CKM norm^2

## Interpretation
Each manifold encodes:
1. Its OWN slope coordinate in its geodesic spectrum (local information)
2. The SHARED inter-sector invariants det=11 and CKM norm^2=29 (shared info)

This is the structure of a communication channel:
  local state + shared channel invariants

## Null test
p=0.0137 (significant) for three slope invariants on phi-lattice.
N=9580 random coprime slope pairs.

## Connection to covering tower
- Prime 11 ~ phi^5: activates at degree 2, geodesic at log(11) in both
- Prime 29 ~ phi^7: activates at degree 5, geodesic at log(29) in both
- Activation degree ~ q/pi is approximate; exact connection TBD

## Next steps
1. Test generalization: do other cusped manifold pairs show same pattern?
2. Prove: for arithmetic hyperbolic 3-manifold M=m(p,q), geodesic lengths
   include log|Norm(p+q*omega)| where omega = fundamental unit of Q(sqrt(5))
3. Connect to Ihara zeta function: zeros at u=exp(-ell) ~ phi^{-q}
4. The phi^4 ~ 7 signal is weak (3-5%) -- check if 7 is truly absent
   or if the relevant geodesic is beyond cutoff=4.5
