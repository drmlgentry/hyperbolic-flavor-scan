# Spectral Gap Analysis — Final Assessment
## Date: 2026-04-21

## Method
Heat kernel trace via Selberg trace formula.
K(t) = vol * exp(-t) / (4*pi*t)^(3/2)
      + sum_gamma ell/(2*sinh(ell/2)) * exp(-t)/(4*pi*t)^(1/2) * exp(-ell^2/(4t))
Large-t linear fit gives mu_1 estimate.

## Results (279 manifolds, vol 0.94-6.36)
Single-exponential fit:  mu_1 - 1 = 0.00955 +/- 0.00004 (CV=0.385%)
Slope extrapolation:     mu_1 - 1 ~ 0.004 (more reliable, less biased)
Two-exponential fit:     mu_1 - 1 ~ 0.006 (unstable due to near-degeneracy)

True mu_1 in range [1.003, 1.007] -- method cannot resolve more precisely.

## Correct interpretation
mu_1 > 1 for all compact H^3/Gamma (strict inequality, well-known)
mu_1 -> 1 as vol -> inf (Benjamini-Schramm convergence to H^3)
Bottom of continuous spectrum of H^3 is exactly 1.

The tight clustering (CV=0.4%) reflects:
1. All census manifolds are small-volume with similar arithmetic structure
2. Benjamini-Schramm: at these volumes, mu_1 is pulled toward 1
3. NOT a new universal constant -- this is known spectral geometry

## What IS notable
- CV < 0.4% across 6x volume range
- Consistent with all manifolds having mu_1 near the LOWER bound of [1, 47.32]
- Small-volume arithmetic manifolds cluster near the Selberg analog bound

## What we cannot claim
- lambda_1 = sqrt(5)/2 exactly: gap is 0.17% but within method uncertainty
- lambda_1 = phi^(1/4): ruled out at 431 sigma
- A new universal constant: the clustering is real but explained by known theory

## Next step for exact values
Use Booker-Strombergsson method as in Bonifacio-Mazac-Pal (2025),
Commun. Math. Phys. 406:51, which computes exact t_1 for census manifolds.
Contact authors or implement their method.

## Connection to framework
The Dirac spectrum does NOT directly encode the phi-lattice.
The phi-lattice arises from geodesic lengths via Mahler measure,
not from Laplacian eigenvalues.
These are structurally distinct objects.
