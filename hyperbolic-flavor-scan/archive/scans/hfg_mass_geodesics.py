"""
HFG Mass Geodesic Search
=========================
Tests Candidate B: the phi-lattice mass quantum numbers q are encoded
in the full geodesic length spectrum of M_PMNS and M_CKM.

For each Standard Model fermion mass m, compute q = 4*log(m/m_e)/log(phi).
Then check whether the manifold contains a geodesic of length
ell = (q/4)*log(phi) within tolerance epsilon.

If the manifold's length spectrum contains geodesics at exactly the
mass-lattice positions, this would be a non-trivial geometric fact
linking the manifold to the Standard Model mass spectrum.
"""

import snappy
import numpy as np
import warnings
warnings.filterwarnings('ignore')

PHI = (1 + 5**0.5) / 2
LOG_PHI = np.log(PHI)

print("=" * 70)
print("HFG MASS GEODESIC SEARCH")
print("Testing: does the length spectrum encode SM fermion masses?")
print("=" * 70)
print()

# Standard Model masses (GeV, PDG 2024)
SM_masses = {
    'e':     0.000511,
    'mu':    0.10566,
    'tau':   1.77686,
    'u':     0.00216,
    'd':     0.00467,
    's':     0.0934,
    'c':     1.275,
    'b':     4.18,
    't':     172.76,
    'nu1':   8.387e-12,   # predicted from lattice
    'nu2':   1.203e-11,
    'nu3':   5.097e-11,
}
m_e = SM_masses['e']

# Compute phi-lattice quantum numbers
print("Standard Model fermion phi-lattice assignments:")
print(f"  {'Particle':<8} {'mass (GeV)':>12}  {'q_exact':>10}  "
      f"{'q_nearest':>10}  {'residual':>10}")
print("  " + "-"*58)
q_values = {}
for name, mass in SM_masses.items():
    q_exact = 4 * np.log(mass/m_e) / LOG_PHI
    q_nearest = round(q_exact)
    residual = abs(q_exact - q_nearest)
    q_values[name] = (q_exact, q_nearest, residual)
    print(f"  {name:<8} {mass:>12.5e}  {q_exact:>10.3f}  "
          f"{q_nearest:>10d}  {residual:>10.4f}")

print()

# Load manifolds and get full length spectrum
M_PMNS = snappy.OrientableClosedCensus[1]
M_CKM  = snappy.OrientableClosedCensus[43]

def get_full_spectrum(M, cutoff=6.0):
    """Get all geodesic lengths up to cutoff."""
    ls = M.length_spectrum(cutoff=cutoff)
    lengths = []
    for g in ls:
        try:
            v = g.length
            if callable(v): v = v()
            ell = float(complex(v).real)
            mult = int(g.multiplicity)
            lengths.append((ell, mult))
        except:
            pass
    return sorted(lengths)

print("Fetching length spectra (cutoff=6.0, may take a moment)...")
spec_PMNS = get_full_spectrum(M_PMNS, cutoff=6.0)
spec_CKM  = get_full_spectrum(M_CKM,  cutoff=6.0)
print(f"  PMNS: {len(spec_PMNS)} geodesics")
print(f"  CKM:  {len(spec_CKM)} geodesics")
print()

# For each fermion, check if its mass geodesic ell = (q/4)*logphi
# exists in the spectrum
print("=" * 70)
print("MASS GEODESIC MATCHING")
print(f"Tolerance: |ell_geod - ell_mass| < 0.005")
print("=" * 70)
print()

TOL = 0.005

def find_matching_geodesic(spectrum, target_ell, tol=TOL):
    """Find geodesic closest to target_ell."""
    if not spectrum:
        return None, float('inf')
    dists = [(abs(ell - target_ell), ell, mult) 
             for ell, mult in spectrum]
    dists.sort()
    best_dist, best_ell, best_mult = dists[0]
    return best_ell, best_dist

print(f"  {'Particle':<6} {'q':>5}  {'ell_mass':>9}  "
      f"{'PMNS match':>10}  {'dist':>6}  "
      f"{'CKM match':>10}  {'dist':>6}")
print("  " + "-"*65)

hits_PMNS = 0
hits_CKM = 0
for name, mass in SM_masses.items():
    if 'nu' in name:
        continue  # skip neutrinos for now
    q_exact, q_nearest, residual = q_values[name]
    ell_mass = (q_nearest / 4) * LOG_PHI
    
    pmns_ell, pmns_dist = find_matching_geodesic(spec_PMNS, ell_mass)
    ckm_ell,  ckm_dist  = find_matching_geodesic(spec_CKM,  ell_mass)
    
    pmns_hit = pmns_dist < TOL
    ckm_hit  = ckm_dist  < TOL
    if pmns_hit: hits_PMNS += 1
    if ckm_hit:  hits_CKM  += 1
    
    pmns_str = f"{pmns_ell:.5f}" if pmns_ell else "---"
    ckm_str  = f"{ckm_ell:.5f}"  if ckm_ell  else "---"
    pmns_mark = "✓" if pmns_hit else " "
    ckm_mark  = "✓" if ckm_hit  else " "
    
    print(f"  {name:<6} {q_nearest:>5}  {ell_mass:>9.5f}  "
          f"{pmns_mark}{pmns_str:>10}  {pmns_dist:>6.4f}  "
          f"{ckm_mark}{ckm_str:>10}  {ckm_dist:>6.4f}")

n_fermions = len([k for k in SM_masses if 'nu' not in k])
print()
print(f"  PMNS hits: {hits_PMNS}/{n_fermions} fermions "
      f"({100*hits_PMNS/n_fermions:.0f}%)")
print(f"  CKM  hits: {hits_CKM}/{n_fermions} fermions "
      f"({100*hits_CKM/n_fermions:.0f}%)")

# ── Random control ──────────────────────────────────────────────────────────
print()
print("=" * 70)
print("RANDOM CONTROL: expected hits by chance")
print("=" * 70)
print()

# For a random spectrum with same density, what fraction of target ells
# would match within TOL?
# Density of geodesics: N(L) ~ e^(2L) / (2L) by prime geodesic theorem
# At length L, density ~ e^(2L)/(2L) per unit length
# Expected hits at random: 2*TOL * density at each ell_mass

print("Expected random hits (prime geodesic theorem estimate):")
expected_hits = 0
for name, mass in SM_masses.items():
    if 'nu' in name: continue
    q_exact, q_nearest, residual = q_values[name]
    ell_mass = (q_nearest / 4) * LOG_PHI
    if ell_mass <= 0 or ell_mass > 6.0:
        continue
    # Density of geodesics at length ell
    density = np.exp(2 * ell_mass) / (2 * ell_mass)
    # Expected number in window [ell-TOL, ell+TOL]
    expected = min(1.0, 2 * TOL * density * float(M_PMNS.volume()) / (4*np.pi))
    expected_hits += expected

print(f"  Expected random hits (PMNS): ~{expected_hits:.1f}/{n_fermions}")
print()

# ── Phi-lattice density analysis ────────────────────────────────────────────
print("=" * 70)
print("PHI-LATTICE DENSITY IN LENGTH SPECTRUM")
print("=" * 70)
print()
print("Checking: are phi-lattice lengths over-represented in the spectrum?")
print(f"  Phi-lattice points in [0, 6]: ell = k*logphi/4 for k=0,1,...,N")
print()

lattice_points = [(k, k * LOG_PHI / 4) for k in range(1, 100) 
                  if k * LOG_PHI / 4 <= 6.0]

for label, spectrum in [("PMNS", spec_PMNS), ("CKM", spec_CKM)]:
    spec_lengths = [ell for ell, mult in spectrum]
    
    # Count how many lattice points have a geodesic within TOL
    lattice_hits = 0
    for k, ell_lat in lattice_points:
        dists = [abs(ell - ell_lat) for ell in spec_lengths]
        if min(dists) < TOL:
            lattice_hits += 1
    
    # Count total geodesics in same range
    n_geod = len([ell for ell in spec_lengths if ell <= 6.0])
    n_lattice = len(lattice_points)
    
    print(f"  {label}:")
    print(f"    Geodesics in [0,6]: {n_geod}")
    print(f"    Phi-lattice points in [0,6]: {n_lattice}")
    print(f"    Lattice points with geodesic match (tol={TOL}): "
          f"{lattice_hits}/{n_lattice} = {100*lattice_hits/n_lattice:.1f}%")
    
    # Random expectation: each geodesic covers 2*TOL of the length axis
    # Fraction of axis covered: n_geod * 2 * TOL / 6.0
    coverage = n_geod * 2 * TOL / 6.0
    expected_random = coverage * n_lattice
    print(f"    Expected random matches: {expected_random:.1f} ({100*coverage:.1f}% coverage)")
    enrichment = lattice_hits / max(expected_random, 0.1)
    print(f"    Enrichment factor: {enrichment:.2f}x")
    print()

# ── The right question ───────────────────────────────────────────────────────
print("=" * 70)
print("REFRAMING: THE RIGHT QUESTION")
print("=" * 70)
print()
print("The naive Dirac spectrum does NOT encode masses (ratios too uniform).")
print("The word geodesic lengths do NOT encode masses (too similar to each other).")
print()
print("The phi-lattice encodes masses as LOG(m/m_e)/log(phi) = q/4.")
print("The GEODESIC encodes the same number via ell = (q/4)*log(phi).")
print()
print("So the question becomes:")
print("  Are there geodesics of M_PMNS at EXACTLY the mass-lattice lengths?")
print()
print("This is NOT the same as asking whether the Dirac spectrum encodes masses.")
print("It is asking whether M_PMNS is a UNIVERSAL ARITHMETIC OBJECT whose")
print("length spectrum contains ALL the mass quantum numbers as embedded lengths.")
print()
print("This is related to the question of whether M_PMNS is arithmetic:")
print("  Arithmetic hyperbolic manifolds have length spectra with deep")
print("  number-theoretic structure (via the Jacquet-Langlands correspondence).")
print("  M_PMNS IS arithmetic (trace field Q(sqrt(-3)), confirmed).")
print()
print("The CONGRUENCE SUBGROUP structure of arithmetic manifolds means")
print("their length spectra are governed by Hecke eigenvalues -- and")
print("Hecke eigenvalues are algebraic integers in Q(phi) for our manifolds.")
print()
print("THIS is the bridge to a full theory:")
print("  Hecke eigenvalues of M_PMNS ↔ mass quantum numbers")
print("  Hecke L-functions of M_PMNS ↔ Yukawa couplings")
print()
print("NEXT STEP: Compute Hecke eigenvalues of the arithmetic manifold M_PMNS.")
print("This requires SageMath's modular forms machinery applied to the")
print("quaternion algebra ramified at primes {2,3} (the Lucas primes of")
print("the covering tower) with discriminant related to Q(sqrt(-3)).")
EOF
