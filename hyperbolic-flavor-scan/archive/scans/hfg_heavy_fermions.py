"""
Heavy fermion geodesic search in M_PMNS and M_CKM
Checks tau (q=68, ell=8.18), c (q=65, ell=7.82),
b (q=75, ell=9.02) at cutoff=10.0

Also checks: does the phi-lattice enrichment improve
when we look at the FULL spectrum including long geodesics?
"""
import snappy, numpy as np, warnings
warnings.filterwarnings('ignore')

PHI = (1+5**0.5)/2
LOG_PHI = np.log(PHI)
m_e = 0.000511

# All SM fermions with their phi-lattice assignments
all_masses = {
    'e':   (0.000511, 0),
    'u':   (0.00216,  12),
    'd':   (0.00467,  18),
    's':   (0.0934,   43),
    'mu':  (0.10566,  44),
    'c':   (1.275,    65),
    'tau': (1.77686,  68),
    'b':   (4.18,     75),
    't':   (172.76,   106),
}

M_PMNS = snappy.OrientableClosedCensus[1]
M_CKM  = snappy.OrientableClosedCensus[43]

print("=" * 70)
print("HEAVY FERMION GEODESIC SEARCH (cutoff=10.0)")
print("=" * 70)
print()
print("Target geodesic lengths:")
for name, (mass, q) in all_masses.items():
    ell = q/4 * LOG_PHI
    print(f"  {name:>4}: q={q:>4}, ell={ell:.5f}, "
          f"phi^(q/4)={PHI**(q/4):.2f}")
print()

print("Fetching PMNS spectrum (cutoff=10.0, this takes ~2-3 min)...")
ls_pmns = M_PMNS.length_spectrum(cutoff=10.0)
L_PMNS = []
for g in ls_pmns:
    try:
        v = g.length
        if callable(v): v = v()
        L_PMNS.append(float(complex(v).real))
    except: pass
L_PMNS.sort()
print(f"  Got {len(L_PMNS)} geodesics")
print()

print("Fetching CKM spectrum (cutoff=8.0)...")
ls_ckm = M_CKM.length_spectrum(cutoff=8.0)
L_CKM = []
for g in ls_ckm:
    try:
        v = g.length
        if callable(v): v = v()
        L_CKM.append(float(complex(v).real))
    except: pass
L_CKM.sort()
print(f"  Got {len(L_CKM)} geodesics")
print()

TOL = 0.008

print("=" * 70)
print(f"FULL MASS GEODESIC MATCHING (tol={TOL})")
print("=" * 70)
print()
print(f"{'name':>5} {'q':>5} {'ell_target':>10}  "
      f"{'PMNS':>10} {'dist':>7} {'hit':>3}  "
      f"{'CKM':>10} {'dist':>7} {'hit':>3}")
print("-"*70)

hits_P = hits_C = 0
in_range_P = in_range_C = 0

for name, (mass, q) in all_masses.items():
    ell_t = q/4 * LOG_PHI
    if ell_t <= 0:
        continue

    # PMNS
    if ell_t <= max(L_PMNS, default=0) + 0.1:
        in_range_P += 1
        dP = [(abs(l-ell_t), l) for l in L_PMNS]
        pD, pL = min(dP) if dP else (999, 0)
        pH = pD < TOL
        if pH: hits_P += 1
        pmns_str = f"{pL:.5f}" if pD < 5 else "---"
        pmns_dist = f"{pD:.5f}" if pD < 5 else "---"
    else:
        pmns_str = "(above cutoff)"
        pmns_dist = "---"
        pH = False

    # CKM
    if ell_t <= max(L_CKM, default=0) + 0.1:
        in_range_C += 1
        dC = [(abs(l-ell_t), l) for l in L_CKM]
        cD, cL = min(dC) if dC else (999, 0)
        cH = cD < TOL
        if cH: hits_C += 1
        ckm_str = f"{cL:.5f}" if cD < 5 else "---"
        ckm_dist = f"{cD:.5f}" if cD < 5 else "---"
    else:
        ckm_str = "(above cutoff)"
        ckm_dist = "---"
        cH = False

    ph_mark = "✓" if pH else " "
    ch_mark = "✓" if cH else " "
    print(f"{name:>5} {q:>5} {ell_t:>10.5f}  "
          f"{pmns_str:>10} {pmns_dist:>7} {ph_mark:>3}  "
          f"{ckm_str:>10} {ckm_dist:>7} {ch_mark:>3}")

print()
print(f"PMNS: {hits_P} hits / {in_range_P} in range  "
      f"({100*hits_P/max(in_range_P,1):.0f}%)")
print(f"CKM:  {hits_C} hits / {in_range_C} in range  "
      f"({100*hits_C/max(in_range_C,1):.0f}%)")

# ── Precision check on any new hits ─────────────────────────────────────────
print()
print("=" * 70)
print("PRECISION CHECK: ell/log(phi) for all matches")
print("=" * 70)
print()
print(f"{'name':>5} {'q':>5} {'ell':>10} {'ell/logphi':>12} "
      f"{'target_ratio':>14} {'dist_to_int':>13}")
print("-"*65)

for name, (mass, q) in all_masses.items():
    ell_t = q/4 * LOG_PHI
    target_ratio = q/4  # should be integer/4

    for label, lengths in [("PMNS", L_PMNS), ("CKM", L_CKM)]:
        if not lengths: continue
        dists = [(abs(l-ell_t), l) for l in lengths]
        best_d, best_l = min(dists)
        if best_d < 0.05:  # within 5% of target
            ratio = best_l / LOG_PHI
            dist_to_int = abs(ratio - round(ratio*4)/4)
            print(f"{name:>5} {q:>5} {best_l:>10.6f} {ratio:>12.6f} "
                  f"{target_ratio:>14.4f} {dist_to_int:>13.6f}  [{label}]")

# ── The phi-lattice enrichment at higher cutoff ──────────────────────────────
print()
print("=" * 70)
print("PHI-LATTICE ENRICHMENT (full spectrum)")
print("=" * 70)
print()

for label, lengths, cutoff in [("PMNS", L_PMNS, 10.0),
                                ("CKM",  L_CKM,  8.0)]:
    if not lengths: continue
    max_ell = max(lengths)
    lattice = [(k, k*LOG_PHI/4) for k in range(1, 200)
               if 0 < k*LOG_PHI/4 <= max_ell]
    hits = sum(1 for k, ell in lattice
               if any(abs(l-ell) < TOL for l in lengths))
    n_geod = len(lengths)
    coverage = n_geod * 2 * TOL / max_ell
    exp_rand = coverage * len(lattice)
    enrichment = hits / max(exp_rand, 0.1)

    print(f"{label} (cutoff={cutoff}):")
    print(f"  Geodesics: {n_geod}")
    print(f"  Phi-lattice points in [0,{max_ell:.1f}]: {len(lattice)}")
    print(f"  Matches (tol={TOL}): {hits}  "
          f"(expected random: {exp_rand:.1f})")
    print(f"  Enrichment: {enrichment:.3f}x  "
          f"({'OVER' if enrichment>1 else 'UNDER'}-represented)")
    print()

# ── The muon geodesic in context ────────────────────────────────────────────
print("=" * 70)
print("THE MUON GEODESIC IN CONTEXT")
print("=" * 70)
print()
print("Confirmed: PMNS contains geodesic at ell=5.29305")
print(f"  ell/log(phi) = 10.99942  (dist to 11: 0.000580)")
print(f"  Relative error from 11*log(phi): 0.00527%")
print()
print("Is this the closest ANY manifold gets to 11*log(phi)?")
print("Checking neighboring census manifolds...")

for idx in [0, 2, 3, 4, 5, 10, 11, 20, 43]:
    try:
        Mn = snappy.OrientableClosedCensus[idx]
        ls = Mn.length_spectrum(cutoff=5.5)
        lengths_n = []
        for g in ls:
            try:
                v = g.length
                if callable(v): v=v()
                lengths_n.append(float(complex(v).real))
            except: pass
        if not lengths_n: continue
        target = 11 * LOG_PHI
        dists = [(abs(l-target), l) for l in lengths_n]
        best_d, best_l = min(dists)
        ratio = best_l / LOG_PHI
        print(f"  idx={idx:>3} ({Mn.name():<12}) H1={str(Mn.homology()):<12} "
              f"nearest={best_l:.5f}  dist={best_d:.5f}  "
              f"ratio={ratio:.4f}")
    except: pass

print()
print("If idx=1 (M_PMNS) has the smallest distance to 11*log(phi)")
print("among all nearby manifolds, that's geometric significance.")
