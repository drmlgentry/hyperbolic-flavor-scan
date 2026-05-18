#!/usr/bin/env sage
# ckm_covering_tower.py
# Run: sage ckm_covering_tower.py

import snappy
import json
import os
from sympy import factorint
import re

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
MANIFOLD = snappy.OrientableClosedCensus[43]   # m006(-5,2)
MAX_DEG = 9
CHECKPOINT_FILE = "ckm_covering_tower.json"

# Lucas primes (for reference)
LUCAS_PRIMES = {2,3,7,11,29,47,199,521}

# ------------------------------------------------------------
# Helper: extract prime factors from homology string
# ------------------------------------------------------------
def primes_from_homology(h_str):
    primes = set()
    # Find all "Z/ n" patterns
    for n_str in re.findall(r'Z/(\d+)', h_str):
        n = int(n_str)
        if n <= 1:
            continue
        for p, _ in factorint(n).items():
            primes.add(p)
    return primes

# ------------------------------------------------------------
# Load or initialize checkpoint
# ------------------------------------------------------------
if os.path.exists(CHECKPOINT_FILE):
    with open(CHECKPOINT_FILE, 'r') as f:
        state = json.load(f)
    print(f"Resuming from checkpoint: completed degrees {list(state['completed'].keys())}")
else:
    state = {'completed': {}, 'all_primes': [], 'non_lucas': []}
    print("Starting fresh")

# ------------------------------------------------------------
# Scan degrees 2..MAX_DEG
# ------------------------------------------------------------
base_h = str(MANIFOLD.homology())
base_primes = primes_from_homology(base_h)
print(f"Base homology: {base_h}, base primes: {sorted(base_primes)}")
print()

all_primes = set(state.get('all_primes', []))
non_lucas = set(state.get('non_lucas', []))

for d in range(2, MAX_DEG+1):
    d_str = str(d)
    if d_str in state['completed']:
        # Use cached data
        cp = set(state['completed'][d_str]['primes'])
        nl = set(state['completed'][d_str]['non_lucas'])
        print(f"Degree {d}: cached – primes={sorted(cp)}, non‑Lucas={sorted(nl)}")
        all_primes.update(cp)
        non_lucas.update(nl)
        continue

    print(f"Degree {d}: computing covers...")
    deg_primes = set()
    deg_non_lucas = set()
    try:
        covers = MANIFOLD.covers(d)
        for i, cov in enumerate(covers):
            h = str(cov.homology())
            primes = primes_from_homology(h)
            # Remove base primes (so we only count new torsion)
            primes -= base_primes
            deg_primes.update(primes)
            for p in primes:
                if p not in LUCAS_PRIMES:
                    deg_non_lucas.add(p)
            if (i+1) % 50 == 0:
                print(f"  processed {i+1} covers so far...")
        # After processing all covers
        print(f"  New primes at degree {d}: {sorted(deg_primes)}")
        print(f"  Non‑Lucas: {sorted(deg_non_lucas)}")
        all_primes.update(deg_primes)
        non_lucas.update(deg_non_lucas)
        # Save checkpoint
        state['completed'][d_str] = {
            'primes': sorted(deg_primes),
            'non_lucas': sorted(deg_non_lucas)
        }
        state['all_primes'] = sorted(all_primes)
        state['non_lucas'] = sorted(non_lucas)
        with open(CHECKPOINT_FILE, 'w') as f:
            json.dump(state, f, indent=2)
    except Exception as e:
        print(f"  Error at degree {d}: {e}")
        # Still save what we have
        state['completed'][d_str] = {'error': str(e)}
        with open(CHECKPOINT_FILE, 'w') as f:
            json.dump(state, f, indent=2)
        break

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------
print("\n" + "="*60)
print("CKM COVERING TOWER SUMMARY (degrees 2..9)")
print("="*60)
print(f"All new torsion primes: {sorted(all_primes)}")
print(f"Non‑Lucas primes: {sorted(non_lucas)}")
if non_lucas:
    print("Lucas‑purity: NO")
else:
    print("Lucas‑purity: YES (observed up to degree 9)")
print("="*60)