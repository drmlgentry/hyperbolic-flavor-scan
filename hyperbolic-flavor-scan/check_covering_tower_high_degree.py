#!/usr/bin/env sage
# check_covering_tower_high_degree.py
# Extend covering tower scan for M_PMNS to degree 12 using low_index_subgroups

import snappy
import json
import os
import re
from sympy import factorint

M = snappy.OrientableClosedCensus[1]   # PMNS manifold
MAX_INDEX = 12
CHECKPOINT = "pmns_covering_tower_high.json"

# Lucas primes (for reference)
LUCAS_PRIMES = {2, 3, 7, 11, 29, 47, 199, 521}
BASE_PRIMES = {5}   # base homology Z/5

def primes_from_homology(h_str):
    pset = set()
    for m in re.findall(r'Z/(\d+)', h_str):
        n = int(m)
        for f, _ in factorint(n).items():
            pset.add(f)
    return pset

# Load existing progress
if os.path.exists(CHECKPOINT):
    with open(CHECKPOINT, 'r') as f:
        state = json.load(f)
    print(f"Resuming from checkpoint: indices {list(state['completed'].keys())}")
else:
    state = {'completed': {}, 'all_primes': [], 'non_lucas': []}
    print("Starting fresh")

all_primes = set(state.get('all_primes', []))
non_lucas = set(state.get('non_lucas', []))

for idx in range(2, MAX_INDEX+1):
    if str(idx) in state['completed']:
        print(f"Index {idx}: already processed")
        continue
    print(f"Processing index {idx}...")
    try:
        subs = M.low_index_subgroups(idx)
        new_primes = set()
        for sub in subs:
            try:
                h_str = str(sub.homology())
                primes = primes_from_homology(h_str)
                # Exclude base primes (only 5)
                primes -= BASE_PRIMES
                new_primes.update(primes)
            except:
                continue
        # Record new primes (not already seen)
        truly_new = new_primes - all_primes
        all_primes.update(truly_new)
        non_lucas.update(p for p in truly_new if p not in LUCAS_PRIMES)
        state['completed'][str(idx)] = {
            'new_primes': list(truly_new),
            'all_primes_so_far': list(all_primes)
        }
        state['all_primes'] = list(all_primes)
        state['non_lucas'] = list(non_lucas)
        print(f"  New primes at index {idx}: {sorted(truly_new)}")
        # Save checkpoint after each index
        with open(CHECKPOINT, 'w') as f:
            json.dump(state, f, indent=2)
    except Exception as e:
        print(f"  ERROR at index {idx}: {e}")
        break

print("\n" + "="*60)
print(f"RESULTS for M_PMNS up to index {MAX_INDEX}:")
print(f"All new torsion primes: {sorted(all_primes)}")
print(f"Non‑Lucas primes: {sorted(non_lucas)}")
if non_lucas:
    print("Lucas‑purity: NO")
else:
    print("Lucas‑purity: YES (so far)")
print("="*60)