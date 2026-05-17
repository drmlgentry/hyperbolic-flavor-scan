"""
tower_deep.py
Extended covering tower for M_PMNS, degrees 20-24.
Saves checkpoint after each degree so progress is not lost.
Run: conda run -n sage python tower_deep.py
"""
import snappy, re, time, json, os

def full_primes(n):
    primes = set()
    temp = n
    p = 2
    while p*p <= temp:
        if temp % p == 0:
            primes.add(p)
            while temp % p == 0:
                temp //= p
        p += 1
    if temp > 1:
        primes.add(temp)
    return primes

CHECKPOINT_FILE = 'tower_checkpoint.json'
LUCAS = {2,3,7,11,29,47,199}

M = snappy.OrientableClosedCensus[1]
print(f"M_PMNS: {M.name()}, H1={M.homology()}")
print()

# Load checkpoint if exists
if os.path.exists(CHECKPOINT_FILE):
    with open(CHECKPOINT_FILE) as f:
        ckpt = json.load(f)
    seen = set(ckpt['seen'])
    completed = ckpt['completed']
    results = ckpt['results']
    print(f"Resuming from checkpoint: completed degrees {completed}")
    print(f"Known new primes so far: {sorted(seen)}")
else:
    seen = {5}  # base prime
    completed = list(range(2, 20))  # already computed through 19
    seen = {2, 3, 5, 11}  # confirmed through deg 19
    results = {
        5: [11], 7: [2], 8: [3]
    }
    print("Starting from degree 20 (degrees 2-19 already computed)")
    print("Known new primes through degree 19: {2, 3, 11}")

print()
START_DEG = max(completed) + 1 if completed else 2
MAX_DEG = 24
base = {5}

for deg in range(START_DEG, MAX_DEG + 1):
    print(f"deg={deg}: computing covers...", flush=True)
    t0 = time.time()

    try:
        covers = M.covers(deg)
        n_covers = len(covers)
        t_enum = time.time() - t0
        print(f"  {n_covers} covers found in {t_enum:.1f}s", flush=True)

        new_at_deg = set()
        for C in covers:
            h1 = str(C.homology())
            if h1 not in ('Z/5', '0', 'trivial'):
                print(f"  H1={h1}", flush=True)
            for m in re.findall(r'Z/(\d+)', h1):
                new_at_deg |= (full_primes(int(m)) - base - seen)

        seen |= new_at_deg
        elapsed = time.time() - t0

        if new_at_deg:
            nl = [p for p in new_at_deg if p not in LUCAS]
            results[deg] = sorted(new_at_deg)
            print(f"  *** NEW PRIMES: {sorted(new_at_deg)}, "
                  f"non-Lucas={nl}, t={elapsed:.1f}s", flush=True)
            if nl:
                print(f"  !!! NON-LUCAS PRIME FOUND: {nl} !!!")
        else:
            results[deg] = []
            print(f"  (none new), t={elapsed:.1f}s", flush=True)

    except Exception as e:
        print(f"  ERROR: {e}")
        elapsed = 0

    completed.append(deg)

    # Save checkpoint
    ckpt = {
        'seen': sorted(seen - base),
        'completed': completed,
        'results': {str(k):v for k,v in results.items()}
    }
    with open(CHECKPOINT_FILE, 'w') as f:
        json.dump(ckpt, f, indent=2)
    print(f"  [checkpoint saved]", flush=True)

print()
print("="*55)
print("FINAL SUMMARY")
print("="*55)
new_primes = sorted(seen - base)
print(f"All new primes through degree {MAX_DEG}: {new_primes}")
non_lucas = [p for p in new_primes if p not in LUCAS]
print(f"Non-Lucas: {non_lucas}")
if 7 in new_primes: print(f"  7 (L_4) FOUND")
else: print(f"  7 NOT found through degree {MAX_DEG}")
if 29 in new_primes: print(f"  29 (L_7) FOUND")
else: print(f"  29 NOT found through degree {MAX_DEG}")
if not non_lucas:
    print("Lucas-purity holds.")
else:
    print("LUCAS PURITY VIOLATED.")
