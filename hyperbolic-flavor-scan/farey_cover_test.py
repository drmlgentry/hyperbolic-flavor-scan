"""
farey_cover_test.py
===================
Tests whether the Farey intersection of two Dehn filling slopes equals
the first new prime introduced in the algebraic covering tower.

Specifically, for a Dehn filling M(p,q) of a cusped manifold, we test:

    det[(p1,p2),(q1,q2)] = first new prime in M.covers(d) for d=2..7?

This is the key conjecture underlying the Lucas structure paper.
The HFG claim: for M_PMNS = m003(-2,3) paired with M_CKM slope (-5,2):
    det[(-2,-5),(3,2)] = -4+15 = 11 = L5 (5th Lucas prime)
    and 11 is the ONLY new prime through degree 6 in M_PMNS.covers()

We test this across a range of Dehn fillings of the m003 cusped manifold.

Run with:
    conda run -n sage python farey_cover_test.py 2>/dev/null

Output: printed table + results saved to data/farey_cover_results.csv
"""

import os
import sys
import csv
from sympy import factorint, isprime

try:
    import snappy
except ImportError:
    print("ERROR: snappy not found. Run in: conda run -n sage python farey_cover_test.py")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

CUSPED_NAME = "m003"          # The cusped manifold (one cusp)
REFERENCE_SLOPE = (-2, 3)     # M_PMNS slope (the "base" filling)
CKM_SLOPE = (-5, 2)           # M_CKM slope (the "paired" filling)
MAX_COVER_DEGREE = 6          # Search covers through this degree
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# Test slopes: a range of Dehn fillings to probe the conjecture
# Include the actual HFG slopes plus neighbors
TEST_SLOPES = [
    (-1, 2),
    (-1, 3),
    (-2, 3),   # M_PMNS — the HFG case
    (-2, 5),
    (-3, 4),
    (-3, 5),
    (-4, 5),
    (-5, 2),   # M_CKM slope (as a single filling of m003, for reference)
    (-1, 4),
    (-3, 7),
    (-4, 7),
    ( 1, 2),
    ( 1, 3),
    ( 2, 3),
    ( 2, 5),
    ( 3, 4),
    ( 3, 5),
    ( 4, 5),
]

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def farey_det(slope1, slope2):
    """
    Compute det[(p1,p2),(q1,q2)] = p1*q2 - p2*q1
    This is the Farey intersection / surgery intersection number.
    """
    p1, q1 = slope1
    p2, q2 = slope2
    return abs(p1 * q2 - p2 * q1)


def get_h1_primes(manifold):
    """
    Extract the prime factors of the torsion part of H1(M; Z).
    Returns a set of primes.
    """
    try:
        h1 = manifold.homology()
        # h1.betti_number() gives free rank; h1.order() gives torsion order
        # For torsion groups like Z/n, order() = n
        order = h1.order()
        if order == 0:
            return set()  # infinite order (free part)
        if order == 1:
            return set()  # trivial
        return set(factorint(order).keys())
    except Exception as e:
        return None


def new_primes_in_cover(base_primes, cover_manifold):
    """
    Find primes in cover's H1 that are not in the base manifold's H1.
    """
    cover_primes = get_h1_primes(cover_manifold)
    if cover_primes is None:
        return None
    return cover_primes - base_primes


def analyze_filling(cusped, slope, reference_slope):
    """
    Fill the cusped manifold at `slope`, find algebraic covers,
    extract new H1 primes, and compare with Farey det.
    """
    p, q = slope
    result = {
        "slope": f"({p},{q})",
        "farey_det": None,
        "base_h1": None,
        "base_h1_primes": None,
        "covers_found": {},
        "first_new_prime": None,
        "farey_equals_cover_prime": None,
        "note": ""
    }

    try:
        M = cusped.copy()
        M.dehn_fill(slope)
        M = M.filled_manifold()

        # Get base H1
        base_primes = get_h1_primes(M)
        result["base_h1"] = str(M.homology())
        result["base_h1_primes"] = sorted(base_primes) if base_primes is not None else []

        # Farey det against reference slope
        fd = farey_det(slope, reference_slope)
        result["farey_det"] = fd

        # Search for covers degree 2..MAX_COVER_DEGREE
        all_new_primes = set()
        first_new_prime = None
        first_new_degree = None

        for deg in range(2, MAX_COVER_DEGREE + 1):
            try:
                covers = M.covers(deg)
                deg_new = []
                for cov in covers:
                    np_set = new_primes_in_cover(base_primes, cov)
                    if np_set:
                        deg_new.extend(sorted(np_set))
                        all_new_primes.update(np_set)
                        if first_new_prime is None:
                            first_new_prime = min(np_set)
                            first_new_degree = deg
                result["covers_found"][deg] = sorted(set(deg_new))
            except Exception as e:
                result["covers_found"][deg] = f"ERROR: {e}"

        result["first_new_prime"] = first_new_prime
        result["first_new_degree"] = first_new_degree
        result["all_new_primes"] = sorted(all_new_primes)

        if first_new_prime is not None and fd is not None:
            result["farey_equals_cover_prime"] = (fd == first_new_prime)
        else:
            result["farey_equals_cover_prime"] = None

    except Exception as e:
        result["note"] = f"FAILED: {e}"

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 72)
    print("FAREY / COVER PRIME CONJECTURE TEST")
    print(f"Cusped manifold: {CUSPED_NAME}")
    print(f"Reference slope: {REFERENCE_SLOPE}  (M_PMNS)")
    print(f"CKM slope:       {CKM_SLOPE}        (used for Farey det)")
    print(f"Max cover degree: {MAX_COVER_DEGREE}")
    print("=" * 72)

    try:
        cusped = snappy.Manifold(CUSPED_NAME)
        print(f"Loaded: {cusped.name()}  cusps={cusped.num_cusps()}\n")
    except Exception as e:
        print(f"ERROR loading cusped manifold: {e}")
        sys.exit(1)

    # HFG reference case first
    print("--- HFG reference case ---")
    print(f"Farey det[{REFERENCE_SLOPE}, {CKM_SLOPE}] = {farey_det(REFERENCE_SLOPE, CKM_SLOPE)}")
    print()

    rows = []
    for slope in TEST_SLOPES:
        print(f"Processing slope {slope} ...", flush=True)
        res = analyze_filling(cusped, slope, REFERENCE_SLOPE)
        rows.append(res)

        # Print compact summary line
        fd = res["farey_det"]
        fnp = res["first_new_prime"]
        fnd = res.get("first_new_degree")
        match = res["farey_equals_cover_prime"]
        all_np = res.get("all_new_primes", [])

        if res["note"]:
            print(f"  {res['slope']:10s}  {res['note']}")
        else:
            match_str = "MATCH" if match else ("MISS " if match is False else "N/A  ")
            print(f"  slope={res['slope']:10s}  "
                  f"H1={res['base_h1']:12s}  "
                  f"Farey={fd:4d}  "
                  f"first_new_prime={str(fnp):6s} (deg {fnd})  "
                  f"all_new={all_np}  "
                  f"[{match_str}]")

    # Summary table
    print()
    print("=" * 72)
    print("SUMMARY TABLE")
    print(f"{'Slope':12s} {'H1':12s} {'Farey':6s} {'FirstPrime':12s} {'Deg':4s} {'Match':6s} {'AllNewPrimes'}")
    print("-" * 72)
    matches = 0
    misses = 0
    na = 0
    for r in rows:
        if r["note"]:
            print(f"{r['slope']:12s} FAILED")
            continue
        fd = r["farey_det"]
        fnp = r["first_new_prime"]
        fnd = r.get("first_new_degree", "?")
        match = r["farey_equals_cover_prime"]
        all_np = r.get("all_new_primes", [])
        match_str = "YES" if match else ("NO" if match is False else "N/A")
        if match is True:
            matches += 1
        elif match is False:
            misses += 1
        else:
            na += 1
        print(f"{r['slope']:12s} {r['base_h1']:12s} {str(fd):6s} "
              f"{str(fnp):12s} {str(fnd):4s} {match_str:6s} {all_np}")

    print("-" * 72)
    print(f"MATCH: {matches}  MISS: {misses}  N/A (no cover prime found): {na}")
    print()

    # Interpretation
    print("INTERPRETATION:")
    if misses == 0 and matches > 0:
        print(f"  Conjecture SUPPORTED on all {matches} testable cases.")
        print("  Farey det = first new algebraic cover prime holds across these fillings.")
        print("  This warrants stating as a conjecture in the paper.")
    elif misses > 0:
        print(f"  Conjecture FAILS on {misses} case(s) — it is slope-specific, not general.")
        print("  The HFG claim should be stated as an observation, not a theorem.")
    else:
        print("  Insufficient data: no filling produced a new cover prime through degree",
              MAX_COVER_DEGREE)
        print("  Try increasing MAX_COVER_DEGREE.")

    # Save CSV
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    outfile = os.path.join(OUTPUT_DIR, "farey_cover_results.csv")
    with open(outfile, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["slope", "h1", "farey_det", "first_new_prime",
                         "first_new_degree", "all_new_primes", "match", "note"])
        for r in rows:
            writer.writerow([
                r["slope"],
                r.get("base_h1", ""),
                r.get("farey_det", ""),
                r.get("first_new_prime", ""),
                r.get("first_new_degree", ""),
                str(r.get("all_new_primes", [])),
                r.get("farey_equals_cover_prime", ""),
                r.get("note", "")
            ])
    print(f"\nResults saved to: {outfile}")
    print("=" * 72)


if __name__ == "__main__":
    main()
