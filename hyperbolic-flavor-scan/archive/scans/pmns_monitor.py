import pandas as pd
import time
import os
import sys

RESULTS_FILE  = "pmns_scan_results.csv"
PROGRESS_FILE = "pmns_scan_progress.txt"
REFRESH       = 30  # seconds

PMNS = [
    [0.821, 0.550, 0.148],
    [0.357, 0.339, 0.871],
    [0.442, 0.762, 0.471]
]

def clear():
    os.system('cls' if os.name == 'nt' else 'clear')

def show():
    clear()
    print("=" * 70)
    print("  PMNS SCAN MONITOR")
    print("=" * 70)

    # Progress
    if os.path.exists(PROGRESS_FILE):
        try:
            line = open(PROGRESS_FILE).read().strip()
            idx, total, best, eta = line.split(',')
            pct = 100 * int(idx) / int(total)
            bar = '#' * int(pct/2) + '-' * (50 - int(pct/2))
            print(f"  Progress: [{bar}] {pct:.1f}%  ({idx}/{total})")
            print(f"  Global best fitness: {float(best):.6f}")
            print(f"  ETA: {float(eta)/60:.1f} min remaining")
        except Exception:
            print("  Waiting for scan to start...")
    else:
        print("  Waiting for scan to start...")

    print()
    print("  PMNS target:")
    for row in PMNS:
        print("  ", "  ".join(f"{v:.3f}" for v in row))
    print()

    # Top results
    if os.path.exists(RESULTS_FILE):
        try:
            df = pd.read_csv(RESULTS_FILE)
            top = df.nsmallest(10, 'fitness')
            print(f"  TOP 10 (of {len(df)} manifolds scanned so far):")
            print(f"  {'Rank':<5} {'Manifold':<10} {'Vol':<8} {'H1':<8} "
                  f"{'Words':<30} {'sigma':<7} {'Fitness':<10} {'Angles'}")
            print("  " + "-" * 100)
            for rank, (_, row) in enumerate(top.iterrows(), 1):
                words = f"{row['word1']}/{row['word2']}/{row['word3']}"
                angles = f"{row['theta12']}/{row['theta13']}/{row['theta23']}"
                print(f"  {rank:<5} {row['manifold']:<10} {row['volume']:<8} "
                      f"{str(row['homology']):<8} {words:<30} "
                      f"{row['sigma']:<7} {row['fitness']:<10.6f} {angles}")
        except Exception as e:
            print(f"  No results yet ({e})")
    else:
        print("  No results yet.")

    print()
    print(f"  [Refreshing every {REFRESH}s — Ctrl+C to exit]")
    print("=" * 70)

print("PMNS Monitor starting... (Ctrl+C to stop)")
while True:
    show()
    time.sleep(REFRESH)