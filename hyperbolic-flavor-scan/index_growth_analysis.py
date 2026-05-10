import snappy
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re

M3 = snappy.OrientableClosedCensus[1]
v3 = float(M3.volume())

print(f"M3/PMNS base: vol={v3:.8f}")
print(f"Scanning degrees 2-8 for covering tower indices...\n")

def prime_factors(n):
    n = abs(n)
    if n <= 1:
        return set()
    factors = set()
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.add(d)
            n //= d
        d += 1
    if n > 1:
        factors.add(n)
    return factors

# Collect all covers by degree
tower = {1: [(1, 'm003', v3, 'Z/5')]}  # base manifold
n_census = len(snappy.OrientableClosedCensus)

for deg in range(2, 9):
    target = deg * v3
    covers = []
    for i in range(n_census):
        if i % 3000 == 0 and deg == 2:
            print(f"  Scanning census... {i}/{n_census}")
        M = snappy.OrientableClosedCensus[i]
        v = float(M.volume())
        if abs(v - target) < 1e-5:
            covers.append((i, M.name(), v, str(M.homology())))
    tower[deg] = covers
    print(f"Degree {deg}: {len(covers)} covers at vol={target:.6f}, "
          f"indices={[c[0] for c in covers]}")

print("\n=== Index Growth Analysis ===")
print(f"{'Degree':>8} {'Count':>7} {'Min idx':>9} {'Max idx':>9} "
      f"{'Mid idx':>9} {'Mean idx':>10}")

degrees = []
counts = []
min_indices = []
max_indices = []
mid_indices = []
mean_indices = []

for deg in sorted(tower.keys()):
    covers = tower[deg]
    if not covers:
        continue
    indices = [c[0] for c in covers]
    count = len(covers)
    mn = min(indices)
    mx = max(indices)
    mid = indices[len(indices)//2]
    mean = np.mean(indices)
    print(f"{deg:>8} {count:>7} {mn:>9} {mx:>9} {mid:>9} {mean:>10.1f}")
    degrees.append(deg)
    counts.append(count)
    min_indices.append(mn)
    max_indices.append(mx)
    mid_indices.append(mid)
    mean_indices.append(mean)

degrees = np.array(degrees)
counts = np.array(counts)
mean_idx = np.array(mean_indices)
min_idx = np.array(min_indices)

print("\n=== Power Law Fit: mean index ~ degree^alpha ===")
# Fit log(mean_index) ~ alpha * log(degree) + const
log_deg = np.log(degrees[1:])   # skip degree 1 (base)
log_mean = np.log(mean_idx[1:])
log_min  = np.log(min_idx[1:])

slope_mean, intercept_mean, r_mean, _, _ = stats.linregress(log_deg, log_mean)
slope_min,  intercept_min,  r_min,  _, _ = stats.linregress(log_deg, log_min)

print(f"Mean index fit:  index ~ deg^{slope_mean:.3f}  (R²={r_mean**2:.4f})")
print(f"Min  index fit:  index ~ deg^{slope_min:.3f}  (R²={r_min**2:.4f})")

print("\n=== Subgroup Growth Analysis ===")
print(f"{'Degree':>8} {'Count':>7} {'Ratio to prev':>15} {'Cumulative':>12}")
cumulative = 0
prev_count = 1
for i, deg in enumerate(degrees):
    cumulative += counts[i]
    ratio = counts[i] / prev_count if prev_count > 0 else float('inf')
    print(f"{deg:>8} {counts[i]:>7} {ratio:>15.3f} {cumulative:>12}")
    prev_count = counts[i]

# Fit subgroup count growth
log_counts = np.log(counts[1:])
slope_cnt, intercept_cnt, r_cnt, _, _ = stats.linregress(log_deg, log_counts)
print(f"\nSubgroup count fit: count ~ deg^{slope_cnt:.3f}  (R²={r_cnt**2:.4f})")

# Compare to known subgroup growth types
print(f"\n=== Comparison to Known Subgroup Growth Types ===")
print(f"Measured exponent (mean index): {slope_mean:.3f}")
print(f"Measured exponent (count):      {slope_cnt:.3f}")
print(f"  Polynomial growth (surface groups): exponent ~ 2-3")
print(f"  Abelian groups Z^n:                 exponent = n")
print(f"  Free groups (exponential):          would not fit power law")
print(f"  Hyperbolic 3-manifold groups:       typically 2-4")

# Prime spectrum evolution across degrees
print("\n=== Prime Spectrum Evolution ===")
print(f"{'Degree':>8}  {'Primes seen':>40}  {'New primes':>15}")
seen_primes = set()
for deg in sorted(tower.keys()):
    covers = tower[deg]
    deg_primes = set()
    for _, _, _, h1 in covers:
        nums = [int(x) for x in re.findall(r'\d+', h1) if int(x) > 1]
        for num in nums:
            deg_primes |= prime_factors(num)
    new = deg_primes - seen_primes
    seen_primes |= deg_primes
    print(f"{deg:>8}  {str(sorted(deg_primes)):>40}  {str(sorted(new)):>15}")

# Plot
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle("Covering Tower of M3/PMNS: Index Growth Analysis", fontsize=14)

# Panel 1: mean index vs degree (log-log)
ax = axes[0,0]
ax.loglog(degrees[1:], mean_idx[1:], 'bo-', label='Mean census index')
ax.loglog(degrees[1:], min_idx[1:], 'rs-', label='Min census index')
fit_y = np.exp(intercept_mean) * degrees[1:]**slope_mean
ax.loglog(degrees[1:], fit_y, 'b--', alpha=0.5,
          label=f'Power law fit α={slope_mean:.2f}')
ax.set_xlabel('Cover degree')
ax.set_ylabel('Census index')
ax.set_title('Index Growth (log-log)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel 2: cover count vs degree
ax = axes[0,1]
ax.plot(degrees, counts, 'go-', linewidth=2)
ax.set_xlabel('Cover degree')
ax.set_ylabel('Number of covers')
ax.set_title('Subgroup Count Growth')
ax.grid(True, alpha=0.3)
for d, c in zip(degrees, counts):
    ax.annotate(str(c), (d, c), textcoords="offset points",
                xytext=(0, 8), ha='center')

# Panel 3: cumulative prime spectrum
ax = axes[1,0]
all_primes_by_deg = []
seen = set()
for deg in sorted(tower.keys()):
    for _, _, _, h1 in tower[deg]:
        nums = [int(x) for x in re.findall(r'\d+', h1) if int(x) > 1]
        for num in nums:
            seen |= prime_factors(num)
    all_primes_by_deg.append((deg, sorted(seen)))

prime_list = [2, 3, 5, 7, 11, 13, 29]
prime_colors = plt.cm.tab10(np.linspace(0, 1, len(prime_list)))
for pi, (p, col) in enumerate(zip(prime_list, prime_colors)):
    first_deg = None
    for deg, primes in all_primes_by_deg:
        if p in primes:
            first_deg = deg
            break
    if first_deg:
        ax.barh(pi, 9 - first_deg, left=first_deg, color=col,
                alpha=0.7, label=f'p={p}')
        ax.text(first_deg + 0.1, pi, f'  first at deg {first_deg}',
                va='center', fontsize=8)
ax.set_yticks(range(len(prime_list)))
ax.set_yticklabels([f'p={p}' for p in prime_list])
ax.set_xlabel('Degree')
ax.set_title('Prime Activation Degrees')
ax.set_xlim(1, 9)
ax.grid(True, alpha=0.3)

# Panel 4: index spacing within each degree group
ax = axes[1,1]
for deg in sorted(tower.keys()):
    covers = tower[deg]
    if len(covers) < 2:
        continue
    indices = sorted([c[0] for c in covers])
    spacings = np.diff(indices)
    ax.scatter([deg]*len(spacings), spacings, alpha=0.7, s=60)
ax.set_xlabel('Cover degree')
ax.set_ylabel('Index spacing within group')
ax.set_title('Intra-group Index Spacings')
ax.grid(True, alpha=0.3)

plt.tight_layout()
outpath = r'C:\dev\hyperbolic-flavor-scan\covering_tower_analysis.png'
plt.savefig(outpath, dpi=150, bbox_inches='tight')
print(f"\nPlot saved to: {outpath}")
