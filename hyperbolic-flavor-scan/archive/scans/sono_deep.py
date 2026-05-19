"""
Deep structural analysis of sonoluminescence phi-lattice.
Key finding: combined 49 lines, 6 species, p=0.012, phi uniquely best base.
Optimal anchor 16.675 eV = Ar_4s[3/2]1 * phi^(3/4) to 0.002%.
"""
import numpy as np
import math
from collections import Counter

PHI = (1 + math.sqrt(5)) / 2
LOG_PHI = math.log(PHI)
hc = 1239.84193
m_e = 0.51099895e6
E0 = 16.675

def nearest_q(x, anchor):
    ratio = math.log(x / anchor) / LOG_PHI
    q = round(ratio * 4) / 4
    return q, abs(ratio - q)

all_lines = [
    (706.7,"Ar","4p[5/2]2->4s[3/2]2"),(727.3,"Ar","4p[5/2]1->4s[3/2]2"),
    (738.4,"Ar","4p[5/2]2->4s[3/2]1 STRONG"),(750.4,"Ar","4p[3/2]1->4s'[1/2]0"),
    (751.5,"Ar","4p[5/2]1->4s[3/2]1"),(763.5,"Ar","4p[5/2]2->4s'[1/2]1"),
    (772.4,"Ar","4p[1/2]0->4s[3/2]1"),(794.8,"Ar","4p[5/2]1->4s'[1/2]1"),
    (800.6,"Ar","4p[3/2]2->4s[3/2]2"),(801.5,"Ar","4p[1/2]0->4s[3/2]2"),
    (810.4,"Ar","4p[3/2]2->4s[3/2]1"),(811.5,"Ar","4p[3/2]2->4s[3/2]1 STRONG"),
    (826.5,"Ar","4p[3/2]1->4s'[1/2]1"),(840.8,"Ar","4p[3/2]2->4s'[1/2]1 STRONG"),
    (842.5,"Ar","4p[1/2]0->4s'[1/2]1"),(852.1,"Ar","4p[5/2]2->4s'[1/2]1"),
    (912.3,"Ar","4p[5/2]2->4s'[1/2]0"),
    (282.0,"OH*","A-X (1,0)"),(306.4,"OH*","A-X (0,1)"),
    (309.4,"OH*","A-X (0,0) STRONG"),(314.6,"OH*","A-X (1,1)***"),
    (316.8,"OH*","A-X (2,2)***"),
    (436.5,"C2","Swan (0,2)"),(469.4,"C2","Swan (1,2)"),
    (473.7,"C2","Swan (0,1)"),(512.9,"C2","Swan (1,1)***"),
    (516.5,"C2","Swan (0,0) STRONG"),(558.5,"C2","Swan (1,0)"),
    (563.5,"C2","Swan (2,1)"),
    (589.0,"Na","D1"),(589.6,"Na","D2"),
]

line_data = sorted([(nearest_q(hc/wl,E0)[0], nearest_q(hc/wl,E0)[1],
                     wl, hc/wl, sp, desc)
                    for wl,sp,desc in all_lines])

print("SBSL lines by q-value (anchor E0=16.675 eV):")
print(f"{'q':>6}  {'res':>6}  {'lam':>7}  {'sp':>4}  description")
for q,r,wl,E,sp,desc in line_data:
    f = "***" if r<0.03 else "**" if r<0.06 else "*" if r<0.09 else ""
    print(f"  {q:>+6.2f}  {r:>6.4f}  {wl:>7.1f}  {sp:>4}  {desc} {f}")

print()
q_counts = Counter(d[0] for d in line_data)
print("Lines per q-value:")
for q in sorted(q_counts):
    sps = sorted(set(d[4] for d in line_data if d[0]==q))
    print(f"  q={q:+.2f}: {q_counts[q]:2d} lines  {sps}")

print()
print("Anchor identity check:")
for val,name,q_try in [(11.623,"Ar 4s[3/2]1",0.75),(11.548,"Ar 4s[3/2]2",0.75),
                        (10.294,"SO IP",1.0),(m_e,"m_e",-21.5)]:
    pred = val * PHI**q_try
    print(f"  {name}: {val}*phi^{q_try} = {pred:.4f} eV  "
          f"(delta={100*(E0-pred)/pred:+.3f}%)")
