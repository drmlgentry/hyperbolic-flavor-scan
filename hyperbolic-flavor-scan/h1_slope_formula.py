import snappy, math

print("H1 SLOPE FORMULA TEST")
print()

def h1_info(name, p, q):
    try:
        M = snappy.Manifold(name)
        M.dehn_fill((p, q))
        sol = M.solution_type()
        h1 = M.homology()
        return str(h1), sol
    except Exception as e:
        return f"ERR:{e}", ""

print("m003 fillings:")
for p,q in [(-2,3),(-5,2),(-1,2),(-1,3),(-3,2),(2,3),(1,2),(1,3),
            (-2,1),(-3,1),(-1,4),(-4,1),(-1,5),(-5,1)]:
    if math.gcd(abs(p),q) != 1: continue
    h1, sol = h1_info("m003", p, q)
    z5 = " <-- Z/5" if h1 == "Z/5" else ""
    hfg = " <-- HFG" if (p,q) in [(-2,3),(-5,2)] else ""
    print(f"  ({p:>3},{q}): H1={h1:>20}  sol={sol[:20]}{z5}{hfg}")

print()
print("Dehn surgery formula for m003:")
print("  H1(m003) = Z + Z/5")
print("  Surgery on slope (p,q) kills the class p*mu + q*lam")
print("  where mu,lam are meridian and longitude")
print("  H1(m003(p,q)) is determined by the linking form")
