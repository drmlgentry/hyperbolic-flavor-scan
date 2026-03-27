path = r"C:\dev\hyperbolic-flavor-scan\analysis\selberg_eigenvalue.py"
with open(path, encoding="utf-8") as f:
    s = f.read()
s = s.replace("MAX_LENGTH = 6.0", "MAX_LENGTH = 10.0")
s = s.replace("MAX_PRIMS  = 2000", "MAX_PRIMS  = 5000")
s = s.replace("MN_MAX     = 4",    "MN_MAX     = 5")
with open(path, "w", encoding="utf-8") as f:
    f.write(s)
print("Done.")
print("MAX_LENGTH:", "10.0" in open(path).read())
