path = r"C:\dev\framework\papers\hyperbolic-flavor-pmns\gentry-pmns-epjc.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Update CKM companion paper reference with EPJC accession number
tex = tex.replace(
    "submitted to Eur. Phys. J. C (2026)",
    "submitted to Eur. Phys. J. C (2026) (EPJC-26-03-349)"
)

# Also catch any remaining PRD reference to CKM
tex = tex.replace(
    "submitted March 2026 (DQ14014)",
    "submitted to Eur. Phys. J. C, EPJC-26-03-349 (2026)"
)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)

print("Done.")
print("EPJC-26-03-349 occurrences:", tex.count("EPJC-26-03-349"))
