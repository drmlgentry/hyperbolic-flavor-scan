path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-torsion-gd.tex"
with open(path, encoding="utf-8") as f:
    lines = f.readlines()

# Show line 24
print(f"Line 24: {repr(lines[23])}")

# Fix: rewrite line 24 cleanly
for i, line in enumerate(lines):
    if "ORCID: 0009-0006-4550-2663" in line and r"\affil" in line:
        lines[i] = r"\affil[1]{Independent Researcher, Seattle, WA, USA. \texttt{drmlgentry@protonmail.com}. ORCID: 0009-0006-4550-2663}" + "\n"
        print(f"Fixed line {i+1}: {repr(lines[i])}")

with open(path, "w", encoding="utf-8") as f:
    f.writelines(lines)
print("Done.")
