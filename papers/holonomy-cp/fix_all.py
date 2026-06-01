path = r"C:\dev\framework\papers\holonomy-cp\gentry-holonomy-cp.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Fix stale IDs
tex = tex.replace("es2026mar11_966", "DQ14014")
tex = tex.replace("es2026mar13_942", "DQ14050")

# Check manifold
if "Weeks manifold" in tex:
    print("WARNING: Weeks manifold still present - need to paste revised file")
else:
    print("Meyerhoff OK")

# Check COI
if r"\section*{Conflict of Interest}" not in tex:
    tex = tex.replace(
        r"\begin{thebibliography}",
        r"\section*{Data Availability}" + "\n"
        r"This manuscript contains no research data." + "\n\n"
        r"\section*{Conflict of Interest}" + "\n"
        r"The author declares no conflicts of interest." + "\n\n"
        r"\begin{thebibliography}"
    )
    print("COI added.")

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
