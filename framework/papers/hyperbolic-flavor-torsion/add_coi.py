path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

coi = (
    r"\section*{Conflict of Interest}" + "\n"
    r"The author declares no conflicts of interest." + "\n\n"
)

tex = tex.replace(r"\begin{thebibliography}", coi + r"\begin{thebibliography}", 1)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("COI statement added.")
