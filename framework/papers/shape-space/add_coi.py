path = r"C:\dev\framework\papers\shape-space\gentry-shape-space.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

if "Conflict of Interest" not in tex:
    tex = tex.replace(
        r"\begin{thebibliography}",
        r"\section*{Data Availability}" + "\n"
        r"This manuscript contains no research data." + "\n\n"
        r"\section*{Conflict of Interest}" + "\n"
        r"The author declares no conflicts of interest." + "\n\n"
        r"\begin{thebibliography}"
    )
    with open(path, "w", encoding="utf-8") as f:
        f.write(tex)
    print("Added.")
else:
    print("Already present.")
