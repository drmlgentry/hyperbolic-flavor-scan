path = r"C:\dev\framework\papers\flavor-mixing\gentry-flavor-mixing.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()
tex = tex.replace(
    r"\section*{Conflict of Interest}",
    r"\section*{Data Availability}" + "\n"
    r"This manuscript contains no research data. No data availability statement is required." + "\n\n"
    r"\section*{Conflict of Interest}"
)
with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Done.")
