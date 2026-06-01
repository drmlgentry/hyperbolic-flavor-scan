path = r"C:\dev\framework\papers\hyperbolic-flavor-pmns\gentry-pmns-epjc.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()
tex = tex.replace(
    r"\bibliography{gentry-pmns}",
    r"\bibliographystyle{spbasic}" + "\n" + r"\bibliography{gentry-pmns}"
)
with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Done.", r"\bibliographystyle{spbasic}" in tex)
