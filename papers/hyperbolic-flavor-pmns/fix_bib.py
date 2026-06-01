path = r"C:\dev\framework\papers\hyperbolic-flavor-pmns\gentry-pmns-rip.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Add bibliographystyle before \bibliography
if "bibliographystyle" not in tex:
    tex = tex.replace(
        r"\bibliography{gentry-pmns}",
        r"\bibliographystyle{unsrt}" + "\n" + r"\bibliography{gentry-pmns}"
    )
    print("Added bibliographystyle")

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Done.")
