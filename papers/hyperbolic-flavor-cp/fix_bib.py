path = r"C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp-epjc.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Check what bib files exist and update accordingly
import os
bib_files = [f.replace(".bib","") for f in os.listdir(
    r"C:\dev\framework\papers\hyperbolic-flavor-cp") if f.endswith(".bib")]
print("Bib files found:", bib_files)

# Update bibliography command to use correct file
for bib in bib_files:
    if "cp" in bib.lower() or "flavor" in bib.lower():
        tex = tex.replace(
            r"\bibliography{gentry-cp}",
            r"\bibliography{" + bib + "}"
        )
        print(f"Updated to: \\bibliography{{{bib}}}")
        break

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Done.")
