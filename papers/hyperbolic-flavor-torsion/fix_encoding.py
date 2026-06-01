path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    content = f.read()

# Remove BOM
content = content.lstrip('\ufeff')

# Replace em-dashes with LaTeX ---
content = content.replace('\u2014', '---')

with open(path, "w", encoding="utf-8") as f:
    f.write(content)
print("Fixed BOM and em-dashes.")
