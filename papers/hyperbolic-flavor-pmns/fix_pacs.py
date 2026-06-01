path = r"C:\dev\framework\papers\hyperbolic-flavor-pmns\gentry-pmns-rip.tex"
with open(path, encoding="utf-8") as f:
    lines = f.readlines()
new_lines = [l for l in lines if not l.strip().startswith(r'\pacs{')
             and not l.strip().startswith(r'\pacs ')]
with open(path, "w", encoding="utf-8") as f:
    f.writelines(new_lines)
print(f"Removed pacs. Lines: {len(new_lines)}")
