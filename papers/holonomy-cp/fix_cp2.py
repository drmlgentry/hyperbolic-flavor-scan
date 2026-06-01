path = r"C:\dev\framework\papers\holonomy-cp\gentry-holonomy-cp-epjc.tex"
with open(path, encoding="utf-8") as f:
    lines = f.readlines()

# Remove lines 13 and 14 (0-indexed: 12 and 13)
# These are the duplicate \texttt{email} and ORCID lines
new_lines = []
for i, line in enumerate(lines):
    if r"\texttt{drmlgentry@protonmail.com}" in line:
        print(f"Removed line {i+1}: {line.rstrip()}")
        continue
    if i > 8 and "ORCID: 0009-0006-4550-2663}" in line and "\\institute" not in lines[i-3]:
        print(f"Removed line {i+1}: {line.rstrip()}")
        continue
    new_lines.append(line)

with open(path, "w", encoding="utf-8") as f:
    f.writelines(new_lines)
print("Done.")
