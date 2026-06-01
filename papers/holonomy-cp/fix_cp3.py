path = r"C:\dev\framework\papers\holonomy-cp\gentry-holonomy-cp-epjc.tex"
with open(path, encoding="utf-8") as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if r"\email{drmlgentry@protonmail.com} \\" in line:
        # Replace the next blank line with ORCID and closing brace
        if i+1 < len(lines) and lines[i+1].strip() == "":
            lines[i] = line.rstrip("\n") + "\n"
            lines[i+1] = "ORCID: 0009-0006-4550-2663}\n"
            print(f"Fixed at line {i+2}")
        break

with open(path, "w", encoding="utf-8") as f:
    f.writelines(lines)
print("Done.")
