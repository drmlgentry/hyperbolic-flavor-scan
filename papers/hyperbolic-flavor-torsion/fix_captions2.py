path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Replace the broken second caption entirely
import re
# Find and replace the broken caption on the fig_coset_floors figure
old = re.search(r'\\caption\{Spectral floors per homology class.*?\}', tex, re.DOTALL)
if old:
    print(f"Found broken caption: {repr(old.group(0)[:80])}")
    tex = tex.replace(old.group(0),
        r"\caption{Spectral floors $\phi_{\rm floor}(k)$ per homology class "
        r"$k\in\mathbb{Z}/5$ for \mfld{m003} (left, Meyerhoff) and \mfld{m006} "
        r"(right), from the word-length-$\leq 7$ census. "
        r"The $637\times$ asymmetry on \mfld{m006} (class-4 floor: $0.005^\circ$; "
        r"next lowest: $1.61^\circ$) produces the spectral gap. "
        r"On \mfld{m003} all five floors lie within a factor of $2\times$.}")

# Remove duplicate figure blocks -- keep only the first occurrence of each
# Find all figure* blocks
blocks = list(re.finditer(r'\\begin\{figure\*\}.*?\\end\{figure\*\}', tex, re.DOTALL))
print(f"Found {len(blocks)} figure* blocks")
for b in blocks:
    print(f"  {b.group(0)[:60]}")

# If there are more than 2, remove duplicates
if len(blocks) > 2:
    # Keep only first two unique figures
    seen = set()
    for b in blocks:
        key = re.search(r'\\includegraphics.*?\{(.*?)\}', b.group(0))
        if key:
            fname = key.group(1)
            if fname in seen:
                tex = tex.replace(b.group(0), '', 1)
                print(f"  Removed duplicate: {fname}")
            else:
                seen.add(fname)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Done.")
