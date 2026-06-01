path = r"C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cpNotes.bib"
with open(path, encoding="utf-8") as f:
    content = f.read()

# Remove revtex @CONTROL entries
import re
content = re.sub(r'@CONTROL\{[^}]*\}\n?', '', content)
content = content.lstrip()

with open(path, "w", encoding="utf-8") as f:
    f.write(content)
print("Done. Remaining @CONTROL:", content.count("@CONTROL"))
print("Total entries:", content.count("@"))
