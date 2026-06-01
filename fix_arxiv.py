import re, os

replacements = {
    "arXiv:????.?????": "",
}

epjc_ids = {
    "Gentry:CKM":   "EPJC-26-03-349, submitted (2026)",
    "Gentry:PMNS":  "EPJC-26-03-351, submitted (2026)",
    "Gentry:CP":    "EPJC-26-03-398, submitted (2026)",
    "Gentry:Twist": "EPJC-26-03-407, submitted (2026)",
}

for fpath in [
    r"C:\dev\framework\papers\hyperbolic-flavor-ckm\gentry-ckm-rip.tex",
    r"C:\dev\framework\papers\hyperbolic-flavor-pmns\gentry-pmns-rip.tex",
]:
    if not os.path.exists(fpath): continue
    with open(fpath, encoding="utf-8") as f:
        tex = f.read()
    for key, note in epjc_ids.items():
        tex = tex.replace(
            f"arXiv:????.?????",
            note
        )
    # Also fix any remaining ????
    tex = re.sub(r"arXiv:\?\?\?\?\.\?\?\?\?\?", "", tex)
    with open(fpath, "w", encoding="utf-8") as f:
        f.write(tex)
    print(f"Fixed: {fpath}")
