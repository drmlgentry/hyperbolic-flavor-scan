path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-torsion-gd.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# 1. Promote covering structure paragraph to subsection in Discussion
tex = tex.replace(
    r"\paragraph{Covering structure and the Z/5 tower.}",
    r"\subsection{Covering structure and $\mathbb{Z}/5$-preserving towers}"
)
tex = tex.replace(
    r"\paragraph{Covering structure and $\mathbb{Z}/5$-preserving towers.}",
    r"\subsection{Covering structure and $\mathbb{Z}/5$-preserving towers}"
)

# 2. Promote SnapPy naming paragraph to subsection and move signal
# First promote to subsection wherever it appears
tex = tex.replace(
    r"\paragraph{SnapPy naming note.}",
    r"\subsection*{SnapPy naming note}"
)

# 3. Unify Discussion subsections - make Floors vanish and Margulis into \subsection*
tex = tex.replace(r"\subsection*{Floors vanish}", r"\subsection*{Floors vanish}")
tex = tex.replace(r"\subsection*{Margulis asymptotics}", r"\subsection*{Margulis asymptotics}")
# These are already subsection* so just ensure covering is promoted too

# 4. Fix the duplicate orphan line if still present
lines = tex.splitlines(keepends=True)
clean = [l for l in lines 
         if not (r"protonmail.com}" in l and "ORCID" in l and l.strip().endswith("}}"))
         and not (l.strip() == r"\texttt{drmlgentry@protonmail.com}. ORCID: 0009-0006-4550-2663.}}")]
tex = "".join(clean)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print(f"Done. Lines: {len(tex.splitlines())}")
print(f"Covering subsection: {'Covering structure' in tex}")
print(f"SnapPy note: {'SnapPy naming note' in tex}")
