import re, os

# The companion sentence to insert
SENTENCE = (
    "\nThe geometric construction developed here provides the theoretical\n"
    "foundation for the companion results in~\\cite{Gentry:CKM,Gentry:PMNS},\n"
    "where the CKM and PMNS mixing matrices are derived explicitly from the\n"
    "$K$- and $N$-factors of the Iwasawa decomposition of compact hyperbolic\n"
    "3-manifold holonomy~\\cite{Gentry:CKM,Gentry:PMNS}.\n"
)

# The two bib entries to append
BIB_ENTRIES = """
@article{Gentry:CKM,
  author  = {Gentry, Marvin L.},
  title   = {{CKM} quark mixing from geodesic axes of hyperbolic 3-manifold holonomy},
  journal = {Phys. Rev. D},
  year    = {2026},
  note    = {submitted March 2026, temporary {ID} es2026mar11\\_966}
}

@article{Gentry:PMNS,
  author  = {Gentry, Marvin L.},
  title   = {Lepton mixing from Borel structure of hyperbolic holonomy},
  journal = {Phys. Rev. D},
  year    = {2026},
  note    = {submitted March 2026, temporary {ID} es2026mar13\\_942}
}
"""

papers = [
    {
        'tex': r'C:\dev\framework\papers\holonomy-cp\gentry-holonomy-cp.tex',
        'bib': r'C:\dev\framework\papers\holonomy-cp\gentry-holonomy-cp.bib',
        'name': 'holonomy-cp (JGP12746)',
    },
    {
        'tex': r'C:\dev\framework\papers\flavor-mixing\gentry-flavor-mixing.tex',
        'bib': r'C:\dev\framework\papers\flavor-mixing\gentry-flavor-mixing.bib',
        'name': 'flavor-mixing (JGP12747)',
    },
    {
        'tex': r'C:\dev\framework\papers\shape-space\gentry-shape-space.tex',
        'bib': r'C:\dev\framework\papers\shape-space\gentry-shape-space.bib',
        'name': 'shape-space (JGP12753)',
    },
]

for p in papers:
    print('Processing: ' + p['name'])

    # --- Update .tex ---
    with open(p['tex'], 'r', encoding='utf-8-sig') as f:
        tex = f.read()

    # Check if already updated
    if 'Gentry:CKM' in tex:
        print('  TEX: already contains Gentry:CKM -- skipping tex update')
    else:
        # Find the comment separator that ends the introduction
        # Pattern: the %--- line right before \section after intro
        # Insert sentence just before that line
        pattern = re.compile(
            r'(\n%-+\n)(\\section(?!\{Intro))',
            re.DOTALL
        )
        m = pattern.search(tex)
        if m:
            tex = tex[:m.start()] + SENTENCE + tex[m.start():]
            with open(p['tex'], 'w', encoding='utf-8') as f:
                f.write(tex)
            print('  TEX: companion sentence inserted before Section 2')
        else:
            # Fallback: insert before first \section that is not Introduction
            pattern2 = re.compile(r'(\\section\{(?!Intro)[^}]+\})')
            m2 = pattern2.search(tex)
            if m2:
                tex = tex[:m2.start()] + SENTENCE + '\n' + tex[m2.start():]
                with open(p['tex'], 'w', encoding='utf-8') as f:
                    f.write(tex)
                print('  TEX: companion sentence inserted (fallback method)')
            else:
                print('  TEX: WARNING -- could not find insertion point')

    # --- Update .bib ---
    if not os.path.exists(p['bib']):
        print('  BIB: file not found at ' + p['bib'])
        continue

    with open(p['bib'], 'r', encoding='utf-8-sig') as f:
        bib = f.read()

    if 'Gentry:CKM' in bib:
        print('  BIB: already contains Gentry:CKM -- skipping bib update')
    else:
        with open(p['bib'], 'a', encoding='utf-8') as f:
            f.write(BIB_ENTRIES)
        print('  BIB: Gentry:CKM and Gentry:PMNS entries appended')

    print()

print('All done. Now recompile each paper.')
