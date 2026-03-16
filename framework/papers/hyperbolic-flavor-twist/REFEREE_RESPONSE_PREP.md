# Paper 8 Referee Response Preparation
# gentry-hyperbolic-flavor-twist (NPB submission)
# Prepared: 2026-03-16

---

## MAJOR COMMENTS

### 1. Selberg zeta computation

**Referee concern:** Truncated product; λ1 values need justification.

**Key fact discovered post-submission:**
m003 = Meyerhoff manifold (vol=0.9814), NOT the Weeks manifold (vol=0.9427).
Inoue (CQG 18, 629, 2001) λ1 ≈ 27.8 applies to the Weeks manifold only.
No published λ1 exists for the Meyerhoff manifold or m006.
Our estimates λ1(m003)≈2.48, λ1(m006)≈2.82 are the FIRST for these manifolds.

**Response text:**
"We clarify that m003 (OrientableClosedCensus[1], vol=0.9814) is the Meyerhoff
manifold, distinct from the Weeks manifold (OrientableClosedCensus[0], vol=0.9427)
for which Inoue (2001) computed λ1 ≈ 27.8. No rigorous λ1 computation for the
Meyerhoff manifold or m006 currently exists in the literature; our estimates are
the first for these specific manifolds. We acknowledge the truncation is not
convergent in the strict sense (as stated in Appendix A) and note that the
Booker-Strömbergsson trace formula method [Lin-Lipnowski, JAMS 35, 2022] provides
a path to rigorous bounds, which we defer to future work."

**New refs:** Lin-Lipnowski (arXiv:1810.06346), Inoue (CQG 18, 2001),
             Booker-Strömbergsson (J. Reine Angew. Math. 607, 2007),
             Meyerhoff (London Math. Soc. 112, 1986)

---

### 2. mb/mc scheme dependence

Already handled in submitted text: m̄_b(m̄_b)/m̄_c(m̄_c) = 4.18/1.27 = 3.2913
(PDG 2024 MSbar running masses at own scales). No revision needed.

---

### 3. MZ/MW -- 6σ issue

Already stated explicitly in submitted Section 4. Add physical mechanism:

**Addition to Section 4:**
"The Hosotani mechanism on compact hyperbolic manifolds [Kaloper et al., PRL 85,
928, 2000; Demir-Shifman, PRD 65, 104002, 2002] provides a natural framework
in which gauge boson masses depend on holonomy data. Since H1(m006)=Z/5, the
manifold admits Z/5-valued Wilson lines; in gauge-Higgs unification on CHMs,
MZ/MW receives geometry-dependent corrections. A quantitative derivation is
left to future work."

**New refs:** Kaloper et al. (hep-ph/0002001), Demir-Shifman (hep-ph/0112090)

---

### 4. θ13_CKM prediction -- soften + update with length-7 result

**New post-submission result (length-7 census):**
- m006 spectral floor drops to φ_fold = 0.0053° at word length 6
  (near-identity geodesic family: aaabab and 25 conjugates, |λ|=3.730)
- Clean spectral gap [0.005°, 1.61°] at lengths ≤ 7
- θ13_CKM = 0.201° sits inside this gap
- Gap is likely arithmetic in origin (H1=Z/5 torsion)
- Original length-12 heuristic does not account for this gap structure

**Replacement text for Section 5:**
"Our length-7 census (completed after submission) reveals that the spectral floor
on m006 drops to φ_fold = 0.0053° at word length 6, arising from a near-identity
geodesic family (canonical representative aaabab, 26 words total, |λ|=3.730).
The target θ13_CKM = 0.201° falls in a clean spectral gap [0.005°, 1.61°] that
persists through length 7. This gap is likely arithmetic in origin, related to the
Z/5 torsion in H1. Whether a geodesic with φ_fold ≈ 0.201° exists at some longer
word length remains open. The original heuristic prediction of ℓ* ≈ 12 was based
on a multiplicative accumulation model that does not account for this gap structure
and should be regarded as superseded."

---

### 5. Table 2 caption -- slash notation

Add to caption: "For ratio rows (fold = r), W1/W2 denotes
|λ(W1)| > |λ(W2)|; the length column gives max(|W1|, |W2|)."

---

### 6. Experimental uncertainty table

Already present as Table 3 in submitted version. Confirm referee has seen it.

---

## MINOR COMMENTS

| Issue | Action |
|---|---|
| Abstract "encodes" | Already "appears to encode" in submitted version |
| eq:lam1 label bug | Check \label{} placement in Sec 6.1 equations |
| delta_CP word triple in Table 1 | Add footnote explaining signed sum |
| θ13_ν 34% error | Add: "listed for completeness; not claimed as coincidence" |
| Section 6.1 method description | Covered by Appendix A |
| Cheeger constant discussion | Add sentence citing Buser (1982) re large gap |

---

## NEW RESULTS AVAILABLE FOR RESPONSE

1. Length-7 census complete -- spectral gap [0.005°, 1.61°] on m006
2. m003 = Meyerhoff manifold (confirmed, not Weeks)
3. 5× torsion pattern: 5×1.687° = 8.435° ≈ θ13_ν = 8.54° (Paper 9)
4. Near-identity geodesic aaabab at length 6 on m006

---

## REFERENCES TO ADD IN REVISION

| Reference | Where |
|---|---|
| Kaloper, March-Russell, Starkman, Trodden (PRL 85, 928, 2000) | Sec 4 MZ/MW |
| Demir, Shifman (PRD 65, 104002, 2002) | Sec 4 MZ/MW |
| Lin, Lipnowski (JAMS 35, 2022; arXiv:1810.06346) | App A Selberg |
| Booker, Strömbergsson (J. Reine Angew. Math. 607, 2007) | App A Selberg |
| Inoue (CQG 18, 629, 2001) | App A -- clarify applies to Weeks |
| Meyerhoff (London Math. Soc. 112, 1986) | Sec 2 manifold ID |
