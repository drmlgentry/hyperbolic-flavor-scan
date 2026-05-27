#!/bin/bash
# HFG Session Capture — WSL/bash version
# Run from: /mnt/c/dev/

OUTFILE="/mnt/c/dev/HFG_SESSION_LOG_$(date +%Y-%m-%d).txt"
SEP="======================================================================"

section() {
    echo "" >> "$OUTFILE"
    echo "$SEP" >> "$OUTFILE"
    echo "  $1" >> "$OUTFILE"
    echo "$SEP" >> "$OUTFILE"
}

echo "HFG COMPLETE SESSION LOG" > "$OUTFILE"
echo "Generated: $(date)" >> "$OUTFILE"
echo "Researcher: Marvin L. Gentry | drmlgentry@protonmail.com" >> "$OUTFILE"
echo "ORCID: 0009-0006-4550-2663" >> "$OUTFILE"

# Git logs
for repo in hyperbolic-flavor-geometry hyperbolic-flavor-scan framework; do
    path="/mnt/c/dev/$repo"
    if [ -d "$path" ]; then
        section "GIT LOG: $repo (last 30 commits)"
        git -C "$path" log --oneline -30 2>&1 >> "$OUTFILE"
        
        section "GIT DIFF STAT: $repo"
        git -C "$path" diff --stat HEAD 2>&1 >> "$OUTFILE"
        
        section "GIT TRACKED FILES: $repo"
        git -C "$path" ls-files 2>&1 >> "$OUTFILE"
    fi
done

# Key canonical values
section "CANONICAL HFG VALUES"
cat >> "$OUTFILE" << 'VALUES'
M_PMNS = m003(-2,3) = OrientableClosedCensus[1]
  vol     = 0.981368828892232
  H_1     = Z/5 (created by surgery)
  CS      = 1/4 (exact)
  ITF     = degree-4, disc=-283, sig=(2,1)  [CFKR 2001]
  Fitness = 0.005087 (PMNS)
  delta   = 195.91 degrees
  sigma_opt = (3/2)*log(sqrt(13/5))  [exact]
  |WRT_r|^2 = 13/r  for r ≡ 1,3 mod 4  [proved May 25 2026]

M_CKM = m006(-5,2) = OrientableClosedCensus[43]
  vol     = 2.028853091474920
  H_1     = Z/5 (inherited from m006 cusped)
  CS      = ~4/15 (provisional)
  ITF     = degree-10, disc=-271488204251=-11*239*103266719
            poly: x^10-7x^8-4x^7+17x^6+14x^5-18x^4-14x^3+8x^2+3x-1
            sig=(8,1), Galois group S_10  [confirmed SnapPy May 25 2026]
  Fitness = 0.009078 (CKM, updated triple)

CUSPED PARENTS:
  m003: ITF = Q(sqrt(-3)), disc=-3, H_1=Z  (no torsion)
  m006: ITF = x^3-2x^2-1, disc=-59, sig=(1,1), H_1=Z/5+Z

CORRECTIONS FROM PRIOR PAPERS:
  Q(sqrt(17)) as ITF(M_CKM): WRONG (real quadratic, 0 complex places)
  phi as sigma_opt: APPROXIMATE (0.69% error, exact = sqrt(13/5))

JOURNAL STATUS (May 26 2026):
  PTEP T06112 / 2605-035: Unified HFG — under review
  PTEP T06113 / 2605-036: CP Letter — under review
  PLB-D-26-01448: PMNS (Kitano editor) — submitted
  PLB (TOPOL transfer): Twist spectrum — completing submission
  PLB (JGP13076 transfer): CKM — completing submission
  PLB-D-26-01341: CP phase — with editor
  PLB-D-26-01006: Charge conjugation — under review
  Universe universe-4362980: Twist spectrum — transferred
  Exp. Math. 266619148: Farey tower — under review
  SSRN 6815721: X_0(11) bridge — posted May 23
  SSRN 6583553: PMNS — updated May 27
VALUES

# Handoff files
section "HANDOFF DOCUMENTS"
for f in /mnt/c/dev/framework/HANDOFF*.md /mnt/c/dev/hyperbolic-flavor-scan/HANDOFF*.md; do
    if [ -f "$f" ]; then
        echo "" >> "$OUTFILE"
        echo "--- $f ---" >> "$OUTFILE"
        cat "$f" >> "$OUTFILE"
    fi
done

# Reproduce scripts
section "REPRODUCE SCRIPTS (canonical)"
for f in /mnt/c/dev/hyperbolic-flavor-scan/reproduce/*.py; do
    if [ -f "$f" ]; then
        echo "" >> "$OUTFILE"
        echo "--- $(basename $f) ---" >> "$OUTFILE"
        cat "$f" >> "$OUTFILE"
    fi
done

# Transcripts index
section "TRANSCRIPT FILES"
ls -la /mnt/transcripts/ 2>/dev/null >> "$OUTFILE" || echo "transcripts not accessible" >> "$OUTFILE"
cat /mnt/transcripts/journal.txt 2>/dev/null >> "$OUTFILE" || true

echo "Done. Log written to: $OUTFILE"
wc -l "$OUTFILE"
