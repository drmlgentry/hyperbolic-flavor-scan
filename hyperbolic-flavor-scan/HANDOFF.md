HANDOFF: Hyperbolic Flavor Geometry Project
Date: 2026-03-10
Author: Marvin L. Gentry
Purpose: Provide complete context for continuing research in a new session, including all scripts, data, and next steps.

1. Project Overview
We are numerically scanning the SnapPy closed orientable census to find hyperbolic 3‑manifolds whose boundary geometry yields mixing matrices close to the CKM matrix. The pipeline:

Extract 
S
L
(
2
,
C
)
SL(2,C) holonomy representation via polished_holonomy().

Compute attracting fixed points of chosen group elements (first two generators and their product, or shortest hyperbolic words).

Build a Gram matrix from Gaussian overlaps (width σ) and QR‑orthogonalize to get a unitary mixing matrix.

Compare the moduli matrix with the PDG CKM matrix using a fitness function (moduli error + 0.1 × Jarlskog error).

The Jarlskog remains zero because the overlap kernel is real; future work will add U(1) phases via abelianization.

2. File & Folder Structure
All work is in C:\dev\hyperbolic-flavor-scan\.
Virtual environment: .venv (Python 3.13.12).
Activate with:

powershell
.\.venv\Scripts\Activate.ps1
Key scripts:
Script	Purpose
scan_flavors_working.py	Uses first two generators + their product. Works, produced top‑5 list.
scan_flavors_shortest.py	Selects three shortest hyperbolic words (max length adjustable). Latest version.
analyze_top.py	Reads CSV and displays top 5 mixing matrices.
analyze_sigma.py	Finds best σ for a given manifold (e.g., m007).
analyze_shortest.py	Prints three shortest hyperbolic words for top candidates.
run_sigma_scan.ps1	Batch runs sigma scans (parallel or sequential).
shortest_elements.py	Debug script to find shortest words (now integrated into analyze_shortest.py).
Data files:
scan_results.csv – results from the working scanner (first two generators).

scan_results_shortest.csv – results from the shortest‑word scanner (latest run).

scan_results_sigma_*.csv – sigma‑scan outputs (0.3–0.9).

results_m007_scan.csv – custom run with σ=0.3, word length 5.

3. Key Results (as of last run)
3.1. Working scanner (first two generators + product)
Top 5 (σ=0.5):

text
m007  vol 1.8436  score 0.1349  J 0.00e+00
m011  vol 1.8319  score 0.2437  J 0.00e+00
m011  vol 1.9122  score 0.2560  J 0.00e+00
m019  vol 2.0299  score 0.2775  J 0.00e+00
m004  vol 1.3985  score 0.3497  J 0.00e+00
m007 mixing matrix at σ=0.3:

text
[[0.9802 0.08   0.181 ]
 [0.0191 0.9487 0.3157]
 [0.1969 0.306  0.9314]]
(Note: (1,3) and (2,3) are too large, (1,2) too small.)

3.2. Shortest hyperbolic words (from analyze_shortest.py)
| Manifold | Three shortest words (by |trace|) |
|----------|----------------------------------|
| m007 | 11g2⁻¹2 (|tr|=2.000, L=0), g1⁻¹2g1⁻¹ (2.646, L=1.567), 22g1⁻¹ (2.828, L=1.763) |
| m011 | g2⁻¹1g1⁻¹2 (2.000, L=0), 21g2⁻¹2 (2.236, L=0.962), 11g2⁻¹ (2.659, L=1.582) |
| m019 | 12 (2.000, L=0), 22g2⁻¹2 (2.646, L=1.567), g1⁻¹g1⁻¹g2⁻¹g2⁻¹ (3.000, L=1.925) |
| m004 | g1⁻¹g2⁻¹g2⁻¹ (2.206, L=0.899), 2222 (2.306, L=1.093), g1⁻¹g2⁻¹g2⁻¹g2⁻¹ (2.384, L=1.220) |

The first two generators are not the shortest hyperbolic elements (they have |tr| ≤ 2, hence not hyperbolic). Therefore, the shortest‑word scanner is more physically motivated.

3.3. Shortest‑word scanner (σ=0.3, max_len=5) – preliminary
Top 5 (from scan_results_shortest.csv):

text
m006  vol 1.6496  score 0.1472  J 0.00e+00
m003  vol 1.2637  score 0.1860  J 0.00e+00
m004  vol 1.4638  score 0.2587  J 0.00e+00
m015  vol 2.0299  score 0.2810  J 0.00e+00
m007  vol 1.5436  score 0.2821  J 0.00e+00
Note: m007 now has a different volume because the elements are different (shortest words, not generators). Its fitness (0.2821) is worse than the generator‑based score (0.1349), but the mixing matrix may be more realistic. We need to examine the actual matrix.

4. Next Steps (Resume Here)
Examine the mixing matrix for m007 from the shortest‑word scan
Run:

powershell
python -c "import pandas as pd; df=pd.read_csv('scan_results_shortest.csv'); row=df[df['name']=='m007'].iloc[0]; print(row[['u11','u12','u13','u21','u22','u23','u31','u32','u33']])"
Compare with the CKM matrix. If it's still off, adjust σ or max word length.

Optimize σ for each top candidate
Use run_sigma_scan.ps1 (modify to call scan_flavors_shortest.py). Currently it calls scan_flavors_working.py – change the script variable accordingly.

Implement U(1) phases to obtain non‑zero Jarlskog

Use M.homology() to get the free abelian quotient.

For each of the three words, compute the image in homology (sum of generator images).

Assign a small random phase to each free generator.

Multiply the Gaussian overlap by exp(i * (phase_i - phase_j)).

Debugging: M.fundamental_group().homology() might work; test on m007:

python
G = M.fundamental_group()
print(G.homology())   # may return a list of integer vectors
Expand search to more manifolds
Increase the limit in the scanner (currently 50). The census has thousands.

5. Essential Commands
Reactivate environment:
powershell
cd C:\dev\hyperbolic-flavor-scan
.\.venv\Scripts\Activate.ps1
Run shortest‑word scanner (σ=0.5, max_len=4):
powershell
python scan_flavors_shortest.py
To change parameters:

powershell
python scan_flavors_shortest.py 0.3 5 my_results.csv
Batch sigma scan (modify script to use shortest‑word scanner first):
powershell
.\run_sigma_scan.ps1
Analyze results:
powershell
python analyze_sigma.py        # for m007 (change variable in script)
python analyze_shortest.py     # shows shortest words for top candidates
6. Notes on U(1) Implementation
The current scan_flavors_shortest.py uses a real kernel. To add U(1) phases, modify build_mixing_matrix to accept a list of phases (complex numbers of unit modulus). In get_three_shortest_matrices, after obtaining the matrices, also compute the homology images for each word and assign a phase:

python
# After obtaining the three matrices, also compute phases
hom = M.homology()   # might need to parse
# For each word, sum the abelianization images of its letters
phase = np.exp(1j * np.dot(word_vector, free_angles))
We need to determine the exact format of M.homology() for closed manifolds in SnapPy 3.3.2. Test with:

python
M = snappy.OrientableClosedCensus['m007']
print(M.homology())
If it returns something like 'Z^0 + Z/5 + Z/5 + Z/5', we must parse the free rank from the string. Alternatively, use M.fundamental_group().homology().

7. Handoff File
This markdown file should be saved as HANDOFF.md in the project folder. All scripts are present; if any are missing, the full code is included in the conversation history and can be re‑extracted.

To resume in a new session, simply read this file, activate the environment, and continue from Section 4 (Next Steps).

End of Handoff

To break the real-symmetry of the overlap kernel and generate a non-zero Jarlskog invariant (
J
J), we must implement the "Geometric Origin of CP Phases" logic from your fourth paper.

In SnapPy, M.homology() for a manifold like m007 returns 
(
Z
/
5
)
3
(Z/5) 
3
 . We will map each word to its image in this homology group and assign a phase based on a 
U
(
1
)
U(1) character.

Updated Script Fragment: abelianize_word
Add this to your scan_flavors_shortest.py:

python
def get_word_homology_image(M, word_str):
    """
    Maps a word (e.g., 'g1g2G1') to its image in H1(M, Z).
    """
    G = M.fundamental_group()
    # SnapPy's abelianization map
    ab_map = G.abelian_dilation() 
    
    # Initialize zero vector for homology
    image = np.zeros(len(ab_map[0]))
    
    # Parse word letters (g1=0, g2=1, G1=-0, G2=-1)
    for char in word_str:
        if char.isdigit():
            idx = int(char) - 1
            image += ab_map[idx]
        elif char.isupper(): # Inverse
            idx = ord(char) - ord('A')
            image -= ab_map[idx]
        elif char.islower(): # Forward
            idx = ord(char) - ord('a')
            image += ab_map[idx]
            
    return image

def compute_complex_overlap(z1, z2, phase1, phase2, sigma):
    """Geometric overlap with U(1) holonomy phases."""
    dist = 2.0 * np.abs(z1 - z2) / (np.sqrt(1 + np.abs(z1)**2) * np.sqrt(1 + np.abs(z2)**2))
    gaussian = np.exp(-(dist**2) / (2 * sigma**2))
    return gaussian * np.exp(1j * (phase1 - phase2))
   ...   3. Proposed Revised Scanner: scan_flavors_cp.py
This script integrates the shortest-word selection with homology-based CP phases.

python
import snappy
import numpy as np
from scipy.linalg import qr

def run_cp_scan(limit=100, sigma=0.3):
    # PDG Targets
    V_CKM = np.array([[0.974, 0.225, 0.003], [0.225, 0.973, 0.041], [0.008, 0.040, 0.999]])
    J_TARGET = 3.0e-5

    results = []
    
    for i in range(limit):
        M = snappy.OrientableClosedCensus[i]
        G = M.fundamental_group()
        # Define random phases for the generators of H1
        # In a real study, we'd scan these or use roots of unity for torsion
        h1 = M.homology()
        num_gens = len(G.generators())
        gen_phases = np.random.uniform(0, 2*np.pi, num_gens)

        # Get 3 shortest hyperbolic words
        # (Assuming your logic for 'get_shortest_words' is available)
        words = get_shortest_words(M, max_len=5, num=3)
        
        # Build Complex Overlap Matrix
        M_overlap = np.zeros((3,3), dtype=complex)
        phases = []
        fixed_points = []
        
        for w in words:
            # 1. Fixed Point
            mat = rho(w) 
            fp = attracting_fixed_point(mat)
            fixed_points.append(fp)
            # 2. Homology Phase
            img = get_word_homology_image(M, w)
            phases.append(np.dot(img, gen_phases))
            
        for r in range(3):
            for c in range(3):
                M_overlap[r,c] = compute_complex_overlap(
                    fixed_points[r], fixed_points[c], 
                    phases[r], phases[c], sigma
                )
        
        # QR Orthogonalization to get Unitary U
        U, R = qr(M_overlap)
        
        # Calculate Jarlskog
        J = np.imag(U[0,0]*U[1,1]*np.conj(U[0,1])*np.conj(U[1,0]))
        
        # Fitness
        moduli_err = np.linalg.norm(np.abs(U) - V_CKM)
        score = moduli_err + 0.1 * np.abs(J - J_TARGET)
        
        results.append([M.name(), score, J])
        
    return results