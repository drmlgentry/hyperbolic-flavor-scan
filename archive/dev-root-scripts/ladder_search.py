# ============================================
# MODULE: ladder_search.py
# ============================================
"""
Structured search for twist-free ladder words in [3,5,3].
Family: W(m,n) = τ^m ρ τ^n ρ^{-1}
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Optional
import mpmath as mp

@dataclass
class LadderCandidate:
    """Candidate ladder word with spectral properties."""
    m: int
    n: int
    translation_length: mp.mpf  # ℓ(W)
    twist_angle: mp.mpf         # θ
    q_phi: mp.mpf              # q = ℓ/log(φ)
    q_nearest_int: int         # Nearest integer to q
    q_residual: mp.mpf         # |q - nearest_int|
    word_digits: str          # Digit representation
    spectral_data: dict       # Full spectral analysis
    
    def __str__(self):
        return (f"W({self.m},{self.n}): ℓ={float(self.translation_length):.6f}, "
                f"θ={float(self.twist_angle):.3e}, q={float(self.q_phi):.8f}, "
                f"|q-{self.q_nearest_int}|={float(self.q_residual):.3e}")

class LadderSearch:
    """Search for twist-free ladder words in Coxeter [3,5,3]."""
    
    def __init__(self, coxeter):
        self.coxeter = coxeter
        self.phi = (1 + mp.sqrt(5)) / 2
        self.log_phi = mp.log(self.phi)
    
    def search_range(self, m_range: Tuple[int, int], n_range: Tuple[int, int],
                    twist_tol: float = 1e-10, ell_min: float = 1e-6,
                    ell_max: float = 20.0) -> List[LadderCandidate]:
        """Search W(m,n) in given ranges."""
        m_min, m_max = m_range
        n_min, n_max = n_range
        
        candidates = []
        
        for m in range(m_min, m_max + 1):
            for n in range(n_min, n_max + 1):
                if m == 0 and n == 0:
                    continue
                
                W = self.coxeter.ladder_word_matrix(m, n)
                ell = self.coxeter.get_translation_length(W)
                theta = self.coxeter.get_twist_angle(W)
                
                # Apply filters
                if ell < ell_min or ell > ell_max:
                    continue
                if abs(theta) > twist_tol:
                    continue
                
                # Compute φ-quantization
                q = ell / self.log_phi
                q_int = int(mp.nint(q))
                q_resid = abs(q - q_int)
                
                # Build word digit representation
                word_digits = self._word_to_digits(m, n)
                
                # Full spectral analysis
                spectral_data = {
                    'ell': ell,
                    'theta': theta,
                    'eigenvalues': mp.eig(W)[0],
                    'trace': mp.trace(W)
                }
                
                candidate = LadderCandidate(
                    m=m, n=n,
                    translation_length=ell,
                    twist_angle=theta,
                    q_phi=q,
                    q_nearest_int=q_int,
                    q_residual=q_resid,
                    word_digits=word_digits,
                    spectral_data=spectral_data
                )
                
                candidates.append(candidate)
        
        # Sort by translation length, then twist
        candidates.sort(key=lambda c: (float(c.translation_length), 
                                      abs(float(c.twist_angle))))
        return candidates
    
    def _word_to_digits(self, m: int, n: int) -> str:
        """Convert W(m,n) to digit string in {0,1,2,3}."""
        # τ = s3 s0 = "30"
        # ρ = s1 s2 = "12"
        # ρ^{-1} = s2 s1 = "21"
        
        def tau_power(k: int) -> str:
            if k > 0:
                return "30" * k
            elif k < 0:
                return "03" * (-k)
            return ""
        
        tau_m = tau_power(m)
        tau_n = tau_power(n)
        
        return tau_m + "12" + tau_n + "21"
    
    def find_best_quantized(self, candidates: List[LadderCandidate],
                           max_residual: float = 0.01) -> List[LadderCandidate]:
        """Filter candidates with good φ-quantization."""
        return [c for c in candidates if float(c.q_residual) < max_residual]