# ============================================
# MODULE: spectral_analysis.py
# ============================================
"""
Spectral analysis and φ-quantization tools.
"""

import mpmath as mp
from dataclasses import dataclass
from typing import List, Tuple, Optional

@dataclass
class SpectralData:
    """Comprehensive spectral analysis of a Coxeter word."""
    matrix: mp.matrix
    eigenvalues: List[complex]
    expanding_ev: Optional[mp.mpf]  # Real > 1
    translation_length: mp.mpf
    twist_angle: mp.mpf
    is_hyperbolic: bool
    is_twist_free: bool
    
    @property
    def q_phi(self) -> mp.mpf:
        """φ-quantization ratio: q = ℓ/log(φ)."""
        phi = (1 + mp.sqrt(5)) / 2
        return self.translation_length / mp.log(phi)

class SpectralAnalyzer:
    """High-precision spectral analysis for Coxeter words."""
    
    def __init__(self, precision: int = 200):
        self.precision = precision
        mp.mp.dps = precision
    
    def analyze(self, matrix: mp.matrix) -> SpectralData:
        """Perform comprehensive spectral analysis."""
        eigvals = mp.eig(matrix)[0]
        
        # Find expanding eigenvalue (real > 1)
        expanding_ev = None
        for ev in eigvals:
            if abs(ev.imag) < 1e-50 and ev.real > 1.0:
                expanding_ev = ev.real
                break
        
        # Compute translation length
        if expanding_ev:
            translation_length = mp.log(expanding_ev)
            is_hyperbolic = True
        else:
            translation_length = mp.mpf(0)
            is_hyperbolic = False
        
        # Find twist angle from complex eigenvalues on unit circle
        twist_angle = mp.mpf(0)
        for ev in eigvals:
            if 0.99 < abs(ev) < 1.01 and abs(ev.imag) > 1e-50:
                angle = mp.arg(ev)
                if angle > mp.pi:
                    angle -= 2 * mp.pi
                twist_angle = angle
                break
        
        is_twist_free = abs(twist_angle) < 1e-12
        
        return SpectralData(
            matrix=matrix,
            eigenvalues=eigvals,
            expanding_ev=expanding_ev,
            translation_length=translation_length,
            twist_angle=twist_angle,
            is_hyperbolic=is_hyperbolic,
            is_twist_free=is_twist_free
        )
    
    def compute_minimal_polynomial(self, eigenvalue: mp.mpc) -> Optional[List[mp.mpf]]:
        """Find minimal polynomial for eigenvalue over Q(√5)."""
        # Try to recognize as (a + b√5)/2 form
        u = eigenvalue + 1/eigenvalue if eigenvalue != 0 else None
        
        if u and abs(u.imag) < 1e-50:
            u_real = u.real
            sqrt5 = mp.sqrt(5)
            
            # Try to express u as (a + b√5)/2
            best_err = mp.inf
            best_coeffs = None
            
            for a in range(-20, 21):
                for b in range(-20, 21):
                    candidate = (a + b * sqrt5) / 2
                    err = abs(u_real - candidate)
                    
                    if err < best_err:
                        best_err = err
                        if err < 1e-40:
                            # Found representation, compute quartic for eigenvalue
                            # x^4 - a x^3 + ((a^2 - 5b^2)/4 + 2) x^2 - a x + 1
                            mid_num = (a**2 - 5 * b**2)
                            if mid_num % 4 == 0:
                                mid = mid_num // 4 + 2
                                best_coeffs = [mp.mpf(1), mp.mpf(-a), mp.mpf(mid), 
                                              mp.mpf(-a), mp.mpf(1)]
            
            if best_coeffs:
                # Verify polynomial
                poly_val = mp.polyval(best_coeffs, eigenvalue)
                if abs(poly_val) < 1e-40:
                    return best_coeffs
        
        return None