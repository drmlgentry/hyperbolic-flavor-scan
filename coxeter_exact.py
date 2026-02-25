# ============================================
# CORE MODULE: coxeter_exact.py
# ============================================
"""
Exact algebraic implementation of hyperbolic Coxeter [3,5,3] group.
Field: Q(√5) with high-precision arithmetic.
"""

import numpy as np
import mpmath as mp
from dataclasses import dataclass
from typing import Tuple, List, Dict, Optional

# Global precision setting
DEFAULT_PRECISION = 200
mp.mp.dps = DEFAULT_PRECISION

@dataclass
class CoxeterConfig:
    """Configuration for the [3,5,3] Coxeter system."""
    phi_lock: bool = True
    precision: int = DEFAULT_PRECISION
    use_exact_gram: bool = True

class Coxeter353:
    """Exact implementation of hyperbolic Coxeter [3,5,3] group."""
    
    def __init__(self, config: Optional[CoxeterConfig] = None):
        self.config = config or CoxeterConfig()
        mp.mp.dps = self.config.precision
        
        self.phi = (1 + mp.sqrt(5)) / 2
        self.log_phi = mp.log(self.phi)
        
        # Core matrices
        self.G = self._build_gram_matrix()
        self.J = self._build_minkowski_metric()
        self.S = self._build_reflections()
        
        # Composite generators
        self.tau = self.S[3] * self.S[0]  # τ = s3 s0
        self.rho = self.S[1] * self.S[2]  # ρ = s1 s2
        self.sigma = self.S[2] * self.S[3]  # σ = s2 s3
        
    def _build_gram_matrix(self) -> mp.matrix:
        """Build exact Gram matrix for [3,5,3] orthoscheme."""
        φ = self.phi
        
        # Angles: cos(π/3) = 1/2, cos(π/5) = φ/2
        g01 = -mp.mpf(1) / 2        # s0-s1: π/3
        g12 = -φ / 2                # s1-s2: π/5  
        g23 = -mp.mpf(1) / 2        # s2-s3: π/3
        g02 = mp.mpf(0)             # non-adjacent
        g13 = mp.mpf(0)             # non-adjacent
        
        # Ultraparallel: cosh(d*) = 1 + φ/2
        cosh_dstar = 1 + φ / 2
        g03 = -cosh_dstar           # s0-s3: ultraparallel
        
        return mp.matrix([
            [1,     g01, g02, g03],
            [g01,   1,   g12, g13],
            [g02,   g12, 1,   g23],
            [g03,   g13, g23, 1]
        ])
    
    def _build_minkowski_metric(self) -> mp.matrix:
        """Build Minkowski metric diag(-1,1,1,1)."""
        return mp.matrix([
            [-1, 0, 0, 0],
            [0,  1, 0, 0],
            [0,  0, 1, 0],
            [0,  0, 0, 1]
        ])
    
    def _build_reflections(self) -> Tuple[mp.matrix, ...]:
        """Build reflection generators s0..s3 exactly."""
        # Compute Lorentz normals via Cholesky-like factorization
        N = mp.zeros(4, 4)
        
        # n0 = (1, 0, 0, 0)
        N[0, 0] = mp.mpf(1)
        
        # n1: orthogonal to n0, satisfy n0·n1 = -g01
        N[1, 0] = -self.G[0, 1]
        N[1, 1] = mp.sqrt(1 - N[1, 0]**2)
        
        # n2: orthogonal to n0,n1
        N[2, 0] = mp.mpf(0)
        N[2, 1] = -self.G[1, 2] / N[1, 1]
        N[2, 2] = mp.sqrt(1 - N[2, 1]**2)
        
        # n3: orthogonal to all previous
        N[3, 0] = -self.G[0, 3]
        N[3, 1] = mp.mpf(0)
        N[3, 2] = -self.G[2, 3] / N[2, 2]
        N[3, 3] = mp.sqrt(1 - N[3, 0]**2 - N[3, 2]**2)
        
        # Build reflection matrices: R_i = I - 2 * (J n_i n_i^T) / <n_i, n_i>
        S = []
        for i in range(4):
            n = N[i,:].transpose()
            denom = (n.T * self.J * n)[0]
            R = mp.eye(4) - 2 * (self.J * n * n.T) / denom
            S.append(R)
        
        return tuple(S)
    
    def word_to_matrix(self, word: List[int]) -> mp.matrix:
        """Convert word in generators [0..3] to matrix product."""
        M = mp.eye(4)
        for d in word:
            M = M * self.S[int(d)]
        return M
    
    def ladder_word_matrix(self, m: int, n: int) -> mp.matrix:
        """Compute W(m,n) = τ^m ρ τ^n ρ^{-1}."""
        tau_m = mp.matrix_pow(self.tau, m)
        tau_n = mp.matrix_pow(self.tau, n)
        rho = self.rho
        rho_inv = mp.inverse(rho)
        
        return tau_m * rho * tau_n * rho_inv
    
    def get_translation_length(self, W: mp.matrix) -> mp.mpf:
        """Compute translation length ℓ(W) = log(spectral radius)."""
        eigvals = mp.eig(W)[0]
        max_abs = max(abs(ev) for ev in eigvals)
        return mp.log(max_abs) if max_abs > 1 else mp.mpf(0)
    
    def get_twist_angle(self, W: mp.matrix) -> mp.mpf:
        """Compute twist angle θ from complex eigenvalues on unit circle."""
        eigvals = mp.eig(W)[0]
        
        for ev in eigvals:
            if 0.99 < abs(ev) < 1.01 and abs(ev.imag) > 1e-50:
                angle = mp.arg(ev)
                # Return angle in radians, normalized to (-π, π]
                if angle > mp.pi:
                    angle -= 2 * mp.pi
                return angle
        
        return mp.mpf(0)  # Twist-free
    
    def is_twist_free(self, W: mp.matrix, tol: float = 1e-12) -> bool:
        """Check if word is twist-free (θ ≈ 0)."""
        return abs(self.get_twist_angle(W)) < tol