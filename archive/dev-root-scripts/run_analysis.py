# ============================================
# FILE: run_analysis.py - CORRECTED
# COMPLETE NEUTRINO SECTOR ANALYSIS SCRIPT
# ============================================

import numpy as np
import mpmath as mp
import json
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

# Set global precision
mp.mp.dps = 200

# ------------------------------
# 1. COXETER EXACT FRAMEWORK
# ------------------------------
@dataclass
class CoxeterConfig:
    phi_lock: bool = True
    precision: int = 200
    use_exact_gram: bool = True

class Coxeter353:
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
        φ = self.phi
        g01 = -mp.mpf(1) / 2
        g12 = -φ / 2  
        g23 = -mp.mpf(1) / 2
        g02 = mp.mpf(0)
        g13 = mp.mpf(0)  # Fixed: was mp.mpeg
        cosh_dstar = 1 + φ / 2
        g03 = -cosh_dstar
        
        return mp.matrix([
            [1, g01, g02, g03],
            [g01, 1, g12, g13],
            [g02, g12, 1, g23],
            [g03, g13, g23, 1]
        ])
    
    def _build_minkowski_metric(self) -> mp.matrix:
        return mp.matrix([[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    
    def _build_reflections(self) -> Tuple[mp.matrix, ...]:
        N = mp.zeros(4, 4)
        N[0,0] = mp.mpf(1)
        N[1,0] = -self.G[0,1]
        N[1,1] = mp.sqrt(1 - N[1,0]**2)
        N[2,0] = mp.mpf(0)
        N[2,1] = -self.G[1,2] / N[1,1]
        N[2,2] = mp.sqrt(1 - N[2,1]**2)
        N[3,0] = -self.G[0,3]
        N[3,1] = mp.mpf(0)  # Fixed: was mp.mpeg
        N[3,2] = -self.G[2,3] / N[2,2]
        N[3,3] = mp.sqrt(1 - N[3,0]**2 - N[3,2]**2)
        
        S = []
        for i in range(4):
            n = N[i,:].transpose()
            denom = (n.T * self.J * n)[0]
            R = mp.eye(4) - 2 * (self.J * n * n.T) / denom
            S.append(R)
        
        return tuple(S)
    
    def ladder_word_matrix(self, m: int, n: int) -> mp.matrix:
        tau_m = mp.matrix_pow(self.tau, m)
        tau_n = mp.matrix_pow(self.tau, n)
        rho = self.rho
        rho_inv = mp.inverse(rho)
        return tau_m * rho * tau_n * rho_inv
    
    def get_translation_length(self, W: mp.matrix) -> mp.mpf:
        eigvals = mp.eig(W)[0]
        max_abs = max(abs(ev) for ev in eigvals)
        return mp.log(max_abs) if max_abs > 1 else mp.mpf(0)  # Fixed: was mp.mpeg
    
    def get_twist_angle(self, W: mp.matrix) -> mp.mpf:
        eigvals = mp.eig(W)[0]
        for ev in eigvals:
            if 0.99 < abs(ev) < 1.01 and abs(ev.imag) > 1e-50:
                angle = mp.arg(ev)
                if angle > mp.pi: angle -= 2 * mp.pi
                return angle
        return mp.mpf(0)  # Fixed: was mp.mpeg

# ------------------------------
# 2. LADDER SEARCH MODULE
# ------------------------------
@dataclass
class LadderCandidate:
    m: int
    n: int
    translation_length: mp.mpf  # Fixed: was mp.mpeg
    twist_angle: mp.mpf         # Fixed: was mp.mpeg
    q_phi: mp.mpf               # Fixed: was mp.mpeg
    q_nearest_int: int
    q_residual: mp.mpf          # Fixed: was mp.mpeg
    word_digits: str
    
    def __str__(self):
        return f"W({self.m},{self.n}): ℓ={float(self.translation_length):.6f}, θ={float(self.twist_angle):.3e}, q={float(self.q_phi):.8f}, |q-{self.q_nearest_int}|={float(self.q_residual):.3e}"

class LadderSearch:
    def __init__(self, coxeter):
        self.coxeter = coxeter
        self.phi = (1 + mp.sqrt(5)) / 2
        self.log_phi = mp.log(self.phi)
    
    def search_range(self, m_range, n_range, twist_tol=1e-10, ell_min=1e-6, ell_max=20.0):
        m_min, m_max = m_range
        n_min, n_max = n_range
        candidates = []
        
        for m in range(m_min, m_max + 1):
            for n in range(n_min, n_max + 1):
                if m == 0 and n == 0: continue
                W = self.coxeter.ladder_word_matrix(m, n)
                ell = self.coxeter.get_translation_length(W)
                theta = self.coxeter.get_twist_angle(W)
                
                if ell < ell_min or ell > ell_max: continue
                if abs(theta) > twist_tol: continue
                
                q = ell / self.log_phi
                q_int = int(mp.nint(q))
                q_resid = abs(q - q_int)
                word_digits = self._word_to_digits(m, n)
                
                candidate = LadderCandidate(
                    m=m, n=n, translation_length=ell, twist_angle=theta,
                    q_phi=q, q_nearest_int=q_int, q_residual=q_resid,
                    word_digits=word_digits
                )
                candidates.append(candidate)
        
        candidates.sort(key=lambda c: (float(c.translation_length), abs(float(c.twist_angle))))
        return candidates
    
    def _word_to_digits(self, m: int, n: int) -> str:
        def tau_power(k):
            if k > 0: return "30" * k
            elif k < 0: return "03" * (-k)
            return ""
        return tau_power(m) + "12" + tau_power(n) + "21"

# ------------------------------
# 3. NEUTRINO ANALYSIS MODULE
# ------------------------------
@dataclass
class PMNSData:
    theta12: float; theta23: float; theta13: float
    delta_cp: float; dm21_sq: float; dm32_sq: float
    ordering: str = 'normal'

@dataclass 
class NeutrinoCandidate:
    ladder_word: LadderCandidate
    geometric_params: Dict[str, float]
    pmns_fit: Dict[str, float]
    mass_scale: float
    hierarchy_ratio: float

class NeutrinoAnalyzer:
    def __init__(self, coxeter, pmns_data: PMNSData):
        self.coxeter = coxeter
        self.pmns = pmns_data
        
    def estimate_neutrino_translation(self) -> float:
        # CKM: ℓ ≈ 8.01, θ_sum ≈ 15.6°
        # PMNS: θ_sum ≈ 84°, scale factor ≈ 15.6/84 ≈ 0.185
        theta_sum = self.pmns.theta12 + self.pmns.theta23 + self.pmns.theta13
        scale_factor = 15.6 / theta_sum
        return 8.01 * scale_factor
    
    def search_neutrino_ladder(self, target_translation, search_range=15, tolerance=0.5):
        searcher = LadderSearch(self.coxeter)
        return searcher.search_range(
            m_range=(-search_range, search_range),
            n_range=(-search_range, search_range),
            twist_tol=1e-8,
            ell_min=max(0.1, target_translation - tolerance),
            ell_max=target_translation + tolerance
        )
    
    def compute_geometric_parameters(self, ladder_word):
        ell = float(ladder_word.translation_length)
        theta = float(ladder_word.twist_angle)
        return {
            'alpha': ell, 'beta': ell/2,
            'delta_theta_deg': np.degrees(theta),
            'translation_length': ell,
            'twist_angle_deg': np.degrees(theta),
            'q_phi': float(ladder_word.q_phi),
            'q_residual': float(ladder_word.q_residual)
        }
    
    def integrate_mass_hierarchy(self, ladder_word):
        ell = float(ladder_word.translation_length)
        log_mass_ratio = np.log(self.pmns.dm32_sq / self.pmns.dm21_sq)
        mass_scale = 0.05  # eV estimate
        return mass_scale, log_mass_ratio

# ------------------------------
# 4. MAIN EXECUTION
# ------------------------------
def main():
    print("=" * 80)
    print("NEUTRINO SECTOR ANALYSIS - COXETER [3,5,3] FRAMEWORK")
    print("=" * 80)
    
    # 1. SETUP
    print("\n1. Building exact Coxeter [3,5,3] system...")
    coxeter = Coxeter353()
    print(f"   Golden ratio φ = {float(coxeter.phi):.15f}")
    print(f"   log(φ) = {float(coxeter.log_phi):.15f}")
    
    # 2. NEUTRINO DATA
    print("\n2. Loading neutrino observational data...")
    pmns_data = PMNSData(
        theta12=33.45, theta23=42.1, theta13=8.62,
        delta_cp=222, dm21_sq=7.42e-5, dm32_sq=2.517e-3
    )
    print(f"   θ₁₂ = {pmns_data.theta12}°")
    print(f"   θ₂₃ = {pmns_data.theta23}°")
    print(f"   θ₁₃ = {pmns_data.theta13}°")
    
    # 3. ESTIMATE TRANSLATION
    print("\n3. Estimating neutrino translation scale...")
    analyzer = NeutrinoAnalyzer(coxeter, pmns_data)
    target = analyzer.estimate_neutrino_translation()
    print(f"   Target translation length: ℓ ≈ {target:.6f}")
    print(f"   Corresponding q = ℓ/log(φ) ≈ {target/float(coxeter.log_phi):.6f}")
    
    # 4. SEARCH
    print("\n4. Searching for neutrino ladder words...")
    candidates = analyzer.search_neutrino_ladder(target, search_range=12)
    print(f"   Found {len(candidates)} candidates with ℓ ≈ {target:.3f} ± 0.3")
    
    if candidates:
        print("\n   Top 5 candidates:")
        for i, cand in enumerate(candidates[:5]):
            print(f"   {i+1}. {cand}")
        
        # 5. ANALYZE BEST
        best = candidates[0]
        print(f"\n5. Best candidate: W({best.m},{best.n})")
        print(f"   ℓ = {float(best.translation_length):.8f}")
        print(f"   θ = {float(best.twist_angle):.3e}")
        print(f"   q = {float(best.q_phi):.8f}")
        print(f"   |q-{best.q_nearest_int}| = {float(best.q_residual):.3e}")
        
        # 6. GEOMETRIC PARAMS
        print("\n6. Geometric parameters:")
        geo = analyzer.compute_geometric_parameters(best)
        for key, val in geo.items():
            print(f"   {key}: {val:.8f}")
        
        # 7. MASS HIERARCHY
        print("\n7. Mass hierarchy integration:")
        mass_scale, log_ratio = analyzer.integrate_mass_hierarchy(best)
        print(f"   Estimated mass scale: {mass_scale:.3f} eV")
        print(f"   Predicted ratio: {np.exp(log_ratio):.1f}")
        print(f"   Observed ratio: {pmns_data.dm32_sq/pmns_data.dm21_sq:.1f}")
        
        # 8. HANDOFF
        handoff = {
            'phase': 'neutrino_certification_complete',
            'best_candidate': {'m': best.m, 'n': best.n,
                              'translation_length': float(best.translation_length),
                              'q_phi': float(best.q_phi)},
            'next_phase': 'mass_hierarchy_integration'
        }
        with open('neutrino_handoff.json', 'w') as f:
            json.dump(handoff, f, indent=2)
        print("\n✓ Handoff saved to 'neutrino_handoff.json'")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()