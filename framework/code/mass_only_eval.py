import json, math, pathlib
from dataclasses import dataclass

PHI = (1+5**0.5)/2
ALPHA = math.log(PHI)/4.0  # per your contract

WEIGHTS = {"a":8, "b":15, "c":24}

PDG_MZ = {
  # GeV (example placeholders; replace with your canonical PDG import)
  "m_u": (0.00216, 0.00026),
  "m_c": (1.27,   0.02),
  "m_t": (172.76, 0.30),
  "m_d": (0.00467,0.00048),
  "m_s": (0.093,  0.011),
  "m_b": (4.18,   0.03),
  "m_e": (0.000511, 2.5e-8),
  "m_mu":(0.10566, 3.5e-7),
  "m_tau":(1.77686, 1.2e-4),
}

def sector(f):
  if f in ("m_u","m_c","m_t"): return "up"
  if f in ("m_d","m_s","m_b"): return "down"
  return "lep"

def W_of(a,b,c): return 8*a + 15*b + 24*c

def all_sites(Wmax=200):
  # small deterministic search domain; expand later
  for a in range(0, Wmax//8 + 1):
    for b in range(0, Wmax//15 + 1):
      for c in range(0, Wmax//24 + 1):
        yield (a,b,c,W_of(a,b,c))

def lexicographic_assign(m_obs, m0, Wmax=200):
  """
  Choose site that minimizes |log m_obs - log(m0) - alpha*W|, tie-break by W.
  """
  target = math.log(m_obs) - math.log(m0)
  best = None
  for a,b,c,W in all_sites(Wmax=Wmax):
    err = abs(target - ALPHA*W)
    key = (err, W)
    if (best is None) or (key < best[0]):
      best = (key, (a,b,c,W))
  return best[1]

def fit_m0_for_sector(obs_list, W_list):
  """
  In log-space: log m_i = log m0 + alpha*W_i  => log m0 = mean(log m_i - alpha*W_i)
  """
  xs = [math.log(m) - ALPHA*W for (m,_sig),W in zip(obs_list,W_list)]
  return math.exp(sum(xs)/len(xs))

def chi2_mass(assignments, sigma_model_frac=0.05):
  chi2=0.0
  for f,(m_obs,sig) in PDG_MZ.items():
    a,b,c,W = assignments[f]
    m0 = assignments[f"_m0_{sector(f)}"]
    m_pred = m0*math.exp(ALPHA*W)
    sig2 = sig*sig + (sigma_model_frac*m_obs)**2
    chi2 += (m_obs-m_pred)**2/sig2
  return chi2

def main():
  in_jsonl = pathlib.Path("out/ladders_output.jsonl")
  if not in_jsonl.exists():
    raise SystemExit("Missing out/ladders_output.jsonl. Run scripts 01 and 02 first.")

  # For now, we just show framework evaluation independent of ladders:
  # ladders become an admissibility filter later (cert_ok gate).
  # This keeps the first pass deterministic and debuggable.

  # Initialize m0 guesses (deterministic)
  m0 = {"up":1e-3, "down":1e-3, "lep":1e-3}

  # One deterministic refinement loop
  for _ in range(3):
    # Assign sites with current m0 per sector
    assignments = {}
    for f,(m_obs,_sig) in PDG_MZ.items():
      assignments[f] = lexicographic_assign(m_obs=m_obs, m0=m0[sector(f)], Wmax=260)

    # Fit m0 per sector given assigned W
    for sec in ("up","down","lep"):
      obs = [(PDG_MZ[f][0], PDG_MZ[f][1]) for f in PDG_MZ if sector(f)==sec]
      Ws  = [assignments[f][3] for f in PDG_MZ if sector(f)==sec]
      m0[sec] = fit_m0_for_sector(obs, Ws)
      assignments[f"_m0_{sec}"] = m0[sec]

  out = {
    "alpha": ALPHA,
    "weights": WEIGHTS,
    "m0": m0,
    "assignments": {k:v for k,v in assignments.items() if not k.startswith("_m0_")},
    "chi2_mass": chi2_mass(assignments),
  }
  pathlib.Path("out/framework_mass_only.json").write_text(json.dumps(out, indent=2, sort_keys=True))
  print("Wrote out/framework_mass_only.json")
  print("chi2_mass =", out["chi2_mass"])

if __name__ == "__main__":
  main()
