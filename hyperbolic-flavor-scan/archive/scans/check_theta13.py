import pandas as pd
import numpy as np

# Load length-7 results for m006
df = pd.read_csv(r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len7_m006.csv")

target = 0.201  # theta13_CKM in degrees
tol    = 0.5    # within 0.5 deg

print(f"theta13_CKM = {target} deg")
print(f"\nAll words with phi_fold within {tol} deg of target:")
near = df[abs(df.phi_fold - target) < tol].sort_values("phi_fold")
if len(near):
    print(near[["word","length","phi_deg","phi_fold","mod_lambda"]].to_string())
else:
    print("  None found at length <= 7.")

print(f"\nAll phi_fold < 0.5 deg (sorted):")
low = df[df.phi_fold < 0.5].sort_values("phi_fold")
print(low[["word","length","phi_deg","phi_fold"]].to_string())

print(f"\nAll phi_fold between 0.05 and 0.40 deg:")
mid = df[(df.phi_fold > 0.05) & (df.phi_fold < 0.40)].sort_values("phi_fold")
if len(mid):
    print(mid[["word","length","phi_deg","phi_fold"]].to_string())
else:
    print("  None -- spectral gap between 0.005 and next value")

# What IS the next phi_fold above 0.005?
print(f"\nSpectral values sorted (first 20):")
print(df.drop_duplicates("phi_fold").sort_values("phi_fold").head(20)[
    ["word","length","phi_fold"]].to_string())
