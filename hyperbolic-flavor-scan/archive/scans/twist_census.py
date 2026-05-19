"""
twist_census.py -- All twist angles up to length 5, m003 and m006
Checks coincidences with SM angles and mass ratios
"""
import snappy, numpy as np
from scipy.linalg import logm
import pandas as pd

def get_eigendata(G,word):
    try:
        mat=np.array([[complex(G.SL2C(word)[i,j]) for j in range(2)] for i in range(2)])
        eigs=np.linalg.eigvals(mat)
        lam=eigs[0] if np.imag(eigs[0])>=0 else eigs[1]
        log_lam=np.log(lam); t=np.real(log_lam); phi=np.imag(log_lam)
        return t,phi,np.exp(abs(t)),np.degrees(phi)
    except: return None,None,None,None

sm_angles={'theta12_CKM':13.04,'theta23_CKM':2.38,'theta13_CKM':0.201,
           'delta_CKM':68.0,'7xtheta12_CKM':91.28,'theta12_NU':33.41,
           'theta23_NU':49.1,'theta13_NU':8.54,'delta_PMNS':197.0,'90deg':90.0}
mass_ratios={'mt/mb':173/4.18,'mb/mc':4.18/1.27,'mc/ms':1.27/0.095,
             'mtau/mmu':1.777/0.106,'mmu/me':0.106/0.000511,'MZ/MW':91.19/80.38}

from itertools import product as iprod
for manifold_idx,manifold_name in [(1,'m003'),(43,'m006')]:
    M=snappy.OrientableClosedCensus[manifold_idx]
    G=M.fundamental_group()
    ngen=len(list(G.generators()))
    letters=[c for c in 'abcd'[:ngen]]+[c.upper() for c in 'abcd'[:ngen]]
    print(f'\n{"="*60}\n{M.name()} vol={float(M.volume()):.4f} H1={M.homology()}\n{"="*60}')

    records=[]
    for length in range(1,6):
        for combo in iprod(letters,repeat=length):
            w=''.join(combo)
            valid=all(w[i].lower()!=w[i+1].lower() or w[i]==w[i+1] for i in range(len(w)-1))
            if not valid: continue
            t,phi,mod_lam,phi_deg=get_eigendata(G,w)
            if t is None: continue
            records.append(dict(word=w,length=length,t=round(t,6),
                phi_deg=round(phi_deg,3),mod_lambda=round(mod_lam,6)))

    df=pd.DataFrame(records).drop_duplicates('word')
    df.to_csv(f'C:\\dev\\hyperbolic-flavor-scan\\twist_census_{manifold_name}.csv',index=False)
    print(f'Words computed: {len(df)}')

    print('\nSorted by |phi_deg| (first 25):')
    df2=df.copy(); df2['abs_phi']=df2['phi_deg'].abs()
    print(df2.sort_values('abs_phi')[['word','length','phi_deg','mod_lambda']].head(25).to_string(index=False))

    print('\nSM angle coincidences (within 3 deg):')
    for _,row in df.iterrows():
        phi=row['phi_deg']
        for sm_name,sm_val in sm_angles.items():
            for cand in [phi,phi%180,abs(phi),180-abs(phi),(phi+360)%360]:
                if abs(cand-sm_val)<3.0:
                    print(f"  {row['word']:8s} phi={phi:8.3f} ~ {sm_name}={sm_val} err={abs(cand-sm_val):.3f}")

    print('\nMass ratio coincidences (within 5%):')
    mods=df['mod_lambda'].values; wds=df['word'].values
    for i in range(len(mods)):
        for j in range(i+1,min(i+50,len(mods))):
            ratio=mods[i]/mods[j] if mods[j]>0 else 0
            for mr_name,mr_val in mass_ratios.items():
                if 0<ratio and abs(ratio/mr_val-1)<0.05:
                    print(f"  {wds[i]}/{wds[j]}={ratio:.4f} ~ {mr_name}={mr_val:.4f} err={abs(ratio/mr_val-1)*100:.1f}%")

print('\nCensus complete.')
