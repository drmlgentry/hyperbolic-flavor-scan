"""
scan_len3x3.py -- Complex Borel len3+len3+len3 scan on m003
Target: simultaneous F<0.020 AND J near J_PDG=-0.00966
Checkpointed output: scan_len3x3_results.csv, scan_len3x3_checkpoint.txt
"""
import snappy, numpy as np
from scipy.linalg import logm, qr
from scipy.optimize import minimize
from itertools import permutations
import pandas as pd, time, os

PMNS=np.array([[0.821,0.550,0.148],[0.357,0.339,0.871],[0.442,0.762,0.471]])
J_PDG=-0.00966; F_TARGET=0.025
ALPHAS=[0.85,0.90,0.95,1.00,1.05,1.10,1.125,1.15]
OUT_CSV=r'C:\dev\hyperbolic-flavor-scan\scan_len3x3_results.csv'
CKPT=r'C:\dev\hyperbolic-flavor-scan\scan_len3x3_checkpoint.txt'

M=snappy.OrientableClosedCensus[1]; G=M.fundamental_group()
print('Manifold:',M.name(),'vol:',M.volume(),'H1:',M.homology())
words3=[a+b+c for a in 'abAB' for b in 'abAB' for c in 'abAB']
print(f'Length-3 words: {len(words3)}')

def get_axis(word,alpha):
    mat=np.array([[complex(G.SL2C(word)[i,j]) for j in range(2)] for i in range(2)])
    L=logm(mat)
    nr=np.array([np.real(L[0,1]+L[1,0]),np.real(1j*(L[1,0]-L[0,1])),np.real(L[0,0]-L[1,1])])
    ni=np.array([np.imag(L[0,1]+L[1,0]),np.imag(1j*(L[1,0]-L[0,1])),np.imag(L[0,0]-L[1,1])])
    nc=nr+alpha*1j*ni; nm=np.sqrt(np.real(np.dot(nc,np.conj(nc))))
    return nc/nm if nm>1e-10 else nc

def evaluate(axes,scale):
    l21=scale*np.dot(np.conj(axes[0]),axes[1])
    l31=scale*np.dot(np.conj(axes[0]),axes[2])
    l32=scale*np.dot(np.conj(axes[1]),axes[2])
    L=np.array([[1+0j,0,0],[l21,1+0j,0],[l31,l32,1+0j]])
    Q,R=qr(L); d=np.diag(R); Q=Q*(d/np.abs(d))[np.newaxis,:]
    U=np.abs(Q); best_f=np.inf; best_p=None
    for p in permutations(range(3)):
        f=np.linalg.norm(U[:,list(p)]-PMNS,'fro')
        if f<best_f: best_f=f; best_p=p
    Qp=Q[:,list(best_p)]
    J=np.imag(Qp[0,0]*Qp[1,1]*np.conj(Qp[0,1])*np.conj(Qp[1,0]))
    return best_f,J

def optimize(words,alpha):
    axes=[get_axis(w,alpha) for w in words]
    def loss(p):
        try: f,J=evaluate(axes,abs(p[0])); return f**2+20*(J-J_PDG)**2
        except: return 99.
    best=(99.,0.,2.0)
    for sc in np.linspace(0.5,5.0,30):
        try:
            f,J=evaluate(axes,sc)
            if f+3*abs(J-J_PDG)<best[0]+3*abs(best[1]-J_PDG): best=(f,J,sc)
        except: pass
    try:
        res=minimize(loss,[best[2]],method='Nelder-Mead',options={'maxiter':800,'xatol':1e-7})
        f,J=evaluate(axes,abs(res.x[0])); return f,J,abs(res.x[0])
    except: return best

start_idx=0
if os.path.exists(CKPT):
    with open(CKPT) as f: start_idx=int(f.read().strip())
    print(f'Resuming from index {start_idx}')

results=[]
if os.path.exists(OUT_CSV):
    results=pd.read_csv(OUT_CSV).to_dict('records')
    print(f'Loaded {len(results)} existing results')

t0=time.time(); n=0; n_saved=0
for w1 in words3:
    for w2 in words3:
        if w2<=w1: continue
        for w3 in words3:
            if w3<=w2: continue
            if n<start_idx: n+=1; continue
            n+=1
            if n%500==0:
                elapsed=time.time()-t0
                print(f'  [{n}] {elapsed:.0f}s saved={n_saved}')
                with open(CKPT,'w') as f: f.write(str(n))
                if results: pd.DataFrame(results).sort_values('joint_loss').to_csv(OUT_CSV,index=False)
            for alpha in ALPHAS:
                try:
                    f,J,sc=optimize([w1,w2,w3],alpha)
                    if f<F_TARGET:
                        joint=f+3*abs(J-J_PDG)
                        results.append(dict(w1=w1,w2=w2,w3=w3,alpha=round(alpha,4),
                            scale=round(sc,5),fitness=round(f,6),J=round(J,6),
                            J_err=round(abs(J-J_PDG),6),joint_loss=round(joint,6)))
                        n_saved+=1
                        if f<0.021 and abs(J-J_PDG)<0.003:
                            print(f'*** EXCELLENT: {w1}/{w2}/{w3} a={alpha} F={f:.5f} J={J:.5f}')
                        elif f<0.021: print(f'** GOOD F: {w1}/{w2}/{w3} a={alpha} F={f:.5f} J={J:.5f}')
                        elif abs(J-J_PDG)<0.002: print(f'** GOOD J: {w1}/{w2}/{w3} a={alpha} F={f:.5f} J={J:.5f}')
                except: pass

with open(CKPT,'w') as f: f.write(str(n))
df=pd.DataFrame(results)
if len(df)>0:
    df=df.sort_values('joint_loss'); df.to_csv(OUT_CSV,index=False)
    print(f'\nDone. {n} triples, {n_saved} saved.')
    print(df.head(20).to_string(index=False))
else:
    print('No results found.')
