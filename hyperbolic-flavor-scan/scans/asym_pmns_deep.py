import numpy as np
from scipy.linalg import qr
from scipy.optimize import minimize
import time

PMNS = np.array([[0.821,0.550,0.148],[0.357,0.339,0.871],[0.442,0.762,0.471]])

def U_qr(O):
    U, _ = qr(O)
    for i in range(3):
        if U[i,i] < 0: U[:,i] *= -1
    return np.abs(U)

def is_pd(O):
    return np.all(np.linalg.eigvalsh(O) > 1e-10)

def fitness_asym(p):
    O = np.array([[1.0,  p[0], p[1]],
                  [p[2], 1.0,  p[3]],
                  [p[4], p[5], 1.0 ]])
    if not is_pd(O): return 1e9
    return float(np.linalg.norm(U_qr(O) - PMNS, 'fro'))

np.random.seed(42)
N = 5000
best, best_p = 1e9, None
t0 = time.time()

print(f'Running {N} Nelder-Mead starts on asymmetric overlap matrix...')
print(f'CKM benchmark: 0.017290')
print(f'Previous best: 0.018968 (500 iters)')
print('='*60)

for i in range(N):
    p0 = np.random.uniform(0.01, 0.99, 6)
    res = minimize(fitness_asym, p0, method='Nelder-Mead',
                   options={'xatol':1e-9,'fatol':1e-9,'maxiter':3000})
    if res.fun < best:
        best, best_p = res.fun, res.x
        O = np.array([[1.0,      best_p[0], best_p[1]],
                      [best_p[2],1.0,       best_p[3]],
                      [best_p[4],best_p[5], 1.0      ]])
        U = U_qr(O)
        elapsed = time.time() - t0
        print(f'[iter {i:4d} | {elapsed:5.1f}s] NEW BEST: {best:.6f}')
        print(f'  Overlap matrix O:')
        print(f'    [ 1.000  {best_p[0]:.4f}  {best_p[1]:.4f} ]')
        print(f'    [{best_p[2]:.4f}   1.000  {best_p[3]:.4f} ]')
        print(f'    [{best_p[4]:.4f}  {best_p[5]:.4f}   1.000 ]')
        print(f'  |U| found:')
        print(f'    {U[0,0]:.4f}  {U[0,1]:.4f}  {U[0,2]:.4f}')
        print(f'    {U[1,0]:.4f}  {U[1,1]:.4f}  {U[1,2]:.4f}')
        print(f'    {U[2,0]:.4f}  {U[2,1]:.4f}  {U[2,2]:.4f}')
        print(f'  PMNS target:')
        print(f'    0.821   0.550   0.148')
        print(f'    0.357   0.339   0.871')
        print(f'    0.442   0.762   0.471')
        print(f'  Row residuals: {np.linalg.norm(U[0]-PMNS[0]):.4f} / {np.linalg.norm(U[1]-PMNS[1]):.4f} / {np.linalg.norm(U[2]-PMNS[2]):.4f}')
        print()

    if (i+1) % 500 == 0:
        elapsed = time.time() - t0
        rate = (i+1)/elapsed
        eta = (N-i-1)/rate
        print(f'--- Progress: {i+1}/{N} | best={best:.6f} | {rate:.0f} it/s | ETA {eta:.0f}s ---')

print()
print('='*60)
print('FINAL RESULT')
print('='*60)
print(f'Best fitness:  {best:.6f}')
print(f'CKM benchmark: 0.017290')
print(f'Improvement over previous: {0.018968 - best:.6f}')
O = np.array([[1.0,      best_p[0], best_p[1]],
              [best_p[2],1.0,       best_p[3]],
              [best_p[4],best_p[5], 1.0      ]])
U = U_qr(O)
print(f'Optimal overlap matrix:')
print(f'  [ 1.000  {best_p[0]:.6f}  {best_p[1]:.6f} ]')
print(f'  [{best_p[2]:.6f}   1.000  {best_p[3]:.6f} ]')
print(f'  [{best_p[4]:.6f}  {best_p[5]:.6f}   1.000 ]')
print(f'Final |U|:')
for row, trow in zip(U, PMNS):
    print(f'  {row[0]:.4f}  {row[1]:.4f}  {row[2]:.4f}   (target: {trow[0]:.3f}  {trow[1]:.3f}  {trow[2]:.3f})')