"""Repeat the full candidate DFS with outward-rounded interval arithmetic.

The analytic pruning and the tails are treated in Section 6 of the paper.
Run from the repository root: sage code/review_dfs_interval.sage
"""
load("code/review_dfs_bounds.sage")
import sys
from pathlib import Path
sys.setrecursionlimit(100000)
PR=list(prime_range(2,3000000))

@cached_function
def logp(p): return I(p).log()

@cached_function
def r1(p,e):
    v=(2+I(p).sqrt()) if p>=5 and e==1 else lam(p,e)
    return v/(p**e+p**(e-1))

def nu(p,e):
    return p**(e//2)+p**(e//2-1) if e%2==0 else 2*p**(e//2)

SF={}
for q in [0,2,3,5,7]:
    S=I(1)
    for p in [3,5,7]:
        if p>q: S*=max_interval([I(1)]+[r1(p,e)*(1+e*logp(p)/2) for e in range(1,60)])
    SF[q]=S

K=24/I.pi()
cands=[]
def G(C,R1,R2,R3): return K*C*R1+7*R2+6*R3

def dfs(N,C,R1,R2,R3,k):
    for idx in range(k,len(PR)):
        q=PR[idx]
        if q>=11:
            b=K*(C+logp(q))*R1*r1(q,1)+(7*R2+6*R3)*I(2)/(q+1)
            if b.upper()<1: return
            for e in range(1,60):
                qe=q**e
                if qe>10**40: break
                psi=qe+q**(e-1)
                C2=C+e*logp(q)
                A=R1*r1(q,e); B=R2*I(2)/psi; D=R3*I(nu(q,e))/psi
                if G(C2,A,B,D).upper()>=1:
                    cands.append(N*qe)
                    dfs(N*qe,C2,A,B,D,idx+1)
        else:
            for e in range(1,60):
                qe=q**e; psi=qe+q**(e-1)
                C2=C+e*logp(q)
                A=R1*r1(q,e); B=R2*I(2)/psi; D=R3*I(nu(q,e))/psi
                if (K*C2*A*SF[q]+7*B+6*D).upper()<1: continue
                if G(C2,A,B,D).upper()>=1: cands.append(N*qe)
                dfs(N*qe,C2,A,B,D,idx+1)
    raise RuntimeError("prime list exhausted")

print("Starting interval DFS.",flush=True)
dfs(ZZ(1),I(2),I(1),I(1),I(1),0)
cands=sorted(set(cands))
old=[ZZ(x) for x in Path("logs/cm_general_candidates.txt").read_text().split()]
extra=sorted(set(cands)-set(old))
missing=sorted(set(old)-set(cands))
print("Interval candidates:",len(cands),"maximum:",max(cands))
print("Additional candidates:",extra[:100],"count:",len(extra))
print("Old candidates not in interval list:",missing[:100],"count:",len(missing))
assert cands==old
print("DFS ALL CHECKS PASSED",flush=True)
