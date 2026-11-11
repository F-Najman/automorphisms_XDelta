"""Exact finite gonality check (the lemma on gonality in Section 7 of the paper).
Run from the repository root: sage code/review_gonality.sage

If the least prime not dividing N is at least 13, the primorial lower
bound already suffices at index 1. Otherwise a failure requires N<=2742.
"""
from pathlib import Path
import re

def psi(n):
    return prod(p**(e-1)*(p+1) for p,e in factor(n))

def least_missing(n):
    ell=ZZ(2)
    while n%ell==0: ell=next_prime(ell)
    return ell

assert 7*psi(2*3*5*7*11)>1600*(13+1)
bad={i:[] for i in [1,2,4]}
for n in range(201,2743):
    ell=least_missing(n)
    for idx in bad:
        if 7*psi(n)*idx<=1600*(ell+1):
            if idx==4 or n>500: bad[idx].append(n)
assert len(bad[1])==105 and max(bad[1])==840
assert bad[2]==[] and bad[4]==[210]
print("Index 1 failures above 500:",len(bad[1]),"; maximum",max(bad[1]))
print("Index 2 failures above 500:",bad[2])
print("Index 4 failures above 200:",bad[4])

cmfails={int(x) for x in Path("logs/cm_exact_fails.txt").read_text().split() if int(x)>500}
needed=set(bad[1])|cmfails
covered=set()
for p in Path("logs/xdelta").glob("*.log"):
    text=p.read_text()
    if "DONE" not in text: continue
    for n in re.findall(r"DELTA N=(\d+) gens=\[[^\]]*\] idx=1 [^\n]* OK",text):
        covered.add(int(n))
extra={int(x) for x in Path("logs/x0_levels.txt").read_text().split()}
assert needed-covered==extra
print("Exactly",len(extra),"additional X_0 certificates are required and listed.")
print("ALL CHECKS PASSED")
