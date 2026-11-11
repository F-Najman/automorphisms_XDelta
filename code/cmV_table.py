# Table of the CM bound V(N, t) for every local type tuple t at level N (primes in increasing order).
# Usage: python3 code/cmV_table.py N > file ;  line format: "code_1 ... code_r V"
import sys, itertools
sys.path.insert(0, 'code')
from cmV import V
from sympy import factorint
N = int(sys.argv[1]); fac = dict(sorted(factorint(N).items()))
def types(p, e):
    if p == 2:
        ts = [0] + ([1] if e >= 2 else [])
        for k in range(1, e - 1): ts += [2 * k, 2 * k + 1]
        return ts
    ts = [0, 3] + ([2] if p > 3 else []) + [2 * a for a in range(2, e + 1)]
    return ts
ps = list(fac)
for tup in itertools.product(*[types(p, fac[p]) for p in ps]):
    print(" ".join(map(str, tup)), "%.6f" % V(fac, dict(zip(ps, tup))))
