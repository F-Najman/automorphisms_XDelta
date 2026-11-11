# Candidates for failure of (P_N), general N (Section 6 of the paper: the finite set of levels).
# (P_N) holds whenever G(N) := (24/pi)(log N + 2) R1(N) + 7 R2(N) + 6 R3(N) < 1, where
#   R1 = prod Lam_p(e)/psi(p^e)  [M(N) <= ((log N + 2)/pi) prod Lam_p(e)],
#   R2 = prod 2/psi(p^e)  [nu_2, nu_3 <= 2^omega],  R3 = prod nu_p(e)/psi(p^e)  [nu_inf = prod nu_p].
# Exhaustive DFS over factorisations (primes increasing):
#  * primes q >= 11 never increase G: max_e r1(q,e)(1 + e log q/2) < 1 (asserted below);
#  * for q >= 11 the single-extension bound B(q) is decreasing in q, so the loop over q stops once B(q) < 1;
#  * for q in {2,3,5,7}: subtree bound (24/pi) C R1 S(q) + 7 R2 + 6 R3, where S(q) = prod over primes q < p < 11 of
#    max(1, max_e r1(p,e)(1 + e log p/2)); subtrees with bound < 1 are skipped.
# Output: logs/cm_general_candidates.txt (one N per line). Run from the repository root: python3 code/cm_general_dfs.py
import math, sys
sys.path.insert(0, 'code')
from cmlocal import Lam, psi_pe
from sympy import primerange

PR = list(primerange(2, 3 * 10**6))
EMAX = 60
lam_cache = {}
def r1(p, e):
    if (p, e) not in lam_cache: lam_cache[(p, e)] = Lam(p, e) / psi_pe(p, e)
    return lam_cache[(p, e)]
def nu_pe(p, e):
    tot = 0
    for j in range(e + 1):
        m = min(j, e - j)
        tot += 1 if m == 0 else p**m - p**(m - 1)
    return tot

for q in [11, 13, 17, 19, 23, 29, 31, 37]:
    assert max(r1(q, e) * (1 + e * math.log(q) / 2) for e in range(1, 12)) < 1, q
    assert all(r1(q, e) * (2 + e * math.log(q)) <= r1(q, 1) * (2 + math.log(q)) for e in range(2, 12)), q
def Sfac(q):
    s = 1.0
    for p in [3, 5, 7]:
        if p > q:
            s *= max(1.0, max(r1(p, e) * (1 + e * math.log(p) / 2) for e in range(1, EMAX)))
    return s
SF = {q: Sfac(q) for q in [0, 2, 3, 5, 7]}
K = 24 / math.pi
cands = []
def G(C, R1, R2, R3): return K * C * R1 + 7 * R2 + 6 * R3

def dfs(N, C, R1, R2, R3, k):
    for idx in range(k, len(PR)):
        q = PR[idx]
        if q >= 11:
            b = K * (C + math.log(q)) * R1 * r1(q, 1) + (7 * R2 + 6 * R3) * 2 / (q + 1)
            if b < 1: return
            for e in range(1, EMAX):
                if q**e > 10**40: break
                C2 = C + e * math.log(q); R1b = R1 * r1(q, e); R2b = R2 * 2 / psi_pe(q, e); R3b = R3 * nu_pe(q, e) / psi_pe(q, e)
                if G(C2, R1b, R2b, R3b) >= 1:
                    cands.append(N * q**e); dfs(N * q**e, C2, R1b, R2b, R3b, idx + 1)
        else:
            for e in range(1, EMAX):
                C2 = C + e * math.log(q); R1b = R1 * r1(q, e); R2b = R2 * 2 / psi_pe(q, e); R3b = R3 * nu_pe(q, e) / psi_pe(q, e)
                if K * C2 * R1b * SF[q] + 7 * R2b + 6 * R3b < 1: continue
                if G(C2, R1b, R2b, R3b) >= 1: cands.append(N * q**e)
                dfs(N * q**e, C2, R1b, R2b, R3b, idx + 1)
    raise RuntimeError("prime list too short")

sys.setrecursionlimit(100000)
dfs(1, 2.0, 1.0, 1.0, 1.0, 0)
cands = sorted(set(cands))
open('logs/cm_general_candidates.txt', 'w').write('\n'.join(map(str, cands)) + '\n')
print("CANDIDATES", len(cands), "max", max(cands), "S factors", SF)
import collections
print("by digits", sorted(collections.Counter(len(str(c)) for c in cands).items()))
print("DONE")
