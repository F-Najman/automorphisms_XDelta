# Number of subgroups of (Z/N)^x containing -1 = number of subgroups of A = (Z/N)^x/<-1>, via Birkhoff's formula
# for each p-primary part (checked against brute force on small groups). Compared with logs/subgroups_magma.log.
import itertools, math, sys
from sympy import factorint, primerange
def unit_group_invariants(N):
    # elementary divisors (p-power cyclic factors) of (Z/N)^x / <-1>, computed from an explicit presentation
    # by brute force on the group: orders of elements via the Smith form of the relation lattice.
    from sympy import Matrix
    from sympy.matrices.normalforms import smith_normal_form
    fac = factorint(N); cyc = []; gens = []
    # generators of (Z/p^e)^x and their orders; -1 expressed in these coordinates
    comps = []
    for p, e in fac.items():
        if p == 2:
            if e == 1: continue
            if e == 2: comps.append([(2, 'm1')]); continue
            comps.append([(2, 'm1'), (2**(e - 2), 'five')])
        else:
            comps.append([((p - 1) * p**(e - 1), 'g')])
    orders = [o for c in comps for (o, _) in c]
    # coordinates of -1: for p=2: (1, 0) on (-1, 5); for odd p: -1 = g^(ord/2)
    m1 = []
    for c in comps:
        for (o, kind) in c:
            m1.append(1 if kind == 'm1' else (0 if kind == 'five' else o // 2))
    n = len(orders)
    if n == 0: return []
    rel = [[orders[i] if j == i else 0 for j in range(n)] for i in range(n)] + [m1]
    S = smith_normal_form(Matrix(rel))
    d = [abs(S[i, i]) for i in range(min(S.shape)) if S[i, i] != 0]
    out = []
    for x in d:
        for p, e in factorint(x).items(): out.append((p, e))
    return out
def gauss_binom(n, k, p):
    if k < 0 or k > n: return 0
    num = den = 1
    for i in range(k):
        num *= p**(n - i) - 1; den *= p**(i + 1) - 1
    return num // den
def conj(part):
    return [sum(1 for x in part if x > i) for i in range(max(part) if part else 0)]
def n_sub_type(lam, mu, p):
    lc, mc = conj(lam), conj(mu)
    L = max(len(lc), 1)
    lc += [0] * (L + 2 - len(lc)); mc += [0] * (L + 2 - len(mc))
    tot = 1
    for i in range(L):
        tot *= p**(mc[i + 1] * (lc[i] - mc[i])) * gauss_binom(lc[i] - mc[i + 1], mc[i] - mc[i + 1], p)
    return tot
def n_sub_pgroup(lam, p):
    lam = sorted(lam, reverse=True); tot = 0
    for mu in itertools.product(*[range(x + 1) for x in lam]):
        mu = sorted([m for m in mu], reverse=True)
        # count each partition mu <= lam once
        pass
    seen = set()
    for mu in itertools.product(*[range(x + 1) for x in lam]):
        m = tuple(sorted([x for x in mu if x > 0], reverse=True))
        if m in seen: continue
        if all(m[i] <= lam[i] for i in range(len(m))): seen.add(m); tot += n_sub_type(lam, list(m), p)
    return tot
def brute(orders):
    # brute-force subgroup count of Z/o1 x ... x Z/ok
    elems = list(itertools.product(*[range(o) for o in orders]))
    def add(a, b): return tuple((x + y) % o for x, y, o in zip(a, b, orders))
    zero = tuple(0 for _ in orders)
    def closure(gs):
        H = {zero}; frontier = [zero]
        while frontier:
            new = []
            for h in frontier:
                for g in gs:
                    x = add(h, g)
                    if x not in H: H.add(x); new.append(x)
            frontier = new
        return frozenset(H)
    subs = {frozenset([zero])}; frontier = list(subs)
    while frontier:
        new = []
        for H in frontier:
            for g in elems:
                if g not in H:
                    K = closure(list(H) + [g]) if len(H) < 50 else None
                    if K is None: continue
                    if K not in subs: subs.add(K); new.append(K)
        frontier = new
    return len(subs)
# self-test of Birkhoff against brute force
for orders, p, lam in [([2, 4], 2, [2, 1]), ([2, 2, 2], 2, [1, 1, 1]), ([4, 8], 2, [3, 2]), ([3, 9], 3, [2, 1]), ([2, 2, 4], 2, [2, 1, 1])]:
    assert brute(orders) == n_sub_pgroup(lam, p), (orders, brute(orders), n_sub_pgroup(lam, p))
print("Birkhoff self-test OK")
Ns = list(range(13, 501)) + [504, 512, 540, 567, 576, 600, 640, 648, 675, 720, 756, 768, 784, 800, 864, 900, 972, 1024, 1080, 1152, 1296, 1440, 1536, 1728, 2304, 2592, 3072, 3456, 4608, 6912]
res = {}
for N in Ns:
    inv = unit_group_invariants(N); byp = {}
    for p, e in inv: byp.setdefault(p, []).append(e)
    c = 1
    for p, lam in byp.items(): c *= n_sub_pgroup(lam, p)
    res[N] = c
mag = {}
for line in open('logs/subgroups_magma.log'):
    if line.startswith('SUB'):
        w = dict(x.split('=') for x in line.split()[1:]); mag[int(w['N'])] = int(w['count'])
bad = [(N, res[N], mag.get(N)) for N in Ns if mag.get(N) != res[N]]
print("compared", len([N for N in Ns if N in mag]), "levels; mismatches:", bad[:20])
