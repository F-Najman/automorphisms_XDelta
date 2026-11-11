# Python port of the exact-conductor CM bound of code/cm_exact.c, for ONE character type tuple.
# types: dict p -> code (odd p: code = 2a + quad; p = 2: code = u + 2k, chi_2 = chi_{-4}^u psi, psi of order 2^k).
# V(N, types) >= CM(N, chi) for every even chi mod N with these local types (Section 6 of the paper).
from math import floor, sqrt, log, pi
import cypari2
pari = cypari2.Pari()
_h = {}
def hcl(D):
    if D not in _h: _h[D] = int(pari.qfbclassno(-D))
    return _h[D]

def phi_pk(p, k): return 1 if k == 0 else p**k - p**(k - 1)

def locinfo(p, tau, x, y):          # tau 0 split, 1 inert, 2 ram
    if tau == 0: v = max(x, y); ne = x + y; r = phi_pk(p, x) * phi_pk(p, y) / phi_pk(p, v)
    elif tau == 1: v = x; ne = 2 * x; r = 1.0 if x == 0 else (p**(2*x) - p**(2*x - 2)) / phi_pk(p, x)
    else: v = (x + 1) // 2; ne = x; r = phi_pk(p, x) / phi_pk(p, v)
    return r, v, ne

def Eloc(p, tau, x, y, s, only2):
    tot = 0.0
    for dx in (0, 1):
        for dy in ((0, 1) if tau == 0 else (0,)):
            mx, my = x - dx, y - dy
            if mx < 0 or my < 0: continue
            r, v, ne = locinfo(p, tau, mx, my)
            if only2 and not (ne <= 1 or (tau == 0 and mx == 1 and my == 1) or (tau == 1 and mx == 1) or (tau == 2 and mx == 2)):
                continue
            if s > v: continue
            tot += (-1) ** (dx + dy) * r
    return tot

def Ehat(p, e, delta, tau, s, corr):
    tot = 0.0
    for x in range(0, 2 * e + 1):
        for y in (range(0, e + 1) if tau == 0 else (0,)):
            r, v, ne = locinfo(p, tau, x, y)
            if delta + ne > e: continue
            w = e - delta - ne + 1
            if not corr: Ev = Eloc(p, tau, x, y, s, False)
            elif p == 2: Ev = Eloc(p, tau, x, y, s, True)
            else:
                sq = (x <= 1 and y <= 1) if tau == 0 else (x <= 1)
                Ev = ((-1) ** (x + y)) if (sq and s <= 0) else 0.0
            tot += Ev * w
    return tot

def s_odd(delta, code):
    a, quad = code // 2, code & 1
    if delta == 0: return a
    if a >= 2: return a
    if a == 1: return 0 if quad else 1
    return 1

def cond2(u, k): return 2 * u if k == 0 else k + 2
def s_two(delta, code, up):
    u, k = code & 1, code >> 1
    if delta == 0: return cond2(u, k)
    if delta == 2: return 2 * ((u + 1) % 2) if k == 0 else k + 2
    if k == 0: return 3
    if k == 1: return 2 * ((u + up) % 2)
    return k + 2

def fund_divisors(fac):
    """fundamental discriminants -D with D | N; fac = {p: e}"""
    odd = [p for p in fac if p != 2]; e2 = fac.get(2, 0); out = []
    for mask in range(1 << len(odd)):
        d = 1; n3 = 0
        for i, p in enumerate(odd):
            if mask >> i & 1: d *= p; n3 += (p % 4 == 3)
        cands = []
        if d % 4 == 3: cands.append(d)
        if e2 >= 2 and d % 4 == 1: cands.append(4 * d)
        if e2 >= 3: cands.append(8 * d)
        for D in cands:
            if D >= 3: out.append((D, 1 if n3 % 2 == 0 else 0))
    return out

def V(fac, types):
    tot = 0.0
    for D, up in fund_divisors(fac):
        prE = 1.0; prC = 1.0
        for p, e in fac.items():
            if p == 2:
                delta = 3 if D % 8 == 0 else (2 if D % 4 == 0 else 0)
                tau = 2 if delta else (0 if (-D) % 8 == 1 else 1)
                s = s_two(delta, types[p], up)
            else:
                delta = 1 if D % p == 0 else 0
                tau = 2 if delta else (0 if int(pari.kronecker(-D, p)) == 1 else 1)
                s = s_odd(delta, types[p])
            prE *= Ehat(p, e, delta, tau, s, False)
            prC *= Ehat(p, e, delta, tau, s, True) if D > 4 else 0.0
        tot += hcl(D) * (prE - prC)
    return tot
