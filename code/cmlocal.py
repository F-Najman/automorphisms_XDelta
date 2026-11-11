# Local factors for the CM bound (Section 6 of the paper), general N.
# For a prime p with p^e || N, a CM field F with delta = v_p(D) and splitting type tau of p in F, and the local type
# t of chi_p, L(p, e, delta, tau, t) = sum over admissible (i, j) (i = v_p(c), j = v_p(N g)) with
#   delta + 2i + j <= e,  i + j >= s(delta, t, eps)   (compatibility: cond(chi eps_F) | c N(g))
# of rho(p^i) a(p^j) (e - delta - 2i - j + 1).
# Odd p: types t = (a, quad) with a = conductor exponent of chi_p, quad = [a == 1 and chi_p is the Legendre symbol].
# p = 2: types t = (u, k), chi_2 = chi_{-4}^u psi with psi of order 2^k on <5>; eps_{F,2} in {1, chi_-4, chi_8, chi_-8}.
from math import sqrt

def rho(p, i, tau):
    if i == 0: return 1.0
    return p**i * {'split': 1 - 1/p, 'inert': 1 + 1/p, 'ram': 1.0}[tau]

def amult(j, tau):
    if j == 0: return 1
    if tau == 'split': return 2
    if tau == 'ram': return 1 if j == 1 else 0
    return 0                      # inert: no primitive part

def L(p, e, delta, tau, s):
    tot = 0.0
    for i in range(0, e // 2 + 1):
        for j in range(0, e - delta - 2 * i + 1):
            if delta + 2 * i + j > e or i + j < s: continue
            am = amult(j, tau)
            if am == 0: continue
            tot += rho(p, i, tau) * am * (e - delta - 2 * i - j + 1)
    return tot

def types_odd(p, e):
    ts = [(0, False)]
    if e >= 1:
        ts.append((1, True))
        if p > 3: ts.append((1, False))
    ts += [(a, False) for a in range(2, e + 1)]
    return ts

def s_odd(delta, t):
    a, quad = t
    if delta == 0: return a
    if a >= 2: return a
    if a == 1: return 0 if quad else 1
    return 1

def types_two(e):
    ts = [(0, 0)]
    if e >= 2: ts.append((1, 0))
    for k in range(1, e - 1):
        ts += [(0, k), (1, k)]
    return ts

def cond2(u, k):
    return 2 * u if k == 0 else k + 2

def s_two(delta, t, up=0):
    # delta 0: eps = 1; delta 2: eps = chi_-4; delta 3: eps = chi_8 (up=0) or chi_-8 (up=1)
    u, k = t
    if delta == 0: return cond2(u, k)
    if delta == 2: return cond2((u + 1) % 2, k) if k == 0 else k + 2
    if k == 0: return 3
    if k == 1: return 2 * ((u + up) % 2)
    return k + 2

def Lam(p, e):
    """max over local types of sum_delta p^(delta/2) * (worst splitting type) L"""
    best = 0.0
    if p == 2:
        for t in types_two(e):
            v = max(L(2, e, 0, 'split', s_two(0, t)), L(2, e, 0, 'inert', s_two(0, t)))
            if e >= 2: v += 2 * L(2, e, 2, 'ram', s_two(2, t))
            if e >= 3: v += sqrt(8) * max(L(2, e, 3, 'ram', s_two(3, t, up)) for up in (0, 1))
            best = max(best, v)
    else:
        for t in types_odd(p, e):
            v = max(L(p, e, 0, 'split', s_odd(0, t)), L(p, e, 0, 'inert', s_odd(0, t)))
            v += sqrt(p) * L(p, e, 1, 'ram', s_odd(1, t))
            best = max(best, v)
    return best

def psi_pe(p, e): return p**e + p**(e - 1)

if __name__ == "__main__":
    for p in [3, 5, 7, 11]:
        print(p, "Lam(p,1) =", round(Lam(p, 1), 3), "expected max(4, 2+sqrt p) =", round(max(4, 2 + sqrt(p)), 3),
              "| Lam(p,2) =", round(Lam(p, 2), 3), "expected p+5+3sqrt p =", round(p + 5 + 3 * sqrt(p), 3))
    for e in range(1, 9):
        print("p=2 e=%d Lam=%.2f ratio=%.3f" % (e, Lam(2, e), Lam(2, e) / psi_pe(2, e)))
    for e in range(1, 7):
        print("p=3 e=%d ratio=%.3f  p=5 ratio=%.3f  p=7 ratio=%.3f" % (e, Lam(3, e)/psi_pe(3, e), Lam(5, e)/psi_pe(5, e), Lam(7, e)/psi_pe(7, e)))
