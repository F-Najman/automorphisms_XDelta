# Numerical test on X(1) of the behaviour of T_l^n(P) at elliptic and CM points (Section 5 of the paper).
# a_n(P) = #supp T_l^n(P) = #{ j(E/C) : C cyclic of order l^(n-2k), 0 <= k <= n/2 },  E = C/(Z + tau Z).
# Lemma E claims a_n / sum_k psi(l^(n-2k)) -> 1/e(P).  On X(1) the elliptic points are j = 1728 (e=2), j = 0 (e=3).
# Test points: tau = i (e=2), tau = (1+sqrt-7)/2 (CM, NOT elliptic; 2 splits in Q(sqrt-7), class number 1),
#              tau = (1+sqrt-23)/2 (2 splits, h = 3), tau = sqrt(-2) (2 ramified), tau = 0.1 + 1.3i (generic).
# Cyclic sublattices of index l^m of Z + tau Z: Hermite forms (a b; 0 d), ad = l^m, 0 <= b < d, gcd(a,b,d) = 1;
# the quotient curve is C/(Z + (a tau + b)/d Z) up to homothety.
import mpmath as mp
mp.mp.dps = 60

def reduce_tau(t):
    for _ in range(10000):
        t = t - mp.nint(mp.re(t))
        if abs(t) < 1 - mp.mpf(10)**-40: t = -1 / t
        else: break
    return t

def jinv(t):
    t = reduce_tau(t)
    q = mp.exp(2j * mp.pi * t)
    # j = E4^3 / Delta via Eisenstein series (q small after reduction: |q| <= e^{-pi sqrt 3})
    E4 = 1 + 240 * mp.nsum(lambda n: n**3 * q**n / (1 - q**n), [1, mp.inf])
    E6 = 1 - 504 * mp.nsum(lambda n: n**5 * q**n / (1 - q**n), [1, mp.inf])
    return 1728 * E4**3 / (E4**3 - E6**2)

def cyclic_quotients(tau, l, m):
    out = []
    for i in range(m + 1):
        a, d = l**i, l**(m - i)
        for b in range(d):
            if mp.mpf(0) == 0 and __import__('math').gcd(__import__('math').gcd(a, b), d) != 1: continue
            out.append((a * tau + b) / d)
    return out

def canon(t):
    # canonical representative of the SL_2(Z)-orbit of t in the standard fundamental domain
    eps = mp.mpf(10)**-40
    for _ in range(100000):
        t = t - mp.nint(mp.re(t))
        if abs(t) < 1 - eps: t = -1 / t
        else: break
    if mp.re(t) < -mp.mpf(1)/2 + eps: t = t + 1
    if abs(abs(t) - 1) < eps and mp.re(t) < 0: t = -mp.conj(t)
    return (mp.nstr(mp.re(t), 30), mp.nstr(mp.im(t), 30))

def a_n(tau, l, n):
    keys = [canon(t) for k in range(n // 2 + 1) for t in cyclic_quotients(tau, l, n - 2 * k)]
    return len(set(keys)), len(keys)

def psi(l, m): return 1 if m == 0 else l**m + l**(m - 1)

pts = {'i (elliptic, e=2)': mp.mpc(0, 1),
       '(1+sqrt-7)/2 (CM, e=1)': (1 + mp.sqrt(-7)) / 2,
       '(1+sqrt-23)/2 (CM, e=1)': (1 + mp.sqrt(-23)) / 2,
       'sqrt-2 (CM, e=1)': mp.sqrt(-2),
       'generic': mp.mpc(0.1, 1.3)}
for l in [2, 3]:
    for name, tau in pts.items():
        row = []
        for n in range(2, 11 if l == 2 else 7):
            a, tot = a_n(tau, l, n)
            S = sum(psi(l, n - 2 * k) for k in range(n // 2 + 1))
            assert tot == S
            row.append("%d:%d/%d=%.3f" % (n, a, S, a / S))
        print("l=%d %-26s %s" % (l, name, "  ".join(row)), flush=True)
print("DONE")

# ---- the repaired invariant: multiplicities in the divisor T_l^n(P) ----
# T_l^n = sum_m c[n][m] T_{l^m} (cyclic parts; on X(1) <l> = 1): T_l T_{l^m} = T_{l^{m+1}} + l T_{l^{m-1}} (m >= 2),
# T_l T_l = T_{l^2} + (l+1) T_1.  Claim: for elliptic P of order e, at most 2(n/2+1) points of T_l^n(P) have
# multiplicity not divisible by e; for non-elliptic P (or order != e) there are >= (l-1) l^(n-1)/3 of them.
from collections import Counter
def coeffs(l, n):
    c = {0: 1}
    for _ in range(n):
        d = Counter()
        for m, v in c.items():
            if m == 0: d[1] += v
            elif m == 1: d[2] += v; d[0] += (l + 1) * v
            else: d[m + 1] += v; d[m - 1] += l * v
        c = dict(d)
    return c
def mult_profile(tau, l, n, e):
    cnt = Counter()
    for m, v in coeffs(l, n).items():
        for t in (cyclic_quotients(tau, l, m) if m > 0 else [tau]):
            cnt[canon(t)] += v
    return sum(1 for x in cnt.values() if x % e), len(cnt)
for l in [2, 3]:
    for name, tau in pts.items():
        row = []
        for n in range(3, 10 if l == 2 else 7):
            row.append("%d:%s/%s" % ((n,) + mult_profile(tau, l, n, 2)))
        print("MULT l=%d e=2 test %-26s #points with odd multiplicity / #points: %s" % (l, name, "  ".join(row)), flush=True)
print("DONE2")
