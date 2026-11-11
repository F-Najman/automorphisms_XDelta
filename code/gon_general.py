# (C2) for general N > 200: (7/800) psi(N) idx > 2(l+1), l = least prime not dividing N.
# Reports N in (200, 10^7] where (C2) fails for idx = 2 (N > 500) or idx = 4 (200 < N <= 500).
from sympy import primerange
M = 10**7
spf = list(range(M + 1))
for p in primerange(2, int(M**0.5) + 1):
    for j in range(p * p, M + 1, p):
        if spf[j] == j: spf[j] = p
bad2, bad4 = [], []
for N in range(201, M + 1):
    n = N; psi = N; ps = []
    while n > 1:
        p = spf[n]; ps.append(p)
        while n % p == 0: n //= p
        psi = psi // p * (p + 1)
    l = 2
    while N % l == 0: l += 1 if l == 2 else 2          # l runs over 2,3,5,7,9?.. -> fix below
    while any(l % q == 0 for q in range(2, int(l**0.5) + 1)) or N % l == 0:
        l += 1
    for idx, lst in ((2, bad2), (4, bad4)):
        if not (7 * psi * idx > 800 * 2 * (l + 1)):
            if (idx == 2 and N > 500) or (idx == 4): lst.append(N)
print("idx=2 failures with N>500:", bad2[:20], len(bad2))
print("idx=4 failures with 200<N:", bad4[:20], len(bad4))
print("DONE")
