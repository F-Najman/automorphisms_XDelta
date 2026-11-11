# Check of the index formula of the theorem on the normalizer of Gamma_Delta(N) (Section 10 of the paper)
# against the sigma-H enumeration of code/aut_orders.py (column AutMod of logs/aut_orders_*.tsv).
import csv, sys
from math import gcd
from fractions import Fraction
from sympy import factorint

def psi(n):
    r = Fraction(n)
    for p in factorint(n): r *= Fraction(p + 1, p)
    return r

def group(N, gens):
    G = {1 % N}; frontier = [1 % N]
    gens = [(-1) % N] + gens
    while frontier:
        x = frontier.pop()
        for g in gens:
            y = x * g % N
            if y not in G: G.add(y); frontier.append(y)
    return G

def formula(N, gens):
    units = [d for d in range(N) if gcd(d, N) == 1]
    D = group(N, gens)
    i_D = len(units) // len(D)
    good = []
    for h in range(1, N + 1):
        if N % (h * h): continue
        if any((d * d - 1) % h for d in D): continue
        if any(d not in D for d in units if (d - 1) % (N // h) == 0): continue
        good.append(h)
    hs = max(good)
    assert all(hs % h == 0 for h in good), (N, gens, good)
    ns = N // hs**2
    M = N // hs
    Dbar = {d % M for d in D}
    fm = factorint(M)
    primes_n = list(factorint(ns))
    count = 0
    for mask in range(1 << len(primes_n)):
        Qp = [primes_n[i] for i in range(len(primes_n)) if mask >> i & 1]
        def iota(d):
            # CRT: invert the p-components for p in Qp, keep the others
            x = 0
            for p, e in fm.items():
                q = p**e; c = pow(d, -1, q) if p in Qp else d % q
                Mq = M // q
                x += c * Mq * pow(Mq, -1, q)
            return x % M
        if {iota(d) for d in Dbar} == Dbar: count += 1
    val = i_D * psi(N) / psi(ns) * count
    assert val.denominator == 1
    return int(val), hs, ns

bad = n = 0
for path in sys.argv[1:]:
    for row in csv.DictReader(open(path), delimiter='\t'):
        N = int(row['N']); gens = [int(x) for x in row['gens'].split(',') if x]
        val, hs, ns = formula(N, gens)
        n += 1
        if val != int(row['AutMod']):
            bad += 1; print('MISMATCH', N, gens, val, row['AutMod'], hs, ns)
print(n, 'curves,', bad, 'mismatches')
