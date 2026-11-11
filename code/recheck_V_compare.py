# Compare V (code/cmV.py), split by CM discriminant D, with the Magma CM proxies of logs/V/N*.log.
# Expect V_D >= proxy_D; V_D = proxy_D for D > 4 when the proxy is exact.
# Run from the repository root: python3 code/recheck_V_compare.py
import sys, re, glob
sys.path.insert(0, 'code')
import cmV
from sympy import factorint
def V_byD(fac, types):
    out = {}
    for D, up in cmV.fund_divisors(fac):
        prE = 1.0; prC = 1.0
        for p, e in fac.items():
            if p == 2:
                delta = 3 if D % 8 == 0 else (2 if D % 4 == 0 else 0)
                tau = 2 if delta else (0 if (-D) % 8 == 1 else 1)
                s = cmV.s_two(delta, types[p], up)
            else:
                delta = 1 if D % p == 0 else 0
                tau = 2 if delta else (0 if int(cmV.pari.kronecker(-D, p)) == 1 else 1)
                s = cmV.s_odd(delta, types[p])
            prE *= cmV.Ehat(p, e, delta, tau, s, False)
            prC *= cmV.Ehat(p, e, delta, tau, s, True) if D > 4 else 0.0
        out[D] = cmV.hcl(D) * (prE - prC)
    return out
rows = 0; bad = []; eq = 0; gt = 0
for f in sorted(glob.glob('logs/V/N*.log')):
    txt = open(f).read()
    m = re.search(r'LEVEL N=(\d+) discs=\[ ?([^\]]*?) ?\]', txt)
    N = int(m.group(1)); discs = [int(x) for x in m.group(2).replace(' ', '').split(',') if x]
    fac = dict(sorted(factorint(N).items())); ps = list(fac)
    for line in txt.splitlines():
        m = re.match(r'VAL N=\d+ type=\[ ?([^\]]*?) ?\] order=(\d+) cond=(\d+) dim=(\d+) cm=\[ ?([^\]]*?) ?\]', line)
        if not m: continue
        t = [int(x) for x in m.group(1).replace(' ', '').split(',')]
        cm = [int(x) for x in m.group(5).replace(' ', '').split(',') if x]
        v = V_byD(fac, dict(zip(ps, t)))
        rows += 1
        for D, c in zip(discs, cm):
            vd = round(v.get(D, 0.0), 6)
            if vd + 1e-9 < c: bad.append((N, t, m.group(3), D, vd, c))
            elif abs(vd - c) < 1e-9: eq += 1
            else: gt += 1
        tot = sum(v.values())
        print("N=%d type=%s ord=%s cond=%s dim=%s proxy=%s V_D=%s" % (N, t, m.group(2), m.group(3), m.group(4), cm,
              [round(v.get(D, 0), 2) for D in discs]))
print("ROWS", rows, "D-entries equal:", eq, " V_D > proxy:", gt, " VIOLATIONS (V_D < proxy):", bad)
