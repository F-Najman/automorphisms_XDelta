"""Interval certificates for the finite inequalities in the CM DFS.

Local factors are recomputed with exact integers; roots/logs use 128-bit
intervals. The infinite tails are treated analytically in Section 6 of the paper.
Run from the repository root: sage code/review_dfs_bounds.sage
"""
I=RealIntervalField(128)

@cached_function
def L(p,e,delta,tau,s):
    ans=ZZ(0)
    for i in range(e//2+1):
        rho=1 if i==0 else (p**i+({'split':-1,'inert':1,'ram':0}[tau])*p**(i-1))
        for j in range(e-delta-2*i+1):
            if i+j < s: continue
            mult=1 if j==0 else (2 if tau=='split' else (1 if tau=='ram' and j==1 else 0))
            ans+=rho*mult*(e-delta-2*i-j+1)
    return ans

def L0(p,e,s):
    return max(L(p,e,0,'split',s),L(p,e,0,'inert',s))

def max_interval(vals):
    return I(max(v.lower() for v in vals),max(v.upper() for v in vals))

def lam(p,e):
    if p!=2:
        # All other odd-prime character types are dominated by these two.
        vals=[I(L0(p,e,a))+I(p).sqrt()*L(p,e,1,'ram',b) for a,b in [(0,1),(1,0)]]
    else:
        types=[(0,2,3)]
        if e>=2: types.append((2,0,3))
        if e>=3: types.append((3,3,0))
        vals=[I(L0(p,e,a))+2*L(p,e,2,'ram',b)+I(8).sqrt()*L(p,e,3,'ram',c)
              for a,b,c in types]
    return max_interval(vals)

def f(p,e):
    return lam(p,e)/(p**e+p**(e-1))*(1+e*I(p).log()/2)

wa=I(0); wb=I(0)
for p in prime_range(11,2000):
    base=f(p,1)
    for e in range(1,26):
        v=f(p,e)
        assert v.upper()<1
        wa=max_interval([wa,v])
        if e>=2:
            ratio=v/base
            assert ratio.upper()<1
            wb=max_interval([wb,ratio])
print("Finite range max f upper:",wa.upper())
print("Finite range max f(p,e)/f(p,1), e>=2:",wb.upper())

def U(p,e):
    return I((e+1)**2+2*e)*I(p)**(-I(e)/2)/(1-I(1)/p)

def tail_a(p,e):
    return U(p,e)*(1+e*I(p).log()/2)

def tail_b(p,e):
    return tail_a(p,e)/((2+I(p).sqrt())/(p+1)*(1+I(p).log()/2))

for p,e in [(11,26),(2000,2)]:
    assert tail_a(p,e).upper()<1 and tail_b(p,e).upper()<1
    print("Tail endpoint",p,e,"a <=",tail_a(p,e).upper(),"b <=",tail_b(p,e).upper())

S=I(1)
for p in [2,3,5,7]:
    v=max_interval([I(1)]+[f(p,e) for e in range(1,60)])
    S*=v
    print("Small-prime factor",p,"<=",v.upper())
assert S.upper()<16
e=60
U2=I((e+1)**2+4*e)*I(2)**(-e/2)/(1-I(1)/2)
ftail=U2*(1+e*I(2).log()/2)
Gtail=(48/I.pi())*16*ftail+7*I(2)**(-59)+6*I(2)**(-30)
assert Gtail.upper()<1
print("Uniform majorant when some exponent >=60:",Gtail.upper())

power_tail=(48/I.pi())*16*I(3718)*I(11)/10*I(10)**(-20)*(1+20*I(10).log())+14*I(10)**(-40)+12*I(10)**(-20)
assert power_tail.upper()<1
print("Majorant for discarded q^e > 10^40, e<60:",power_tail.upper())
print("ALL CHECKS PASSED")
