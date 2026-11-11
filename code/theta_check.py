# Numerical checks of the constants in the proof of the lemma bounding theta(q, e) (Section 6 of the paper).
# Run: python3 code/theta_check.py
from mpmath import mp, mpf, sqrt, log, diff
mp.dps = 40
def theta1(q): q = mpf(q); return (2+sqrt(q))/(q+1)*(1+log(q)/2)
def R(q,e): q = mpf(q); return mpf(11)/10*e*((e+1)**2+2*e)*q**(-mpf(e-1)/2)
print("theta(11,1) =", theta1(11))
# theta(q,1) decreasing for q >= 11: derivative negative (check sign analytically below, numerically on a grid)
import sympy as sp
x = sp.symbols('x', positive=True)   # x = sqrt(q)
f = (2+x)*(1+sp.log(x))/(x**2+1)    # theta(q,1) with log q = 2 log x
num = sp.simplify(sp.diff(f,x)*(x**2+1)**2)
print("numerator of d/dx:", sp.expand(num))
xs = [sp.sqrt(11)+k*0.01 for k in range(0,20000)]
print("derivative < 0 on grid x in [sqrt 11, sqrt11+200]:", all(float(num.subs(x,v))<0 for v in xs[::50]))
# (2+sqrt q)/(q+1) >= q^{-1/2}
print("(2+sqrt q)/(q+1) >= 1/sqrt q  <=> 2 sqrt q >= 1: trivially")
# (1+e y)/(1+y) <= e for e>=1, y>=0: trivially
# P(e+1)/P(e) with P(e)=e(e^2+4e+1)
P = lambda e: e*(e*e+4*e+1)
print("P(e+1)/P(e), e=2..6:", [sp.Rational(P(e+1),P(e)) for e in range(2,7)], "33/13 =", 33/13)
print("33/13 / sqrt(11) =", 33/13/float(sp.sqrt(11)))
print("R(2000,2) =", R(2000,2), " R(11,26) =", R(11,26))
# Compare R with the true bound U(q,e)(1+e log q/2)/theta(q,1) at a few points, and sanity: R >= that
def U(q,e): q=mpf(q); return ((e+1)**2+2*e)/(1-1/q)*q**(-mpf(e)/2)
for (q,e) in [(2000,2),(11,26),(2003,3),(5000,2),(11,30)]:
    print(q,e, "U*(1+e log q/2)/theta1 =", U(q,e)*(1+e*log(q)/2)/theta1(q), "<= R =", R(q,e))
