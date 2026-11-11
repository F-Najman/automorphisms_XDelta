load "lib.m";
// ---------- X_1(13): genus 2 ----------
N := 13; Delta := SubgroupList(N, [1]);
M, B := QBasis(N, Delta, 200);
print "N=13 dim", #B;
PS<q> := PowerSeriesRing(Rationals());
f1 := PS!B[1]; f2 := PS!B[2];
x := f1/f2; // careful with valuations
dx := Derivative(x);
y := q*dx/f2;
// find relation y^2 = sum c_i x^i, i<=6
prec := 150;
mons := [ x^i : i in [0..6] ] cat [ -y^2 ];
// x has a pole of order? compute Laurent series
L<t> := LaurentSeriesRing(Rationals());
xL := L!f1/L!f2; yL := L!q*Derivative(xL)/(L!f2);
vals := [ xL^i : i in [0..6] ] cat [ yL^2 ];
v0 := Min([Valuation(v) : v in vals]);
mat := Matrix(Rationals(), [ [Coefficient(v, i) : i in [v0..v0+80]] : v in vals ]);
K := Kernel(mat); print "kernel dim", Dimension(K);
b := Basis(K)[1]; b := b/b[8];
Px<X> := PolynomialRing(Rationals());
f := -&+[ b[i+1]*X^i : i in [0..6] ];
print "y^2 =", f;
C := HyperellipticCurve(f);
print "Igusa-Clebsch", IgusaClebschInvariants(C);
G := GeometricAutomorphismGroup(C);
print "Geometric Aut X_1(13): order", #G, GroupName(G);
A := AutomorphismGroup(C); print "Aut over Q order", #A;
quit;
