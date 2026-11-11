// Section 9 of the paper. Run from code/: magma -b review_xdelta63.m
// Petri certifies the model: full quadratic relations, 1-dimensional intersection, degree 16.
// X_{Delta2}(63): is the deck involution (diamond of order 2) normal in Aut, with Aut/<delta> of order 48 = |Aut X_0(63)|?
load "lib.m";
N := 63; Delta := SubgroupList(N, [37, 59]);
chars := CharsTrivialOn(N, Delta);
M := CuspForms(chars, 2); g := Dimension(M);
ell := 127; k := 2; Fk := GF(ell^k);
prec := 250;
B := qExpansionBasis(M, prec);
PSF<q> := PowerSeriesRing(Fk, prec);
Bf := [ PSF!b : b in B ];
Q2, R := Relations(Bf, 2, Fk);
C := Curve(ProjectiveSpace(R), Q2);
assert g eq 9 and #Q2 eq 21;
assert Dimension(C) eq 1 and Degree(C) eq 16;
print "Canonical quadratic model: dimension 1, degree 16, 21 quadrics.";
printf "X_Delta2(63): genus %o, quadrics %o\n", g, #Q2;
// deck involution: acts as +1 on the trivial-character part S_2(Gamma_0(63)) and -1 on the chi-part.
S0 := CuspForms(Gamma0(N), 2);
B0 := [ PSF!qExpansion(f, prec) : f in Basis(S0) ];
Mq := Matrix(Fk, [ [Coefficient(b, n) : n in [1..prec-1]] : b in Bf ]);
coords := function(f) ok, s := IsConsistent(Mq, Vector(Fk, [Coefficient(f, n) : n in [1..prec-1]])); assert ok; return s; end function;
V0 := sub< VectorSpace(Fk, g) | [ coords(b) : b in B0 ] >;
printf "dim of Gamma_0 part: %o\n", Dimension(V0);
// the chi-part: forms f with f|<d> = -f; use complement given by the chi-character space
chi := [c : c in chars | Order(c) eq 2][1];
S1 := CuspForms(chi, 2);
B1 := [ PSF!qExpansion(f, prec) : f in Basis(S1) ];
V1 := sub< VectorSpace(Fk, g) | [ coords(b) : b in B1 ] >;
printf "dim of chi part: %o\n", Dimension(V1);
T := Matrix(Fk, [ Eltseq(v) : v in Basis(V0) ] cat [ Eltseq(v) : v in Basis(V1) ]);  // rows = new basis in old coords
Dg := DiagonalMatrix(Fk, [1 : i in [1..Dimension(V0)]] cat [-1 : i in [1..Dimension(V1)]]);
Adelta := T^-1 * Dg * T;   // row convention: row i = delta^*(f_i) in basis f
act := func< A, v | Vector(Fk, Eltseq(A * Transpose(Matrix(Fk, [Eltseq(v)])))) >;
AmbR := CoordinateRing(Ambient(C));
// check delta preserves the curve
imgs := [ &+[ Adelta[i][j]*AmbR.j : j in [1..g] ] : i in [1..g] ];
printf "delta preserves the canonical ideal: %o\n", &and[ Evaluate(f, imgs) in Ideal(C) : f in Q2 ];
delta := iso< C -> C | imgs, [ &+[ (Adelta^-1)[i][j]*AmbR.j : j in [1..g] ] : i in [1..g] ] >;
auts := Automorphisms(C);
printf "|Aut X_Delta2(63)| = %o\n", #auts;
// Coerce the explicit diamond into Magma's exact curve automorphism group.
A := AutomorphismGroup(C);
dg := A!delta;
assert Order(dg) eq 2 and #A eq 96;
assert &and[a*dg eq dg*a : a in Generators(A)];
printf "Exact group operations: diamond central in group of order %o.\n",#A;
print "ALL CHECKS PASSED";
quit;
