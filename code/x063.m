// Aut X_0(63) acting on its 8 cusps; which elements preserve the branch cusps (denominators 3, 21) of X_1(63) -> X_0(63)?
load "lib.m";
N := 63; ell := 127; F := GF(ell); k := 2; Fk := GF(ell^k);
S := CuspForms(Gamma0(N), 2); g := Dimension(S);
prec := 200;
B := [ qExpansion(f, prec) : f in Basis(S) ];
PSF<q> := PowerSeriesRing(Fk, prec);
Bf := [ PSF!b : b in B ];
Q2, R := Relations(Bf, 2, Fk);
C := Curve(ProjectiveSpace(R), Q2);
printf "X_0(63): genus %o, %o quadrics, genus of model %o\n", g, #Q2, Genus(C);
// matrix of a linear operator given by its action on q-expansions, in the basis Bf
function OpMatrix(images)
  // images[i] = q-expansion of op(f_i); express in basis Bf
  M := Matrix(Fk, [ [Coefficient(b, n) : n in [1..prec-1]] : b in Bf ]);
  rows := [];
  for im in images do
    v := Vector(Fk, [Coefficient(im, n) : n in [1..prec-1]]);
    ok, sol := IsConsistent(M, v); assert ok;
    Append(~rows, Eltseq(sol));
  end for;
  return Matrix(Fk, rows);  // row i = coefficients of op(f_i)
end function;
// Atkin-Lehner matrices (rational)
W := AssociativeArray();
for Q in [7, 9, 63] do
  A := AtkinLehnerOperator(S, Q);
  W[Q] := ChangeRing(A, Fk);
end for;
// S_3: f(z) -> f(z + 1/3): a_n -> zeta3^n a_n
z3 := RootOfUnity(3, Fk);
S3 := OpMatrix([ &+[ Coefficient(b, n)*z3^n*q^n : n in [1..prec-1] ] + O(q^prec) : b in Bf ]);
S3i := S3^-1;
e1 := Vector(Fk, [1] cat [0 : i in [2..g]]);
pt := func< v | C!Eltseq(v) >;
// phi(sigma x) = A^T-type action: with row convention (row i = sigma^* f_i), phi(sigma x) = A * phi(x) as column vector
act := func< A, v | Vector(Fk, Eltseq(A * Transpose(Matrix(Fk, [Eltseq(v)])))) >;
Pinf := e1;
P0 := act(W[63], Pinf); P17 := act(W[9], Pinf); P19 := act(W[7], Pinf);
cusps_good := [ Pinf, P0, P17, P19 ];            // denominators 1, 63(=0), 7, 9
branch := [ act(S3, P0), act(S3i, P0), act(S3, P17), act(S3i, P17) ];   // 1/3, 2/3, 10/21, -4/21
allc := cusps_good cat branch;
cpts := [ pt(v) : v in allc ];
printf "cusp points distinct: %o\n", #SequenceToSet(cpts) eq 8;
A := AutomorphismGroup(C);
auts := Automorphisms(C);
printf "|Aut X_0(63) over F_(127^2)| = %o\n", #auts;
bset := SequenceToSet(cpts[5..8]);
cset := SequenceToSet(cpts);
npres_cusps := 0; npres_branch := 0;
for a in auts do
  im := [ a(P) : P in cpts ];
  pc := SequenceToSet(im) eq cset;
  pb := SequenceToSet(im[5..8]) eq bset;
  npres_cusps +:= pc select 1 else 0;
  npres_branch +:= pb select 1 else 0;
end for;
printf "automorphisms preserving the set of cusps: %o ; preserving the branch cusps {1/3,2/3,10/21,-4/21}: %o\n", npres_cusps, npres_branch;
quit;
