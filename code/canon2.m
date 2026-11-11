// canon2.m = canon.m plus: genus 3 (plane quartic) and no genus check for g >= 12 (point count still checked).
// Canonical model of X_Delta(N) over F_ell, Frobenius data, automorphism group over F_ell^k.
// usage: magma -b N:=17 D:=1 ell:=103 canon.m     (D = generator(s) of Delta, comma separated; -1 is always added)
load "lib.m";
if not assigned N then N := 17; elif Type(N) eq MonStgElt then N := StringToInteger(N); end if;
if not assigned D then D := "1"; end if;
if not assigned ell then ell := 0; elif Type(ell) eq MonStgElt then ell := StringToInteger(ell); end if;
if not assigned K then Kmax := 0; else Kmax := StringToInteger(K); end if;
gens := [ StringToInteger(s) : s in Split(D, ",") ];
Delta := SubgroupList(N, gens);
mu := MuDelta(N, Delta);
printf "N=%o Delta=%o |Delta|=%o mu=%o\n", N, Delta, #Delta, mu;
t0 := Cputime();
chars := CharsTrivialOn(N, Delta);
M := CuspForms(chars, 2);
g := Dimension(M);
printf "genus g=%o\n", g;
// |B_Delta|: diamonds times number of Q with iota_Q(Delta)=Delta
Zn := Integers(N);
function IotaQ(a, Q) // a in Z/N: inverse mod Q, same mod N/Q
  a := Integers()!a;
  inv := Integers()!((Integers(Q)!a)^-1);
  return Integers()!CRT([inv, a mod (N div Q)], [Q, N div Q]);
end function;
Qs := [ Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1 ];
stableQ := [ Q : Q in Qs | Q eq 1 or Sort([ IotaQ(a,Q) mod N : a in Delta ]) eq Delta ];
Bord := (EulerPhi(N) div #Delta) * #stableQ;
printf "Atkin-Lehner Q with iota_Q(Delta)=Delta: %o ; |B_Delta| = %o\n", stableQ, Bord;
if ell eq 0 then
  ell := NextPrime(100);
  while ell mod N ne 1 do ell := NextPrime(ell); end while;
end if;
printf "ell=%o (ell mod N = %o)\n", ell, ell mod N;
// ---- Frobenius char poly from Hecke operator (ell = 1 mod N so <ell> = 1)
assert ell mod N eq 1;
MS := CuspidalSubspace(ModularSymbols(chars, 2, 1));
assert Dimension(MS) eq g;
T := HeckeOperator(MS, ell);
h := CharacteristicPolynomial(T);
Px<X> := PolynomialRing(Rationals());
PP<Y,Z> := PolynomialRing(Rationals(), 2);
// P(X) = prod (X^2 - t X + ell) over roots t of h  = Res_t(h(t), X^2 - tX + ell)
Pt<tt> := PolynomialRing(Px);
hP := Pt![ Px!c : c in Coefficients(h) ];
PF := Resultant(hP, X^2 - tt*X + ell);
PF := Px!PF;
printf "Tr(T_ell) = %o, predicted #X(F_ell) = %o\n", Trace(T), ell + 1 - Trace(T);
// root-of-unity ratios of Frobenius eigenvalues
Pred := &*[ fa[1] : fa in Factorization(PF) ];
Pxy<xx> := PolynomialRing(Px);
R := Resultant( Pxy![ Px!c : c in Coefficients(Pred) ], Evaluate(Pred, X*xx) ); // roots: ratios
// careful: Evaluate(Pred, X*xx) as polynomial in xx with coeffs in Px
R := Px!R;
cycl := [ n : n in [2..2000] | EulerPhi(n) le Degree(R) and R mod (Px!CyclotomicPolynomial(n)) eq 0 ];
kk := #cycl eq 0 select 1 else LCM(cycl);
printf "root-of-unity orders among Frobenius ratios: %o ; endomorphisms of J over Fbar_ell defined over F_{ell^%o}\n", cycl, kk;
// ---- q-expansions mod ell, canonical model
prec := Ceiling(4*mu/12) + 20;
Bq := qExpansionBasis(M, prec);
F := GF(ell);
PSF<q> := PowerSeriesRing(F, prec);
Bf := [ PSF!b : b in Bq ];
Q2, R2 := Relations(Bf, 2, F);
printf "number of quadrics: %o  (non-hyperelliptic expects %o)\n", #Q2, (g-2)*(g-3) div 2;
Pg := ProjectiveSpace(R2);
if #Q2 eq (g-1)*(g-2) div 2 then
  print "HYPERELLIPTIC (quadrics = rational normal curve)";
  quit;
end if;
Z := Scheme(Pg, Q2);
printf "dim of quadric intersection: %o\n", Dimension(Z);
eqs := Q2;
if Dimension(Z) gt 1 then
  prec3 := Ceiling(6*mu/12) + 20;
  Bq3 := qExpansionBasis(M, prec3);
  PSF3 := PowerSeriesRing(F, prec3);
  Q3, R3 := Relations([PSF3!b : b in Bq3], 3, F);
  Q3 := [ Evaluate(c, [R2.i : i in [1..g]]) : c in Q3 ];
  eqs := Q2 cat Q3;
  printf "added %o cubics\n", #Q3;
end if;
if g eq 3 then
  prec4 := Ceiling(8*mu/12) + 20;
  Bq4 := qExpansionBasis(M, prec4);
  PSF4 := PowerSeriesRing(F, prec4);
  Q4, R4 := Relations([PSF4!b : b in Bq4], 4, F);
  eqs := [ Evaluate(c, [R2.i : i in [1..g]]) : c in Q4 ];
  printf "genus 3: %o quartic relation(s)\n", #eqs;
  assert #eqs eq 1;
end if;
C := Curve(Pg, eqs);
if g lt 12 then
  printf "curve built; checking genus... (%o s)\n", Cputime(t0);
  gC := Genus(C);
  printf "genus of model over F_ell: %o\n", gC;
  assert gC eq g;
end if;
npts := #Points(C);
printf "#C(F_ell) = %o (predicted %o)\n", npts, ell + 1 - Trace(T);
kuse := Max(kk, Kmax);
if kuse gt 1 then
  C := BaseChange(C, GF(ell^kuse));
end if;
printf "computing AutomorphismGroup over GF(%o^%o)... (%o s)\n", ell, kuse, Cputime(t0);
A := AutomorphismGroup(C);
try
  Ap := PermutationGroup(A); ordA := #Ap; gname := GroupName(Ap);
catch e
  printf "PermutationGroup failed (%o); using #A\n", e`Object;
  ordA := #A; gname := "?";
end try;
printf "RESULT N=%o Delta=%o g=%o ell=%o k=%o : |Aut(X_Fbar)| = %o, |B_Delta| = %o, group %o (%o s)\n", N, Delta, g, ell, kuse, ordA, Bord, gname, Cputime(t0);
quit;
