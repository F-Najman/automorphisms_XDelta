// Exhaustive Frobenius root-of-unity checks, with a proved order cutoff.
// Run from code/: magma -b review_frobenius_fields.m
load "lib.m";
cases := [ <25,[1]>, <63,[37,59]>,
           <63,[a : a in [1..62] | GCD(a,63) eq 1]> ];
for c in cases do
  N := c[1]; Delta := SubgroupList(N,c[2]);
  ell := N eq 25 select 101 else 127;
  chars := CharsTrivialOn(N,Delta);
  MS := CuspidalSubspace(ModularSymbols(chars,2,1));
  T := HeckeOperator(MS,ell); h := CharacteristicPolynomial(T);
  Q<x> := PolynomialRing(Rationals());
  Qy<y> := PolynomialRing(Q);
  PF := Q!Resultant(Qy![Q!a : a in Coefficients(h)],x^2-y*x+ell);
  Pred := &*[fa[1] : fa in Factorization(PF)];
  Rat := Q!Resultant(Qy![Q!a : a in Coefficients(Pred)],Evaluate(Pred,x*y));
  d := Degree(Rat);
  // phi(n)^2 >= n/2, hence phi(n) <= d implies n <= 2*d^2.
  ns := [n : n in [2..2*d^2] | EulerPhi(n) le d];
  orders := [n : n in ns | Rat mod (Q!CyclotomicPolynomial(n)) eq 0];
  kk := #orders eq 0 select 1 else LCM(orders);
  printf "N=%o |Delta|=%o g=%o ratio_degree=%o max_candidate_order=%o cyclotomic_orders=%o field_degree=%o\n",
         N,#Delta,Dimension(MS),d,Max(ns),orders,kk;
  assert kk eq (N eq 25 select 1 else 2);
  if N eq 63 and #Delta eq 36 then
    M := CuspForms(chars,2); prec := 150; F := GF(127^2);
    PS<q> := PowerSeriesRing(F,prec);
    Q2,R := Relations([PS!b : b in qExpansionBasis(M,prec)],2,F);
    C := Curve(ProjectiveSpace(R),Q2);
    assert Dimension(C) eq 1 and Degree(C) eq 8;
    A := AutomorphismGroup(C);
    assert #A eq 48;
    print "X_0(63): full geometric special-fibre automorphism group has order 48.";
  end if;
end for;
print "ALL CHECKS PASSED";
quit;
