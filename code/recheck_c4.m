// The conditions of the descent lemma (Section 8 of the paper) for X_Delta(N) and every Delta' containing Delta.
// Problem found: coprimality of the Q-charpolys of T_p^2 <p>^-1 for ONE prime p only excludes common NON-CM
// isogeny factors (Ribet). Two CM forms with the same CM field are Qbar-isogenous without being twists of
// Galois conjugates, and T_p^2<p>^-1 separates them when p splits in the CM field.
// Valid certificates, for Delta' ⊋ Delta with g' = g(X_Delta') >= 2:
//   FULL : coprime at a prime p that is inert in every Q(sqrt -D), D | N  (then at most one side has CM forms)
//          => Hom_Qbar(J_Delta', P) = 0, the original (C4).
//   NCM  : coprime at some p (=> no common non-CM factor) AND (C1) for X_Delta': 2 g'_CM < g' - 1.
//          Then Fix((J^H)^nCM) = H by Riemann-Hurwitz, so H is normal (refined lemma (C4')).
// Also prints, per CM discriminant D, the CM-proxy dimension on each side (joint kernel of T_q, q inert, q <= 150).
// Usage: magma -b N:=144 D:=11,71 recheck_c4.m
SetColumns(0);
if not assigned N then N := "144"; end if; N := StringToInteger(N);
if not assigned D then D := "11,71"; end if;
gens := [StringToInteger(s) : s in Split(D, ",")] cat [-1];
G := DirichletGroup(N, CyclotomicField(EulerPhi(N)));
reps := [chi : chi in GaloisConjugacyRepresentatives(G) | &and[Evaluate(chi, a) eq 1 : a in gens]];
discs := [d : d in Divisors(N) | d ge 3 and IsFundamentalDiscriminant(-d)];
ps := [p : p in PrimesUpTo(600) | N mod p ne 0];
pinert := [p : p in ps | &and[KroneckerSymbol(-d, p) eq -1 : d in discs]];
printf "N=%o Delta=%o discs=%o primes inert in all CM fields (<=600): %o\n", N, D, discs, pinert;
Px<x> := PolynomialRing(Rationals());
data := [* *];
for chi in reps do
  chi0 := MinimalBaseRingCharacter(chi);
  M := CuspidalSubspace(ModularSymbols(chi0, 2, 1));
  K := BaseRing(M); deg := (K cmpeq Rationals()) select 1 else Degree(K);
  d := Dimension(M);
  pols := [];
  for p in ps do
    if d eq 0 then Append(~pols, Px!1); continue; end if;
    f := CharacteristicPolynomial(HeckeOperator(M, p)^2 / Evaluate(chi0, p));
    Append(~pols, (K cmpeq Rationals()) select Px!f else Px!Norm(f));
  end for;
  cm := [];
  for dd in discs do
    if d eq 0 then Append(~cm, 0); continue; end if;
    V := VectorSpace(K, d); W := V;
    for q in [q : q in PrimesUpTo(150) | N mod q ne 0 and KroneckerSymbol(-dd, q) eq -1] do
      W := W meet Kernel(HeckeOperator(M, q)); if Dimension(W) eq 0 then break; end if;
    end for;
    Append(~cm, deg * Dimension(W));
  end for;
  Append(~data, <chi0, deg * d, pols, cm>);
  printf "ORBIT cond=%o order=%o size=%o dim=%o cmproxy_by_D=%o\n", Conductor(chi0), Order(chi0), deg, deg*d, cm;
end for;
U, mU := UnitGroup(Integers(N));
S := sub<U | [(Integers(N)!a) @@ mU : a in gens]>;
for rec in Subgroups(U) do
  S2 := rec`subgroup;
  if not (S subset S2) or S2 eq S then continue; end if;
  g2 := [Integers()!mU(y) : y in Generators(S2)];
  in2 := [ &and[Evaluate(t[1], a) eq 1 : a in g2] : t in data ];
  gD2 := &+[Integers() | data[j][2] : j in [1..#data] | in2[j]];
  if gD2 lt 2 then continue; end if;
  cmA := [ &+[Integers() | data[j][4][i] : j in [1..#data] | in2[j]] : i in [1..#discs]];
  cmB := [ &+[Integers() | data[j][4][i] : j in [1..#data] | not in2[j]] : i in [1..#discs]];
  overlap := [discs[i] : i in [1..#discs] | cmA[i] gt 0 and cmB[i] gt 0];
  gcm := &+cmA;   // upper bound for g'_CM (sum over fields, may overcount)
  C1 := 2 * gcm lt gD2 - 1;
  anyp := 0; fullp := 0;
  for i in [1..#ps] do
    A := &*[Px | data[j][3][i] : j in [1..#data] | in2[j]];
    B := &*[Px | data[j][3][i] : j in [1..#data] | not in2[j]];
    if Degree(GCD(A, B)) eq 0 then
      if anyp eq 0 then anyp := ps[i]; end if;
      if fullp eq 0 and ps[i] in pinert then fullp := ps[i]; end if;
    end if;
    if anyp gt 0 and fullp gt 0 then break; end if;
  end for;
  status := fullp gt 0 select "FULL" else (anyp gt 0 and C1 select "NCM+C1" else "FAIL");
  printf "C4CHECK N=%o Delta=%o Delta'=%o index=%o g'=%o cm'(bound)=%o cm_complement=%o overlapD=%o coprime_any_p=%o coprime_inert_p=%o C1(X_Delta')=%o : %o\n",
    N, D, g2, #S2 div #S, gD2, cmA, cmB, overlap, anyp, fullp, C1, status;
end for;
print "DONE";
quit;
