// Check of the CM bound V(type) of code/cmV.py (Section 6 of the paper), run at all levels 11 <= N <= 200 (for
// N <= 10 and N = 12 there are no cusp forms of weight 2 on Gamma_1(N)) and at 24 levels above 200.
// For each Galois orbit of even characters chi mod N of order <= MAXORD, prints the local type tuple (same coding
// as perdelta.m / cmV.py), dim S_2(N,chi) and the CM proxy per CM discriminant (joint kernel of T_q over primes
// q <= B inert in Q(sqrt -D), q not | N; an upper bound for the CM-by-D dimension, per character).
// Usage: magma -b N:=256 MAXORD:=8 B:=200 recheck_V.m
SetColumns(0);
N := StringToInteger(N);
if not assigned MAXORD then MAXORD := "4"; end if; MAXORD := StringToInteger(MAXORD);
if not assigned B then B := "200"; end if; B := StringToInteger(B);
fac := Factorization(N); ps := [f[1] : f in fac]; es := [f[2] : f in fac];
discs := [D : D in Divisors(N) | D ge 3 and IsFundamentalDiscriminant(-D)];
G := DirichletGroup(N, CyclotomicField(EulerPhi(N)));
reps := [chi : chi in GaloisConjugacyRepresentatives(G) | Evaluate(chi, -1) eq 1 and Order(chi) le MAXORD];
function TypeOf(chi)
  comps := Decomposition(chi); t := [];
  for i in [1..#ps] do
    p := ps[i]; e := es[i];
    c := [x : x in comps | Modulus(x) eq p^e][1];
    a := Valuation(Conductor(c), p);
    if p eq 2 then
      if e eq 1 then Append(~t, 0); continue; end if;
      u := Evaluate(c, -1) eq -1 select 1 else 0;
      k := 0;
      if e ge 3 then z := Evaluate(c, 5); while z ne 1 do z := z^2; k +:= 1; end while; end if;
      Append(~t, u + 2*k);
    else
      quad := (a eq 1 and Order(c) eq 2) select 1 else 0;
      Append(~t, 2*a + quad);
    end if;
  end for;
  return t;
end function;
printf "LEVEL N=%o discs=%o orbits=%o\n", N, discs, #reps;
for chi in reps do
  chi0 := MinimalBaseRingCharacter(chi);
  M := CuspidalSubspace(ModularSymbols(chi0, 2, 1));
  d := Dimension(M); K := BaseRing(M);
  cm := [];
  for D in discs do
    if d eq 0 then Append(~cm, 0); continue; end if;
    W := VectorSpace(K, d);
    for q in [q : q in PrimesUpTo(B) | N mod q ne 0 and KroneckerSymbol(-D, q) eq -1] do
      W := W meet Kernel(HeckeOperator(M, q)); if Dimension(W) eq 0 then break; end if;
    end for;
    Append(~cm, Dimension(W));
  end for;
  printf "VAL N=%o type=%o order=%o cond=%o dim=%o cm=%o\n", N, TypeOf(chi), Order(chi0), Conductor(chi0), d, cm;
end for;
print "DONE";
quit;
