// Per-Delta check of (C1_Delta) at level N using the CM bound table from code/cmV_table.py.
// Test: sum over chi in Delta^perp (orbit reps x orbit size) of (2 V(type chi) - dim S_2(N,chi)) < -1.
// Usage: magma -b N:=192 T:=table_file perdelta.m
SetColumns(0);
N := StringToInteger(N);
tab := AssociativeArray();
for line in Split(Read(T), "\n") do
  w := Split(line, " "); if #w lt 2 then continue; end if;
  tab[[StringToInteger(x) : x in w[1..#w-1]]] := eval w[#w];
end for;
fac := Factorization(N); ps := [f[1] : f in fac]; es := [f[2] : f in fac];
G := DirichletGroup(N, CyclotomicField(EulerPhi(N)));
reps := [chi : chi in GaloisConjugacyRepresentatives(G) | Evaluate(chi, -1) eq 1];
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
vals := [];
for chi in reps do
  chi0 := MinimalBaseRingCharacter(chi);
  sz := EulerPhi(Order(chi));
  Append(~vals, <chi, sz * (2 * tab[TypeOf(chi)] - DimensionCuspForms(chi0, 2))>);
end for;
U, mU := UnitGroup(Integers(N)); m1 := (Integers(N)!(-1)) @@ mU;
nf := 0; ns := 0;
for rec in Subgroups(U) do
  S := rec`subgroup;
  if not m1 in S then continue; end if;
  idx := #U div #S; if idx lt 2 then continue; end if;
  ns +:= 1;
  gens := [Integers()!mU(x) : x in Generators(S)];
  tot := &+[RealField() | v[2] : v in vals | &and[Evaluate(v[1], a) eq 1 : a in gens]];
  if not (tot lt -1) then nf +:= 1; printf "PERDELTA_FAIL N=%o gens=%o idx=%o sum=%o\n", N, gens, idx, RealField(8)!tot; end if;
end for;
printf "PERDELTA N=%o subgroups=%o fail=%o\n", N, ns, nf;
quit;
