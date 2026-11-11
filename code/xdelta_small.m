// The scan of level N (Section 7 of the paper): the conditions (CM) [C1], a) [C3], b) [C3'] and c) [C2] of the
// proposition on correspondences, for all X_Delta(N) at one level N.
// Per Galois orbit of even characters chi: d = dim S_2(N,chi), W = CM proxy
//   (sum over imaginary quadratic F with d_F | N of the joint kernel of T_p, p inert in F, p <= B, p not | N;
//    W contains the CM part and is Hecke-stable, so dim W >= CM dim and tr(T_l T_l^* | W) >= S_l^CM),
//   s_l = tr(chi(l)^-1 T_l^2) on S_2(N,chi) and on W (absolute traces, i.e. summed over the orbit).
// Then every subgroup Delta of (Z/N)^* containing -1 is tested.
// xdelta_small.m: only characters of order <= MAXIDX and Delta of index <= MAXIDX.
// Usage: magma -b N:=63 L:=37 B:=60 MAXIDX:=3 xdelta_small.m      Output lines: ORBIT, DELTA ... OK|FAIL, SUMMARY, DONE
SetColumns(0);
if not assigned N then N := "63"; end if; N := StringToInteger(N);
if not assigned L then L := "37"; end if; L := StringToInteger(L);
if not assigned B then B := "60"; end if; B := StringToInteger(B);
if not assigned MAXIDX then MAXIDX := "3"; end if; MAXIDX := StringToInteger(MAXIDX);
ls := [l : l in PrimesUpTo(L) | N mod l ne 0];
discs := [D : D in Divisors(N) | D ge 3 and IsFundamentalDiscriminant(-D)];
G := DirichletGroup(N, CyclotomicField(EulerPhi(N)));
reps := [chi : chi in GaloisConjugacyRepresentatives(G) | Evaluate(chi, -1) eq 1 and Order(chi) le MAXIDX];
AbsTr := func<K, x | (K cmpeq Rationals()) select Rationals()!x else Rationals()!Trace(K!x)>;
data := [* *];
for chi in reps do
  chi0 := MinimalBaseRingCharacter(chi);
  M := CuspidalSubspace(ModularSymbols(chi0, 2, 1));
  d := Dimension(M); K := BaseRing(M); deg := (K cmpeq Rationals()) select 1 else Degree(K);
  s := [Rationals()!0 : l in ls]; sW := [Rationals()!0 : l in ls]; dW := 0;
  if d gt 0 then
    V := VectorSpace(K, d); W := sub<V | >;
    for D in discs do
      WD := V;
      for p in [p : p in PrimesUpTo(B) | N mod p ne 0 and KroneckerSymbol(-D, p) eq -1] do
        WD := WD meet Kernel(HeckeOperator(M, p));
        if Dimension(WD) eq 0 then break; end if;
      end for;
      W := W + WD;
    end for;
    dW := Dimension(W);
    for i in [1..#ls] do
      T := HeckeOperator(M, ls[i]); c := Evaluate(chi0, ls[i]);
      s[i] := AbsTr(K, Trace(T^2) / c);
      if dW gt 0 then
        Wb := BasisMatrix(W); X := Solution(Wb, Wb * T);
        sW[i] := AbsTr(K, Trace(X^2) / c);
      end if;
    end for;
  end if;
  Append(~data, <chi0, deg, d, dW, s, sW>);
  printf "ORBIT N=%o order=%o conductor=%o size=%o dim=%o cmproxy=%o\n", N, Order(chi0), Conductor(chi0), deg, d, dW;
end for;
g1 := &+[Integers() | data[j][2]*data[j][3] : j in [1..#data]];
printf "CHECK N=%o (partial genus, characters of order <= %o) %o\n", N, MAXIDX, g1;
psiN := N * &*[Rationals() | 1 + 1/p : p in PrimeDivisors(N)];
lmin := [l : l in PrimesUpTo(1000) | N mod l ne 0][1];
U, mU := UnitGroup(Integers(N));
m1 := (Integers(N)!(-1)) @@ mU;
nok := 0; nfail := 0;
for rec in Subgroups(U) do
  S := rec`subgroup;
  if not m1 in S then continue; end if;
  gens := [Integers()!mU(x) : x in Generators(S)];
  idx := #U div #S;
  if idx gt MAXIDX then continue; end if;
  inH := [ &and[Evaluate(data[j][1], a) eq 1 : a in gens] : j in [1..#data] ];
  g := &+[Integers() | data[j][2]*data[j][3] : j in [1..#data] | inH[j]];
  if g lt 2 then continue; end if;
  gW := &+[Integers() | data[j][2]*data[j][4] : j in [1..#data] | inH[j]];
  SH := [ &+[Rationals() | data[j][5][i] : j in [1..#data] | inH[j]] : i in [1..#ls]];
  SWH := [ &+[Rationals() | data[j][6][i] : j in [1..#data] | inH[j]] : i in [1..#ls]];
  C1 := 2*gW lt g - 1;
  l3 := [ls[i] : i in [1..#ls] | C1 and SH[i] gt (ls[i]+1)^2];
  l3p := [ls[i] : i in [1..#ls] | SH[i] - 2*SWH[i] gt (ls[i]+1)^2];
  C2 := C1 and (7/800) * psiN * idx gt 2*(lmin + 1);
  ok := #l3 gt 0 or #l3p gt 0 or C2;
  if ok then nok +:= 1; else nfail +:= 1; end if;
  printf "DELTA N=%o gens=%o idx=%o genus=%o cmproxy=%o C1=%o C3_l=%o C3p_l=%o C2=%o %o\n", N, gens, idx, g, gW,
    C1, #l3 gt 0 select l3[1] else 0, #l3p gt 0 select l3p[1] else 0, C2, ok select "OK" else "FAIL";
end for;
printf "SUMMARY N=%o ok=%o fail=%o\n", N, nok, nfail;
print "DONE";
quit;
