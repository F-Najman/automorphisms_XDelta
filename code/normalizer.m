// |N(Gamma_Delta) ∩ N(Gamma_0(N)) / ±Gamma_Delta| for Gamma_1(N) ⊆ Gamma_Delta ⊆ Gamma_0(N).
// N(Gamma_0(N))/Gamma_0(N) is enumerated from the Atkin-Lehner matrices w_Q (Q || N) and S_h = (1 1/h; 0 1),
// h the largest divisor of 24 with h^2 | N (Atkin-Lehner-Newman). g normalises Gamma_Delta iff g x g^-1 is in
// Gamma_Delta for generators x of Gamma_Delta (Gamma_1(N) generators and one matrix per generator of Delta).
// Usage: magma -b N:=36 D:=23 normalizer.m      (D as for canon.m; -1 is always added)
load "lib.m";
if not assigned N then N := "36"; end if; N := StringToInteger(N);
if not assigned D then D := "23"; end if;
gens := [StringToInteger(s) : s in Split(D, ",")];
Delta := SubgroupList(N, gens);
MQ := MatrixAlgebra(Rationals(), 2);
function Prim(z) // primitive integer matrix proportional to z (rational 2x2), up to sign
  e := Eltseq(z); den := LCM([Denominator(t) : t in e]); ei := [Integers()!(t*den) : t in e];
  g := GCD(ei); return [t div g : t in ei];
end function;
function InG0(z)
  p := Prim(z); return p[1]*p[4]-p[2]*p[3] eq 1 and p[3] mod N eq 0;
end function;
function InGD(z)
  p := Prim(z); return p[1]*p[4]-p[2]*p[3] eq 1 and p[3] mod N eq 0 and (p[4] mod N) in Delta;
end function;
GD := [MQ!Eltseq(x) : x in Generators(Gamma1(N))];
for d in gens cat [-1] do
  a := Integers()!((Integers(N)!d)^-1); b := (a*d - 1) div N;
  Append(~GD, MQ![a, b, N, d]);
end for;
assert &and[InGD(x) : x in GD];
// generators of N(Gamma_0(N)) modulo Gamma_0(N)
NG := [];
for Q in [Q : Q in Divisors(N) | Q gt 1 and GCD(Q, N div Q) eq 1] do
  _, u, v := XGCD(Q, N div Q);          // Q u + (N/Q) v = 1
  Append(~NG, MQ![Q*u, -v, N, Q]);      // det = Q^2 u + N v = Q
end for;
h := Max([h : h in Divisors(24) | N mod h^2 eq 0]);
if h gt 1 then Append(~NG, MQ![1, 1/h, 0, 1]); end if;
// enumerate N(Gamma_0(N))/Gamma_0(N)
reps := [MQ!1]; queue := [MQ!1];
while #queue gt 0 do
  x := queue[1]; Remove(~queue, 1);
  for s in NG do
    y := x * s;
    if not &or[InG0(y * r^-1) : r in reps] then Append(~reps, y); Append(~queue, y); end if;
  end for;
end while;
k := #[g : g in reps | &and[InGD(g * x * g^-1) : x in GD]];
Dord := EulerPhi(N) div #Delta;
function IotaQ(a, Q) a := Integers()!a; inv := Integers()!((Integers(Q)!a)^-1);
  return Integers()!CRT([inv, a mod (N div Q)], [Q, N div Q]); end function;
stableQ := [Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1 and (Q eq 1 or Sort([IotaQ(a,Q) mod N : a in Delta]) eq Delta)];
printf "NORMALIZER N=%o Delta=%o h=%o |N(G0)/G0|=%o normalising=%o |D_Delta|=%o |Norm/Gamma_Delta|=%o |B_Delta|=%o\n",
  N, Delta, h, #reps, k, Dord, k*Dord, Dord*#stableQ;
quit;
