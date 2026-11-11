load "lib.m";
if not assigned Nmax then Nmax := 50; elif Type(Nmax) eq MonStgElt then Nmax := StringToInteger(Nmax); end if;
function IotaQ(N, a, Q)
  a := Integers()!a;
  inv := Integers()!((Integers(Q)!a)^-1);
  return Integers()!CRT([inv, a mod (N div Q)], [Q, N div Q]);
end function;
for N in [5..Nmax] do
  if not IsSquarefree(N) then continue; end if;
  U, m := UnitGroup(Integers(N));
  mone := (Integers(N)!(-1)) @@ m;
  subs := [ H`subgroup : H in Subgroups(U) | mone in H`subgroup ];
  for H in subs do
    Delta := Sort([ Integers()!(m(h)) : h in H ]);
    if #Delta eq EulerPhi(N) then continue; end if; // skip X_0(N)
    gens := [ Integers()!(m(H.i)) : i in [1..Ngens(H)] ];
    chars := CharsTrivialOn(N, Delta);
    g := Dimension(CuspForms(chars, 2));
    if g lt 2 then continue; end if;
    Qs := [ Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1 and Q gt 1 ];
    stab := [ Q : Q in Qs | Sort([ IotaQ(N,a,Q) mod N : a in Delta ]) eq Delta ];
    nonstab := [ Q : Q in Qs | Q notin stab ];
    printf "N=%o |Delta|=%o index=%o gens=%o g=%o stableQ=%o NONSTABLE=%o\n", N, #Delta, EulerPhi(N) div #Delta, gens, g, stab, nonstab;
  end for;
end for;
quit;
