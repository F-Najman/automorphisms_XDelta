// For the non-iota-stable Delta at N=65: do the characters chi (trivial on Delta) with S_2(65,chi) != 0 cut out Delta?
load "lib.m";
N := 65;
for gens in [[28],[33],[8],[18]] do
  Delta := SubgroupList(N, gens);
  chars := CharsTrivialOn(N, Delta);
  occ := [ chi : chi in chars | Dimension(CuspForms(chi,2)) gt 0 ];
  U, m := UnitGroup(Integers(N));
  ker := [ Integers()!(m(u)) : u in U | &and[ Evaluate(chi, Integers()!(m(u))) eq 1 : chi in occ ] ];
  printf "Delta=<%o> |Delta|=%o g=%o : %o of %o characters have S_2 != 0 ; intersection of kernels = Delta ? %o\n",
    gens, #Delta, Dimension(CuspForms(chars,2)), #occ, #chars, Sort(ker) eq Delta;
  printf "   dims: %o\n", [ Dimension(CuspForms(chi,2)) : chi in chars ];
end for;
quit;
