load "lib.m";
N := 21; Delta := SubgroupList(N, [8]);
print "Delta =", Delta, " mu =", MuDelta(N, Delta);
M, B := QBasis(N, Delta, 300);
print "genus", #B;
for b in B do print b; end for;
Q2, R := Relations(B, 2, Rationals());
print "quadric relations:", Q2;
Q4, R := Relations(B, 4, Rationals());
print "#quartic relations:", #Q4;
// hyperelliptic model: find conic, parametrize
if #Q2 eq 1 then
  P2 := ProjectiveSpace(R);
  Con := Conic(P2, Q2[1]);
  print "conic has rational point:", HasRationalPoint(Con);
  phi := Parametrization(Con);
  print phi;
end if;
quit;
