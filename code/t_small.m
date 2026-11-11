load "lib.m"; load "hyp.m";
// X_1(16), X_1(18): genus 2, hyperelliptic; X_1(20): genus 3
for N in [16, 18] do
  Delta := SubgroupList(N, [1]);
  M, B := QBasis(N, Delta, 300);
  f, x, y := HypModelFromBasis(B, 2);
  C := HyperellipticCurve(f);
  G := GeometricAutomorphismGroup(C);
  printf "X_1(%o): y^2 = %o ; genus %o ; |Aut_Qbar| = %o (%o) ; |Aut_Q| = %o ; |B| = %o\n", N, f, Genus(C), #G, GroupName(G), #AutomorphismGroup(C), (EulerPhi(N) div 2) * 2^#PrimeDivisors(N);
end for;
N := 20; Delta := SubgroupList(N, [1]);
M, B := QBasis(N, Delta, 300);
Q2, R := Relations(B, 2, Rationals());
Q4, R4 := Relations(B, 4, Rationals());
printf "X_1(20): genus %o, #quadrics %o, #quartics %o\n", #B, #Q2, #Q4;
C := Curve(ProjectiveSpace(R4), Q4[1]);
printf "plane quartic: %o ; genus %o\n", Q4[1], Genus(C);
G := GeometricAutomorphismGroup(C);
printf "X_1(20): |Aut_Qbar| = %o (%o); |Aut_Q| = %o ; |B| = %o\n", #G, GroupName(G), #AutomorphismGroup(C), (EulerPhi(N) div 2) * 2^#PrimeDivisors(N);
quit;
