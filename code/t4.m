load "lib.m"; load "hyp.m";
N := 21; Delta := SubgroupList(N, [8]);
M, B := QBasis(N, Delta, 400);
// check B[3] = x^2 B[1] relation form: B1*B3 = B2^2
f, x, y := HypModelFromBasis(B, 3);
print "y^2 =", f; print "factorization:", Factorization(f);
C := HyperellipticCurve(f);
print "genus", Genus(C);
G := GeometricAutomorphismGroup(C);
print "Geometric Aut X_{+-8}(21): order", #G, GroupName(G);
print "Aut over Q:", #AutomorphismGroup(C);
// Hecke data: newforms and a_3, a_7
for chi in CharsTrivialOn(N, Delta) do
  S := CuspForms(chi, 2);
  nf := Newforms(S);
  for orb in nf do
    f0 := orb[1];
    printf "char order %o cond %o : newform dim %o level %o, a3=%o a7=%o a2=%o\n", Order(chi), Conductor(chi), #orb, Level(f0), Coefficient(f0,3), Coefficient(f0,7), Coefficient(f0,2);
  end for;
end for;
quit;
