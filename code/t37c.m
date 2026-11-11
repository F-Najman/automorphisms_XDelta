load "lib.m";
N := 37; Delta := SubgroupList(N, [8]);
chars := CharsTrivialOn(N, Delta);
M := CuspForms(chars, 2); g := Dimension(M);
mu := MuDelta(N, Delta);
prec := Ceiling(6*mu/12) + 30;
B := qExpansionBasis(M, prec);
Q2, R := Relations(B, 2, Rationals());
Q3, R3 := Relations(B, 3, Rationals());
Q3 := [ Evaluate(c, [R.i : i in [1..g]]) : c in Q3 ];
P3 := ProjectiveSpace(R);
I2 := ideal<R | Q2>;
cub := [ c : c in Q3 | c notin I2 ];
C := Curve(P3, Q2 cat [cub[1]]);
pts := [ C!P : P in [ [0,-1,1,1], [0,0,1,1], [-1,-2,2,1], [1,0,0,0], [0,-1,1,0], [1,1,0,0] ] ];
Pinf := C![1,0,0,0];
cusps := [ C![1,0,0,0], C![0,0,1,1], C![-1,-2,2,1] ];
noncusps := [ C![0,-1,1,1], C![0,-1,1,0], C![1,1,0,0] ];
auts := Automorphisms(C);
inv := [ a : a in auts | a(Pinf) eq C![0,-1,1,1] ][1];
print "involution v:"; print DefiningEquations(inv);
// check it is an involution
comp := inv*inv;
printf "v^2 acts trivially on points: %o\n", &and[ comp(P) eq P : P in pts ];
// quotient by v
phi := iso<C -> C | DefiningEquations(inv), DefiningEquations(inv)>;
G := AutomorphismGroup(C, [phi]);
CG, prj := CurveQuotient(G);
printf "genus of X/<v> = %o\n", Genus(CG);
if Genus(CG) eq 1 then
  pt := prj(Pinf);
  E := EllipticCurve(CG, pt);
  Emin := MinimalModel(E);
  printf "X/<v> = elliptic curve %o, conductor %o, Cremona %o, torsion %o, rank %o\n", aInvariants(Emin), Conductor(Emin), CremonaReference(Emin), Invariants(TorsionSubgroup(Emin)), Rank(Emin);
end if;
// fixed points of v over Qbar: count via the fixed-point scheme

eqv := DefiningEquations(inv);
AmbR := CoordinateRing(Ambient(C));
fixeqs := [ eqv[i]*AmbR.j - eqv[j]*AmbR.i : i,j in [1..4] | i lt j ];
FixS := Scheme(C, fixeqs);
printf "fixed point scheme of v: dim %o, degree %o\n", Dimension(FixS), Degree(ReducedSubscheme(FixS));
// cuspidal class orders modulo good primes
for ell in [3,5,7,11,13] do
  Cl := ChangeRing(C, GF(ell));
  if not IsNonsingular(Cl) then printf "ell=%o singular model\n", ell; continue; end if;
  red := func< P | Cl![GF(ell)!c : c in Eltseq(P)] >;
  pl := func< P | Places(red(P))[1] >;
  function ord(D)
    for k in [1..500] do if IsPrincipal(k*D) then return k; end if; end for; return 0;
  end function;
  b2 := pl(cusps[2]) - pl(cusps[1]);
  b3 := pl(cusps[3]) - pl(cusps[1]);
  x2 := pl(noncusps[1]) - pl(cusps[1]);
  vb2 := pl(inv(cusps[2])) - pl(inv(cusps[1]));
  printf "ell=%o #X(F_ell)=%o : ord(b2)=%o ord(b3)=%o ord(x1-c1)=%o ord(b2+v(b2))=%o\n", ell, #Points(Cl), ord(b2), ord(b3), ord(x2), ord(b2+vb2);
end for;


quit;
