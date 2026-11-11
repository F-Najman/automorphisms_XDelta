// For each canon2.m run of Section 8.1 of the paper (the 28 curves of Table 1 and X_1(N), N = 21, 22, 24, 30), check that
// the q-expansion basis used by canon2.m has integral coefficients and stays linearly independent modulo ell,
// so that its reduction is a basis of the differentials of X_{F_ell}.
// Run from code/:  magma -b review_manuscript_reduction.m < /dev/null
SetColumns(0);
load "lib.m";
// N, ell, g, Delta (from the RESULT lines of logs/canon2/; see code/review_manuscript_logs.py)
runs := [
<24,193,3,[1,7,17,23]>, <24,193,3,[1,5,19,23]>, <25,101,4,[1,7,18,24]>, <26,131,4,[1,3,9,17,23,25]>,
<26,131,4,[1,5,21,25]>, <28,113,4,[1,13,15,27]>, <28,113,4,[1,3,9,19,25,27]>, <30,151,5,[1,11,19,29]>,
<31,311,6,[1,5,6,25,26,30]>, <32,193,5,[1,15,17,31]>, <33,199,5,[1,2,4,8,16,17,25,29,31,32]>,
<36,109,7,[1,17,19,35]>, <36,109,3,[1,11,13,23,25,35]>, <40,241,7,[1,3,9,13,27,31,37,39]>,
<40,241,5,[1,9,11,19,21,29,31,39]>, <42,127,9,[1,5,17,25,37,41]>, <44,353,8,[1,5,7,9,19,25,35,37,39,43]>,
<45,181,5,[1,4,11,14,16,19,26,29,31,34,41,44]>, <48,193,7,[1,5,19,23,25,29,43,47]>,
<48,193,7,[1,7,17,23,25,31,41,47]>, <48,193,5,[1,11,13,23,25,35,37,47]>,
<49,197,3,[1,6,8,13,15,20,22,27,29,34,36,41,43,48]>, <50,101,4,[1,9,11,19,21,29,31,39,41,49]>,
<54,109,10,[1,17,19,35,37,53]>, <56,113,9,[1,3,9,19,25,27,29,31,37,47,53,55]>,
<64,193,5,[1,7,9,15,17,23,25,31,33,39,41,47,49,55,57,63]>, <72,433,9,[1,11,13,23,25,35,37,47,49,59,61,71]>,
<37,149,4,[1,6,8,10,11,14,23,26,27,29,31,36]>, <21,127,5,[1,20]>, <22,199,6,[1,21]>, <30,151,9,[1,29]>,
<24,193,5,[1,23]>
];
nbad := 0;
for r in runs do
  N, ell, g, Delta := Explode(r);
  M := CuspForms(CharsTrivialOn(N, Delta), 2);
  assert Dimension(M) eq g;
  prec := Ceiling(2*MuDelta(N, Delta)/12) + 20;
  B := qExpansionBasis(M, prec);
  coeffs := [ [ Coefficient(b, i) : i in [0..prec-1] ] : b in B ];
  integral := &and[ IsCoercible(Integers(), c) : c in row, row in coeffs ];
  rk := Rank(Matrix(GF(ell), [ [ GF(ell)!c : c in row ] : row in coeffs ]));
  ok := integral and rk eq g;
  if not ok then nbad +:= 1; end if;
  printf "N=%o ell=%o g=%o |Delta|=%o integral=%o rank_mod_ell=%o %o\n", N, ell, g, #Delta, integral, rk, ok select "OK" else "FAIL";
end for;
printf "RUNS %o FAILURES %o\n", #runs, nbad;
if nbad eq 0 then print "ALL REDUCTION CHECKS PASSED"; end if;
quit;
