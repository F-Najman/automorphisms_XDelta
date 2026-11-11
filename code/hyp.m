// hyperelliptic model from basis f1,f2,f3... with canonical image rational normal curve f_i = x^{i-1} f_1
function HypModelFromBasis(B, g)
  L<t> := LaurentSeriesRing(Rationals());
  f1 := L!B[1]; f2 := L!B[2];
  x := f2/f1;
  y := t*Derivative(x)/f1;
  vals := [ x^i : i in [0..2*g+2] ] cat [ y^2 ];
  v0 := Min([Valuation(v) : v in vals]);
  mat := Matrix(Rationals(), [ [Coefficient(v, i) : i in [v0..v0+200]] : v in vals ]);
  K := Kernel(mat);
  assert Dimension(K) eq 1;
  b := Eltseq(Basis(K)[1]); b := [c/b[#b] : c in b];
  Px<X> := PolynomialRing(Rationals());
  f := -&+[ b[i+1]*X^i : i in [0..2*g+2] ];
  return f, x, y;
end function;
