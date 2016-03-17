
GenericEquation := function(M)
  n := M[1]; M[1] := 0;
  ExpChar := [ M[i]/n : i in [2..#M] ];
  N := [ i gt 1 select Gcd(Self(i-1), M[i]) else n : i in [1..#M] ];

  if (N[#N] ne 1) or (N[#N] eq N[#N-1]) then
    print "Not a valid set of characteristic exponents.";
  end if;

  Kn := CyclotomicField(n);
  zeta_n := Kn.1;

  S := [ Floor(M[i+1]/N[i]) - Ceiling(M[i]/N[i]) + 1: i in [2..#M-1] ];
  if #S eq 0 then num_params := 1; else num_params := &+S + 1; end if;
  A := PolynomialRing(Kn, num_params);
  AssignNames(~A, [ "a" * IntegerToString(i) : i in [1..num_params] ]);

  P<x> := PolynomialRing(A);
  s := 0; j := 1;
  for i in [2..#M-1] do
    hi := (M[i+1] - M[i]) div N[i];
    for l in [0..hi] do
      s +:= A.j*x^(M[i] + l*N[i]);
      j +:= 1;
    end for;
  end for;
  s +:= A.j*x^M[#M];
  S := [Evaluate(s, zeta_n^i*x) : i in [1..n]];

  Q<x> := PuiseuxSeriesRing(A, 100);
  R<y> := PolynomialRing(Q);
  S := [y - Composition(Q!s, x^(1/n)) : s in S];

  return &*S;

end function;

WeierstrassEquation := function(s)
  n := ExponentDenominator(s);
  m := IntegerRing()!(n*Valuation(s));
  P<z> := PolynomialRing(CoefficientRing(Parent(s)));
  zetas := Roots(z^n - 1);
  C := Coefficients(s);
  t := Parent(s).1;
  S := [];
  for zeta in zetas do
    Append(~S, &+[C[i]*zeta[1]^(m + i - 1)*t^((m + i -1)/n) : i in [1..#C]]);
  end for;
  Q<y> := PolynomialRing(Parent(s));
  f := &*[y - s : s in S];
  return f;
end function;
