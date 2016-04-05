import "SemiGroup.m": Euclides, TailExponentMatrix;

intrinsic WeierstrassEquation(s::RngSerPuisElt) -> RngMPolElt
{ Computes the Weierstrass equation associated to a Puiseux series }
require Type(CoefficientRing(Parent(s))) eq FldAC:
  "Argument's coefficient ring must be an algebraically closed field";
  A := CoefficientRing(Parent(s)); Q<x, y> := PolynomialRing(A, 2);
  n := ExponentDenominator(s); ZZ := IntegerRing();
  m := Valuation(s) eq Infinity() select ZZ!0 else ZZ!(n*Valuation(s));
  zeta := RootOfUnity(n, A); C := Coefficients(s); S := [];
  for j in [1..n] do
    Append(~S, &+[Q | C[i]*zeta^(j*(m + i - 1))*x^(m + i - 1) : i in [1..#C]]);
  end for;
  return &+[LeadingCoefficient(T) * x^(Exponents(T)[1] div n)
   * y^(Exponents(T)[2]) : T in Terms(&*[y - s : s in S])];
end intrinsic;

ClusterFactorization := function(P, v, c)
  n := Ncols(P); B := [* *]; ZZ := IntegerRing(); N := Transpose(P)*P;
  // For each point with strictly positive excess.
  for i in [i : i in [1..n] | (v*N)[1][i] gt 0] do
    p := i; I := [p]; // Traverse the cluster back to the origin
    while p ne 1 do
      p := [j : j in Reverse([1..n]) | P[p][j] eq -1][1]; I := [p] cat I;
    end while;
    e := ZeroMatrix(ZZ, 1, #I); e[1][#I] := (v*N)[1][i];
    Q := Submatrix(P, I, I); e := e*Q^-1;
    B cat:= [* <Q, e*Transpose(Q^-1), c[I]> *];
  end for; return B;
end function;

SharplyCurve := function(P, v, c)
  m := Gcd(Eltseq(v)); v := v div m; // Mult. of the irreducible cluster
  G := SemiGroup(P, v); M := CharExponents(G); T := TailExponentMatrix(P, v);
  if #c gt 1 and Type(c[2]) eq Infty then error "Inverted curve!"; end if; //TEMP!

  P<t> := PuiseuxSeriesRing(Parent(c[1])); s := P!0; k := 1; n := M[1][2];
  if T ne -1 then M cat:= [<T, 1>]; else s +:= 1*t^(M[#M][1]/n); end if;
  for i in [2..#M] do
    mj := M[i-1][1]; nj := M[i-1][2]; mi := M[i][1]; h0 := (mi - mj) div nj;
    s +:= &+[P | c[k + l]*t^((mj + l*nj)/n) : l in [0..h0]];
    k +:= &+Euclides(mi - mj, nj)[1];
  end for; return WeierstrassEquation(s)^m;
end function;

intrinsic FooBar(X::<>) -> []
{ FooBar }
  //print ClusterFactorization(X[1], X[2], X[3]);
  //print Parent(X[3][1]);
  return [* SharplyCurve(Y[1], Y[2], Y[3]) :
    Y in ClusterFactorization(X[1], X[2], X[3]) *];
end intrinsic;
