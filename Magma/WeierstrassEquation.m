
intrinsic WeierstrassEquation(s::RngSerPuisElt) -> RngMPolElt
{ Computes the Weierstrass equation associated to a Puiseux series }
require Type(CoefficientRing(Parent(s))) eq FldAC:
  "Argument's coefficient ring must be an algebraically closed field";
  // Make a hard-copy of the AlgebraicClosure field.
  A := CoefficientRing(Parent(s)); B := AlgebraicClosure(BaseField(A));
  P<z> := PolynomialRing(B);
  for g in Generators(Ideal(A)) do
    Roots(Evaluate(UnivariatePolynomial(g), z));
  end for;
  // Perform all the computations in the new field B.
  n := ExponentDenominator(s); m := IntegerRing()!(n*Valuation(s));
  Q<x0, y0> := PolynomialRing(B, 2); zeta := RootOfUnity(n, B);
  C := Coefficients(s); S := [];
  for j in [1..n] do
    Append(~S, &+[B!C[i]*zeta^(j*(m + i - 1))*x0^(m + i - 1) : i in [1..#C]]);
  end for;
  // Go back to the original AlgebraicClosure A.
  Q<x, y> := PolynomialRing(A, 2);
  return &+[A!LeadingCoefficient(T) * x^(Exponents(T)[1] div n)
   * y^(Exponents(T)[2]) : T in Terms(&*[y0 - s : s in S])];
end intrinsic;
