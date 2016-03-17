
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

xFactor := function(f)
  return Min([Exponents(t)[1] : t in Terms(f) | t ne 0]);
end function;

yFactor := function(f)
  return Min([Exponents(t)[2] : t in Terms(f) | t ne 0]);
end function;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

FacePolynomial := function(F)
  f := FaceFunction(F);
  P<Z> := PolynomialRing(CoefficientRing(Parent(f)));
  f := Evaluate(f, <1, Z>);
  E, _ := Support(f); C := Coefficients(f);
  n := GradientVector(F)[1]; b := Integers()!EndVertices(F)[2][2];
  return &+[C[e + 1]*Z^((e - b) div n) : e in E];
end function;

forward NewtonPuiseuxAlgorithmLoop;

intrinsic NewtonPuiseuxAlgorithm(f::RngMPolElt : Terms := -1) -> [ ]
{ Computes the Puiseux expansion of any bivariate polynomial }
require Rank(Parent(f)) eq 2: "Argument must be a bivariate polynomial";
  // If Nf start on the right of the x-axis, we have an x-factor.
  yBranch := (xFactor(f) gt 0) select [*<Parent(f).1, xFactor(f)>*] else [* *];

  P<x, y> := PolynomialRing(AlgebraicClosure(CoefficientRing(Parent(f))), 2);
  S := NewtonPuiseuxAlgorithmLoop(P!SquarefreePart(f),
    [<P!g[1], g[2], 1> : g in SquarefreeFactorization(f)], 1, Terms - 1);
  return yBranch cat SequenceToList([<s[1], s[2][1][1]> : s in S]);
end intrinsic;

intrinsic NewtonPuiseuxAlgorithm(L::[RngMPolElt] : Terms := -1) -> []
{ Computes the Puiseux expansion for the product of all the elements of L }
require #L gt 0: "Argument must be a non-empty list";
require &and[Rank(Parent(f)) eq 2 : f in L]:
  "Elements of L must be bivariate polynomials";

  f := &*L;
  // If Nf start on the right of the x-axis, we have an x-factor.
  yBranch := (xFactor(f) gt 0) select [*<Parent(f).1,
    [<xFactor(L[i]), i> : i in [1..#L] | xFactor(L[i]) ne 0]>*] else [* *];

  P<x, y> := PolynomialRing(AlgebraicClosure(CoefficientRing(Parent(f))), 2);
  sqFreePart := P!SquarefreePart(f);
  sqFreeFact := [[<P!g[1], g[2], i>: g in SquarefreeFactorization(L[i])][1] :
    i in [1..#L]];
  S := NewtonPuiseuxAlgorithmLoop(sqFreePart, sqFreeFact, 1, Terms - 1);
  return yBranch cat SequenceToList(S);
end intrinsic;

NewtonPuiseuxAlgorithmLoop := function(f, L, ord, terms)
  A := CoefficientRing(Parent(f)); x0 := Parent(f).1; y0 := Parent(f).2;
  Q<x> := PuiseuxSeriesRing(A: Precision := ord);
  // Step (i.a): Select only those factors containing the 0 branch.
  S := yFactor(f) gt 0 select [<Q!0, [<g[2], g[3]> : g in L
    | yFactor(g[1]) ne 0]>] else [];
  // Step (i.b): For each side...
  for F in Faces(NewtonPolygon(f)) do
    n := GradientVector(F)[1]; m := GradientVector(F)[2];
    P<x> := PuiseuxSeriesRing(A: Precision := ord*n);
    // Apply the change of variables (1).
    C := Reverse(Coefficients(n eq 1 select f else Evaluate(f, 1, x0^n), 2));
    CL := [ <Reverse(Coefficients(n eq 1 select g[1] else
      Evaluate(g[1], 1, x0^n), 2)), g[2], g[3]> : g in L];
    // For each root...
    for a in [<Root(a[1], n), a[2]> : a in Roots(FacePolynomial(F))] do
      // Apply the change of variables (2) & get the sub-solution recursively.
      ff := [i gt 1 select C[i] + Self(i-1)*x0^m*(a[1] + y0) else C[1] :
        i in [1..#C]][#C];
      LL := [ <[i gt 1 select Cj[1][i] + Self(i-1)*x0^m*(a[1] + y0) else
        Cj[1][1] : i in [1..#Cj[1]]][#Cj[1]], Cj[2], Cj[3]> : Cj in CL ];
      // Select only those factors that contain the current branch.
      LL := [g : g in LL | Vertices(NewtonPolygon(g[1]:
        Faces := "All"))[1][2] ne 0];
      // If the mult. of a is greater than 1 continue.
      R := a[2] ne 1 or terms gt 0 select
        NewtonPuiseuxAlgorithmLoop(ff, LL, ord*n, terms - 1) else
          [<P!0, [<g[2], g[3]> : g in LL]>];
      // Undo the change of variables.
      S cat:= [<x^(m/n)*(a[1] + Composition(s[1], x^(1/n))), s[2]> :
        s in R];
    end for;
  end for;

  return S;
end function;

forward NewtonPuiseuxAlgorithmReducedLoop;

intrinsic NewtonPuiseuxAlgorithmReduced(f::RngMPolElt : Terms := -1) -> [ ]
{ Computes the Puiseux expansion of a reduced bivariate polynomial }
require Rank(Parent(f)) eq 2: "Argument must be a bivariate polynomial";
  // If Nf start on the right of the x-axis, we have an x-factor.
  yBranch := (xFactor(f) gt 0) select [* Parent(f).1 *] else [* *];

  P := PolynomialRing(AlgebraicClosure(CoefficientRing(Parent(f))), 2);
  S := NewtonPuiseuxAlgorithmReducedLoop(P!SquarefreePart(f), 1, Terms - 1);
  return yBranch cat SequenceToList(S);
end intrinsic;

NewtonPuiseuxAlgorithmReducedLoop := function(f, ord, terms)
  A := CoefficientRing(Parent(f)); x0 := Parent(f).1; y0 := Parent(f).2;
  Q<x> := PuiseuxSeriesRing(A: Precision := ord);
  // Step (i.a): Select only those factors containing the 0 branch.
  S := yFactor(f) gt 0 select [Q!0] else [];
  // Step (i.b): For each side...
  for F in Faces(NewtonPolygon(f)) do
    n := GradientVector(F)[1]; m := GradientVector(F)[2];
    P<x> := PuiseuxSeriesRing(A: Precision := ord*n);
    // Apply the change of variables (1).
    C := Reverse(Coefficients(n eq 1 select f else Evaluate(f, 1, x0^n), 2));
    // For each root...
    for a in [<Root(a[1], n), a[2]> : a in Roots(FacePolynomial(F))] do
      // Apply the change of variables (2) & get the sub-solution recursively.
      R := a[2] ne 1 or terms gt 0 select
        NewtonPuiseuxAlgorithmReducedLoop([i gt 1 select
          C[i] + Self(i-1)*x0^m*(a[1] + y0) else C[1] : i in [1..#C]][#C],
        ord*n, terms - 1) else [P!0];
      // Undo the change of variables.
      S cat:= [x^(m/n)*(a[1] + Composition(s, x^(1/n))) : s in R];
    end for;
  end for;

  return S;
end function;
