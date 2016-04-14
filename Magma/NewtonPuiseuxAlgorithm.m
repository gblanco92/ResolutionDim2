
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

xFactor := function(f)
  return Min([IntegerRing() | Exponents(t)[1] : t in Terms(f) | t ne 0]);
end function;

yFactor := function(f)
  return Min([IntegerRing() | Exponents(t)[2] : t in Terms(f) | t ne 0]);
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

intrinsic NewtonPuiseuxAlgorithm(f::RngMPolElt : Terms := -1,
                                                 Polynomial := false) -> [ ]
{ Computes the Puiseux expansion of any bivariate polynomial }
require Rank(Parent(f)) eq 2: "Argument must be a bivariate polynomial";
  // If Nf start on the right of the x-axis, we have an x-factor.
  yBranch := (xFactor(f) gt 0) select [*<Parent(f).1,
    [<xFactor(f), 1>], Parent(f).1>*] else [* *];

  P<x, y> := PolynomialRing(AlgebraicClosure(CoefficientRing(Parent(f))), 2);
  S := yBranch cat SequenceToList(NewtonPuiseuxAlgorithmLoop(P!SquarefreePart(f),
    [<P!g[1], g[2], 1> : g in SquarefreeFactorization(f)], 1, Terms - 1));
print S;
  if not Polynomial then return [* <s[1], s[2][1][1]> : s in S *];
  else return [* <s[1], s[2][1][1], s[3]> : s in S *]; end if;
end intrinsic;

intrinsic NewtonPuiseuxAlgorithm(L::[RngMPolElt] : Terms := -1,
                                                   Polynomial := false) -> []
{ Computes the Puiseux expansion for the product of all the elements of L }
require #L gt 0: "Argument must be a non-empty list";
require &and[Rank(Parent(f)) eq 2 : f in L]:
  "Elements of L must be bivariate polynomials";
  f := &*L;
  P<x, y> := PolynomialRing(AlgebraicClosure(CoefficientRing(Parent(f))), 2);
  // If Nf start on the right of the x-axis, we have an x-factor.
  yBranch := (xFactor(f) gt 0) select [*<Parent(f).1,
    [<xFactor(L[i]), i> : i in [1..#L] | xFactor(L[i]) ne 0], x>*] else [* *];

  sqFreePart := P!SquarefreePart(f); sqFreeFact := [];
  for i in [1..#L] do
  sqFreeFact cat:= [<P!g[1], g[2], i>: g in SquarefreeFactorization(L[i]) |
    Evaluate(L[i], <0, 0>) eq 0];
  end for;
  S := yBranch cat SequenceToList(NewtonPuiseuxAlgorithmLoop(sqFreePart,
    sqFreeFact, 1, Terms - 1));

  // Return the polynomial residue if requested.
  if not Polynomial then return [*<s[1], s[2]> : s in S*];
  else return S; end if;
end intrinsic;

NewtonPuiseuxAlgorithmLoop := function(f, L, ord, terms)
  Q<x> := PuiseuxSeriesRing(CoefficientRing(Parent(f)));
  x0 := Parent(f).1; y0 := Parent(f).2;
  // Step (i.a): Select only those factors containing the 0 branch.
  S := yFactor(f) gt 0 select [<Q!0, [<g[2], g[3]> : g in L
    | yFactor(g[1]) ne 0], y0>] else [];
  // Step (i.b): For each side...
  for F in Faces(NewtonPolygon(f)) do
    n := GradientVector(F)[1]; m := GradientVector(F)[2];
    // Apply the change of variables (1).
    C := Reverse(Coefficients(n eq 1 select f else Evaluate(f, 1, x0^n), 2));
    CL := [<Reverse(Coefficients(n eq 1 select g[1] else
      Evaluate(g[1], 1, x0^n), 2)), g[2], g[3]> : g in L];
    // For each root...
    for a in [<Root(a[1], n), a[2]> : a in Roots(FacePolynomial(F))] do
      // Apply the change of variables (2) & get the sub-solution recursively.
      ff := [i gt 1 select C[i] + Self(i-1)*x0^m*(a[1] + y0) else C[1] :
        i in [1..#C]][#C];
      LL := [<[i gt 1 select Cj[1][i] + Self(i-1)*x0^m*(a[1] + y0) else
        Cj[1][1] : i in [1..#Cj[1]]][#Cj[1]], Cj[2], Cj[3]> : Cj in CL];
      // Select only those factors that contain the current branch.
      LL := [g : g in LL | Vertices(NewtonPolygon(g[1]:
        Faces := "All"))[1][2] ne 0];
      // If the mult. of a is greater than 1 continue.
      R := (a[2] ne 1 and terms lt -1) or terms gt 0 select
        NewtonPuiseuxAlgorithmLoop(ff, LL, ord*n, terms - 1) else
          [<Q!0, [<g[2], g[3]> : g in LL], ff>];
      // Undo the change of variables.
      S cat:= [<x^(m/n)*(a[1] + ChangePrecision(Composition(s[1], x^(1/n)),
        Infinity())), s[2], s[3]> : s in R];
    end for;
  end for;

  return S;
end function;

forward NewtonPuiseuxAlgorithmReducedLoop;

intrinsic NewtonPuiseuxAlgorithmReduced(f::RngMPolElt : Terms := -1,
                                        Polynomial := false) -> [ ]
{ Computes the Puiseux expansion of a reduced bivariate polynomial }
require Rank(Parent(f)) eq 2: "Argument must be a bivariate polynomial";
  P := PolynomialRing(AlgebraicClosure(CoefficientRing(Parent(f))), 2);
  // If Nf start on the right of the x-axis, we have an x-factor.
  yBranch := (xFactor(f) gt 0) select [* <Parent(f).1, P!0> *] else [* *];

  S := yBranch cat SequenceToList(NewtonPuiseuxAlgorithmReducedLoop(
    P!SquarefreePart(f), 1, Terms - 1));
  if Polynomial then return S; else return [* s[1] : s in S *]; end if;
end intrinsic;

intrinsic NewtonPuiseuxAlgorithmExpandReduced(s::RngSerPuisElt, f::RngMPolElt
                                              : Terms := 1) -> [ ]
{ Expands the Puiseux expansion s of a reduced bivariate polynomial }
require Rank(Parent(f)) eq 2: "Argument f must be a bivariate polynomial";

  n := ExponentDenominator(s); x := Parent(s).1;
  m := s eq 0 select 0 else Degree(s);
  S := Terms gt 0 select NewtonPuiseuxAlgorithmReducedLoop(f, n, Terms - 1)
       else [<PuiseuxSeriesRing(CoefficientRing(Parent(f)))!0, f>];
  P<x> := PuiseuxSeriesRing(CoefficientRing(Parent(f)));
  return
    [s + x^m*Composition(ChangePrecision(si[1], Infinity()), x^(1/n)): si in S];
end intrinsic

intrinsic NewtonPuiseuxAlgorithmExpandReduced(x::RngMPolElt, f::RngMPolElt
                                              : Terms := 1) -> [ ]
{ Expands the Puiseux expansion s of a reduced bivariate polynomial }
require Rank(Parent(f)) eq 2: "Argument f must be a bivariate polynomial";
  return [x];
end intrinsic

NewtonPuiseuxAlgorithmReducedLoop := function(f, ord, terms)
  Q<x> := PuiseuxSeriesRing(CoefficientRing(Parent(f)));
  x0 := Parent(f).1; y0 := Parent(f).2;
  // Step (i.a): Select only those factors containing the 0 branch.
  S := yFactor(f) gt 0 select [<Q!0, Parent(f)!0>] else [];
  // Step (i.b): For each side...
  for F in Faces(NewtonPolygon(f)) do
    n := GradientVector(F)[1]; m := GradientVector(F)[2];
    // Apply the change of variables (1).
    C := Reverse(Coefficients(n eq 1 select f else Evaluate(f, 1, x0^n), 2));
    // For each root...
    for a in [<Root(a[1], n), a[2]> : a in Roots(FacePolynomial(F))] do
      // Apply the change of variables (2) & get the sub-solution recursively.
      R := (a[2] ne 1 and terms lt -1) or terms gt 0 select
        NewtonPuiseuxAlgorithmReducedLoop([i gt 1 select
          C[i] + Self(i-1)*x0^m*(a[1] + y0) else C[1] : i in [1..#C]][#C],
        ord*n, terms - 1) else [<Q!0, f>];
      // Undo the change of variables.
      S cat:= [<x^(m/n)*(a[1] + ChangePrecision(Composition(s[1], x^(1/n)),
        Infinity())), s[2]> : s in R];
    end for;
  end for; return S;
end function;
